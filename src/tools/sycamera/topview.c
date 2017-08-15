/* Topview output */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "domain.h"
#include "sycamera.h"
#include "sycout.h"

#ifdef USE_MPI
#	include <mpi.h>
#	include "smpi.h"
#endif

camera_image *sycout_topview_result, *sycout_topview_camim;
char *sycout_topview_filename;
enum sycout_topview_type sycout_topview_brightness=1;
double sycout_topview_majorradius=0.;

#pragma omp threadprivate(sycout_topview_camim)

void sycout_topview_init(struct general_settings *settings) {
    int i, j;
    sycout_topview_result = malloc(sizeof(camera_image));
    sycout_topview_result->pixels = 0;

    /* Load settings */
    for (i = 0; i < settings->n; i++) {
        if (!strcmp(settings->setting[i], "brightness")) {
            if (!strcmp(settings->value[i], "intensity"))
                sycout_topview_brightness = SYCOUT_TOPVIEW_TYPE_INTENSITY;
            else if (!strcmp(settings->value[i], "bw"))
                sycout_topview_brightness = SYCOUT_TOPVIEW_TYPE_BW;
			else if (!strcmp(settings->value[i], "histogram"))
				sycout_topview_brightness = SYCOUT_TOPVIEW_TYPE_HIST;
            else {
                fprintf(stderr, "WARNING: Unrecognized choice for sycout topview option 'brightness': %s\n", settings->value[i]);
            }
        } else if (!strcmp(settings->setting[i], "name")) {
            sycout_topview_filename = settings->value[i];
        } else if (!strcmp(settings->setting[i], "pixels")) {
			sycout_topview_result->pixels = atoi(settings->value[i]);
        } else {
            fprintf(stderr, "WARNING: Unrecognized sycout topview setting '%s'\n", settings->setting[i]);
        }
    }

    if (sycout_topview_result->pixels <= 0) {
		fprintf(stderr, "ERROR: Invalid number of pixels selected (%d)!\n", sycout_topview_result->pixels);
		exit(-1);
    }

    /* Initialize resulting image */
	sycout_topview_result->canvas = malloc(sizeof(double*)*sycout_topview_result->pixels);
	for (i = 0; i < sycout_topview_result->pixels; i++) {
		sycout_topview_result->canvas[i] = malloc(sizeof(double)*sycout_topview_result->pixels);
		for (j = 0; j < sycout_topview_result->pixels; j++) {
			sycout_topview_result->canvas[i][j] = 0;
		}
	}

    sycout_topview_majorradius = domain_get_major_radius();
}
void sycout_topview_init_run(void) {
    int i, j;
	sycout_topview_camim = malloc(sizeof(camera_image));
	sycout_topview_camim->pixels = sycout_topview_result->pixels;

	sycout_topview_camim->canvas = malloc(sizeof(double*)*sycout_topview_camim->pixels);
	for (i = 0; i < sycout_topview_camim->pixels; i++) {
		sycout_topview_camim->canvas[i] = malloc(sizeof(double)*sycout_topview_camim->pixels);
		for (j = 0; j < sycout_topview_camim->pixels; j++) {
			sycout_topview_camim->canvas[i][j] = 0;
		}
	}
}
void sycout_topview_init_particle(particle *p) {}
void sycout_topview_deinit_run(void) {
	#pragma omp critical
	{
		int i, j;
		for (i = 0; i < sycout_topview_camim->pixels; i++) {
			for (j = 0; j < sycout_topview_camim->pixels; j++) {
                if (sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_INTENSITY ||
					sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_HIST)
    				sycout_topview_result->canvas[i][j] += sycout_topview_camim->canvas[i][j];
                else if (sycout_topview_camim->canvas[i][j] > .5)
                    sycout_topview_result->canvas[i][j] = 1;
			}
		}
	}
}
void sycout_topview_step(struct sycout_data *data) {
    int i = (int)round(data->sd->x/sycout_topview_majorradius * sycout_topview_camim->pixels/2.0)+sycout_topview_camim->pixels/2;
    int j = (int)round(data->sd->y/sycout_topview_majorradius * sycout_topview_camim->pixels/2.0)+sycout_topview_camim->pixels/2;

    if (sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_INTENSITY)
        sycout_topview_camim->canvas[j][i] += data->brightness * data->differential;
	else if (sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_HIST)
		sycout_topview_result->canvas[j][i] += 1;
    else
        sycout_topview_camim->canvas[j][i] = 1.;
}

/*
void sycout_topview_combine(FILE *f, camera_image *ci) {
    int i, j;
    double tempval = 0.;
    for (i = 0; i < ci->pixels; i++) {
        for (j = 0; j < ci->pixels; j++) {
            if (fscanf(f, "%le", &tempval) != 1) {
                fprintf(stderr, "Error reading back image during sycout topview output!\n");
                exit(-1);
            }

            if (sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_INTENSITY ||
				sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_HIST)
                ci->canvas[i][j] += tempval;
            else if (tempval >= .5)
                ci->canvas[i][j] = 1.;
        }
    }
}
void sycout_topview_output(FILE *f, camera_image *ci) {
	int i, j;
	for (i = 0; i < ci->pixels; i++) {
		fprintf(f, "%.15e", ci->canvas[0][i]);
		for (j = 1; j < ci->pixels; j++) {
			fprintf(f, "  %.15e", ci->canvas[i][j]);
		}
		fprintf(f, "\n");
	}
}
*/
void sycout_topview_combine(sFILE *sf, camera_image *ci) {
	int i, j;
	sfilesize_t dims[2];
	double **img = sf->get_doubles(sf, "image", dims);
	
    for (i = 0; i < ci->pixels; i++) {
        for (j = 0; j < ci->pixels; j++) {
			if (sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_INTENSITY ||
				sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_HIST)
				ci->canvas[i][j] += img[i][j];
			else if (sycout_topview_brightness == SYCOUT_TOPVIEW_TYPE_BW &&
				img[i][j] >= .5)
				ci->canvas[i][j] = 1;
		}
	}
}
void sycout_topview_output(sFILE *sf, camera_image *ci) {
	double *wall[2];
	domain *d = domain_get();
	wall[0] = d->r;
	wall[1] = d->z;

	sf->write_list(sf, "detectorPosition", Rdet->val, 3);
	sf->write_list(sf, "detectorDirection", ddet->val, 3);
	sf->write_list(sf, "detectorVisang", &visang, 1);
	sf->write_image(sf, "image", ci->canvas, ci->pixels);
	sf->write_array(sf, "wall", wall, 2, d->n);
}
void sycout_topview_write(int mpi_rank, int nprocesses) {
	//FILE *f;
	sFILE *sf;
	enum sfile_type ftype;

#ifdef USE_MPI
	printf("[%d] (sycout topview) Waiting for 'output ready' from previous process.\n", mpi_rank);
	smpi_wor(SYCOUT_MPIID_TOPVIEW);
	printf("[%d] (sycout topview) Received 'output ready' signal from previous process.\n", mpi_rank);
#endif

	ftype = sfile_get_filetype(sycout_topview_filename);
	if (ftype == FILETYPE_UNKNOWN) {
		ftype = FILETYPE_SDT;
		if (mpi_rank == 0)
			printf("[%d] (sycout image) WARNING: Unable to determine filetype of output. Defaulting to SDT.\n", mpi_rank);
	}

	sf = sfile_init(ftype);

	/* First, if this is not the root process
	 * (which always writes first), read in
	 * the image that has already been written
	 * and add it to our image.
	 */
	if (mpi_rank > 0) {
		if (sf->open(sf, sycout_topview_filename, SFILE_MODE_READ)) {
			sycout_topview_combine(sf, sycout_topview_result);
			sf->close(sf);
		} else {
    		fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'\n", mpi_rank, sycout_topview_filename);
		}
	}

	if (!sf->open(sf, sycout_topview_filename, SFILE_MODE_WRITE)) {
        fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'\n", mpi_rank, sycout_topview_filename);
        exit(1);
	}

	sycout_topview_output(sf, sycout_topview_result);
	sf->close(sf);
/*
    if (mpi_rank > 0) {
    	f = fopen(sycout_topview_filename, "r");
    	if (f != NULL) {
            sycout_topview_combine(f, sycout_topview_result);
            fclose(f);
    	} else {
    		fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'\n", mpi_rank, sycout_topview_filename);
    		perror("WARNING");
        }
    }

    f = fopen(sycout_topview_filename, "w");
    if (f == NULL) {
        fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'\n", mpi_rank, sycout_topview_filename);
        perror("ERROR");
        exit(1);
    }

    sycout_topview_output(f, sycout_topview_result);
    fclose(f);
*/

#ifdef USE_MPI
    printf("[%d] (sycout topview) Done, sending 'output ready' to next process.\n", mpi_rank);
    smpi_sor(SYCOUT_MPIID_IMAGE);

	if (mpi_rank == nprocesses-1) {
#endif
		printf("Wrote camera topview to '%s'\n", sycout_topview_filename);
#ifdef USE_MPI
	}
#endif
}
