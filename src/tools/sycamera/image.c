/* Camera image handler */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "config.h"
#include "domain.h"
#include "sfile.h"
#include "sycamera.h"
#include "sycout.h"

#ifdef USE_MPI
#   include <mpi.h>
#	include "smpi.h"
#endif

camera_image *sycout_image_result, *sycout_image_camim;
char *sycout_image_filename;
enum sycout_image_type sycout_image_brightness=SYCOUT_IMAGE_TYPE_INTENSITY;

#pragma omp threadprivate(sycout_image_camim)

void sycout_image_init(struct general_settings *settings) {
    int i, j;
    sycout_image_result = malloc(sizeof(camera_image));
    sycout_image_result->pixels = 0;

    /* Load settings */
    for (i = 0; i < settings->n; i++) {
		if (!strcmp(settings->setting[i], "brightness")) {
			if (!strcmp(settings->value[i], "intensity"))
				sycout_image_brightness = SYCOUT_IMAGE_TYPE_INTENSITY;
			else if (!strcmp(settings->value[i], "bw"))
				sycout_image_brightness = SYCOUT_IMAGE_TYPE_BW;
			else if (!strcmp(settings->value[i], "histogram"))
				sycout_image_brightness = SYCOUT_IMAGE_TYPE_HIST;
			else {
				fprintf(stderr, "WARNING: Unrecognized choice for sycout image option 'brightness': %s\n", settings->value[i]);
			}
        } else if (!strcmp(settings->setting[i], "pixels")) {
			sycout_image_result->pixels = atoi(settings->value[i]);
        } else if (!strcmp(settings->setting[i], "name")) {
            sycout_image_filename = settings->value[i];
        } else {
            fprintf(stderr, "WARNING: Unrecognized sycout image setting '%s'\n", settings->setting[i]);
        }
    }

    if (sycout_image_result->pixels <= 0) {
		fprintf(stderr, "ERROR: Invalid number of pixels selected (%d)!\n", sycout_image_result->pixels);
		exit(-1);
    }

    /* Initialize resulting image */
	sycout_image_result->canvas = malloc(sizeof(double*)*sycout_image_result->pixels);
	for (i = 0; i < sycout_image_result->pixels; i++) {
		sycout_image_result->canvas[i] = malloc(sizeof(double)*sycout_image_result->pixels);
		for (j = 0; j < sycout_image_result->pixels; j++) {
			sycout_image_result->canvas[i][j] = 0;
		}
	}
}
void sycout_image_init_run(void) {
    int i, j;
	sycout_image_camim = malloc(sizeof(camera_image));
	sycout_image_camim->pixels = sycout_image_result->pixels;

	sycout_image_camim->canvas = malloc(sizeof(double*)*sycout_image_camim->pixels);
	for (i = 0; i < sycout_image_camim->pixels; i++) {
		sycout_image_camim->canvas[i] = malloc(sizeof(double)*sycout_image_camim->pixels);
		for (j = 0; j < sycout_image_camim->pixels; j++) {
			sycout_image_camim->canvas[i][j] = 0;
		}
	}
}
void sycout_image_init_particle(particle *p) {}
void sycout_image_init_step(void) {}
void sycout_image_deinit_run(void) {
	#pragma omp critical
	{
		int i, j;
		for (i = 0; i < sycout_image_camim->pixels; i++) {
			for (j = 0; j < sycout_image_camim->pixels; j++) {
				if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_INTENSITY ||
					sycout_image_brightness == SYCOUT_IMAGE_TYPE_HIST)
					sycout_image_result->canvas[i][j] += sycout_image_camim->canvas[i][j];
				else if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_BW) {
					if (sycout_image_camim->canvas[i][j] > 0.5)
						sycout_image_result->canvas[i][j] = 1;
				}
			}
		}
	}
}
void sycout_image_step(struct sycout_data *data) {
    int i = (int)(data->i*sycout_image_result->pixels);
    int j = (int)(data->j*sycout_image_result->pixels);

	//if (i >= sycout_image_result->pixels) i = sycout_image_result->pixels-1;
	//if (j >= sycout_image_result->pixels) j = sycout_image_result->pixels-1;
	//if (j < 586 || j > 599) return;
	//if (j < 490 || j > 510) return;

	if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_INTENSITY)
		sycout_image_camim->canvas[j][i] += data->brightness * data->differential;
	else if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_HIST)
		sycout_image_camim->canvas[j][i] += data->differential;
	else if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_BW)
		sycout_image_camim->canvas[j][i] = 1;
}

/*
void sycout_image_interp_pixels_algorithm(
	double l0x, double l0y, double l1x, double l1y,
	double img_side, double **sollist, int *nsol
) {
	double ax = sycamera_lastlx - l0x,
		   ay = sycamera_lastly - l0y,
		   n0, nlast, t0, tlast,
		   h = img_side / 
	
	/ * Determine t0 and tlast * /
	if (l0x < l1x) {
		if (l0x > 0) t0 = l0x;
		else t0 = 0;

		if (l1x < img_side) tlast = l1x;
		else tlast = img_side;
	} else {
		if (l1x > 0) t0 = l1x;
		else t0 = 0;

		if (l0x < img_side) tlast = l0x;
		else tlast = img_side;
	}

	/ * Compute n0 and nlast * /
	if (ax != 0) {
		double f = ay/ax;
		if (f > 0) {
			n0 = ceil()
		}
	} else {
	
	}
}
*/

/**
 * Interpolate between pixels in order to even out
 * the registered emission between pixels.
 *
 * l1x: x-component of l vector at new particle position
 * l1y: y-component of l vector at new particle position
 * di: Difference in the i pixel index between last and
 *     current particle position.
 * dj: Difference in the j pixel index between last and
 *     current particle position.
 * img_side: Sidelength of image
 * sd: Step-data object
 * intensity: Intensity reaching the camera
 * dt: Length of timestep
 * dPhi: Length of toroidal angle step
 */
/*
void sycout_image_interp_pixels(
	double l1x, double l1y, int di, int dj, double img_side,
	step_data *sd, double intensity, double dt, double dPhi
) {
	double solutions[2][di+dj];/ * x and y components of solutions * /
	int nsol = 0;

	sycamera_interp_pixels_algorithm(l0x, l0y, l1x, l1y, solutions, &nsol);
	sycamera_interp_pidels_algorithm(l0y, l0x, l1y, l1x, solutions, &nsol);

	sycamera_register_radiation(i, j, sd, intensity, dt, dPhi);
}
*/

/*
void sycout_image_combine(FILE *f, camera_image *ci) {
    int i, j;
    double tempval = 0.;
    for (i = 0; i < ci->pixels; i++) {
        for (j = 0; j < ci->pixels; j++) {
            if (fscanf(f, "%le", &tempval) != 1) {
                fprintf(stderr, "Error reading back image during sycout image output!\n");
                exit(-1);
            }

			if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_INTENSITY ||
				sycout_image_brightness == SYCOUT_IMAGE_TYPE_HIST)
				ci->canvas[i][j] += tempval;
			else if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_BW &&
				tempval >= .5)
				ci->canvas[i][j] = 1;
        }
    }
}
void sycout_image_output(FILE *f, camera_image *ci) {
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
void sycout_image_combine(sFILE *sf, camera_image *ci) {
	int i, j;
	sfilesize_t dims[2];
	double **img = sf->get_doubles(sf, "image", dims);
	
    for (i = 0; i < ci->pixels; i++) {
        for (j = 0; j < ci->pixels; j++) {
			if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_INTENSITY ||
				sycout_image_brightness == SYCOUT_IMAGE_TYPE_HIST)
				ci->canvas[i][j] += img[i][j];
			else if (sycout_image_brightness == SYCOUT_IMAGE_TYPE_BW &&
				img[i][j] >= .5)
				ci->canvas[i][j] = 1;
		}
	}
}
void sycout_image_output(sFILE *sf, camera_image *ci) {
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
void sycout_image_write(int mpi_rank, int nprocesses) {
	//FILE *f;
	sFILE *sf;
	enum sfile_type ftype;

#ifdef USE_MPI
	printf("[%d] (sycout image) Waiting for 'output ready' from previous process.\n", mpi_rank);
	smpi_wor(SYCOUT_MPIID_IMAGE);
	printf("[%d] (sycout image) Received 'output ready' signal from previous process.\n", mpi_rank);
#endif

	ftype = sfile_get_filetype(sycout_image_filename);
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
		if (sf->open(sf, sycout_image_filename, SFILE_MODE_READ)) {
			sycout_image_combine(sf, sycout_image_result);
			sf->close(sf);
		} else {
    		fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'\n", mpi_rank, sycout_image_filename);
		}
	}

	if (!sf->open(sf, sycout_image_filename, SFILE_MODE_WRITE)) {
        fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'\n", mpi_rank, sycout_image_filename);
        exit(1);
	}

	sycout_image_output(sf, sycout_image_result);
	sf->close(sf);

/*
    if (mpi_rank > 0) {
    	f = fopen(sycout_image_filename, "r");
    	if (f != NULL) {
            sycout_image_combine(f, sycout_image_result);
            fclose(f);
    	} else {
    		fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'\n", mpi_rank, sycout_image_filename);
    		perror("WARNING");
        }
    }

    f = fopen(sycout_image_filename, "w");
    if (f == NULL) {
        fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'\n", mpi_rank, sycout_image_filename);
        perror("ERROR");
        exit(1);
    }

    sycout_image_output(f, sycout_image_result);
    fclose(f);
*/

#ifdef USE_MPI
    printf("[%d] (sycout image) Done, sending 'output ready' to next process.\n", mpi_rank);
    smpi_sor(SYCOUT_MPIID_POLSPECTROMETER);

	if (mpi_rank == nprocesses-1) {
#endif
		printf("Wrote camera image to '%s'\n", sycout_image_filename);
#ifdef USE_MPI
	}
#endif
}
