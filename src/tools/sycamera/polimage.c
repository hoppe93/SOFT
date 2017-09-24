/* Camera polarized image handler */

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

camera_polimage *sycout_polimage_result, *sycout_polimage_camim;
char *sycout_polimage_filename;

#pragma omp threadprivate(sycout_polimage_camim)

void sycout_polimage_init(struct general_settings *settings) {
    int i, j;
    sycout_polimage_result = malloc(sizeof(camera_polimage));
    sycout_polimage_result->pixels = 0;

    /* Load settings */
    for (i = 0; i < settings->n; i++) {
        if (!strcmp(settings->setting[i], "pixels")) {
			sycout_polimage_result->pixels = atoi(settings->value[i]);
        } else if (!strcmp(settings->setting[i], "name")) {
            sycout_polimage_filename = settings->value[i];
        } else {
            fprintf(stderr, "WARNING: Unrecognized sycout polimage setting '%s'\n", settings->setting[i]);
        }
    }

    if (sycout_polimage_result->pixels <= 0) {
		fprintf(stderr, "ERROR: Invalid number of pixels selected (%d)!\n", sycout_polimage_result->pixels);
		exit(-1);
    }

    /* Initialize resulting image */
	sycout_polimage_result->AlrRe = malloc(sizeof(double*)*sycout_polimage_result->pixels);
	sycout_polimage_result->AlrIm = malloc(sizeof(double*)*sycout_polimage_result->pixels);
	sycout_polimage_result->AudRe = malloc(sizeof(double*)*sycout_polimage_result->pixels);
	sycout_polimage_result->AudIm = malloc(sizeof(double*)*sycout_polimage_result->pixels);

	sycout_polimage_result->AlrRe[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);
	sycout_polimage_result->AlrIm[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);
	sycout_polimage_result->AudRe[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);
	sycout_polimage_result->AudIm[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);

	for (i = 0; i < sycout_polimage_result->pixels; i++) {
		if (i > 0) {
			sycout_polimage_result->AlrRe[i] = sycout_polimage_result->AlrRe[i-1] + sycout_polimage_result->pixels;
			sycout_polimage_result->AlrIm[i] = sycout_polimage_result->AlrIm[i-1] + sycout_polimage_result->pixels;
			sycout_polimage_result->AudRe[i] = sycout_polimage_result->AudRe[i-1] + sycout_polimage_result->pixels;
			sycout_polimage_result->AudIm[i] = sycout_polimage_result->AudIm[i-1] + sycout_polimage_result->pixels;
		}

		for (j = 0; j < sycout_polimage_result->pixels; j++) {
			sycout_polimage_result->AlrRe[i][j] = 0;
			sycout_polimage_result->AlrIm[i][j] = 0;
			sycout_polimage_result->AudRe[i][j] = 0;
			sycout_polimage_result->AudIm[i][j] = 0;
		}
	}
}
void sycout_polimage_init_run(void) {
    int i, j;
	sycout_polimage_camim = malloc(sizeof(camera_polimage));
	sycout_polimage_camim->pixels = sycout_polimage_result->pixels;
	sycout_polimage_camim->AlrRe = malloc(sizeof(double*)*sycout_polimage_camim->pixels);
	sycout_polimage_camim->AlrIm = malloc(sizeof(double*)*sycout_polimage_camim->pixels);
	sycout_polimage_camim->AudRe = malloc(sizeof(double*)*sycout_polimage_camim->pixels);
	sycout_polimage_camim->AudIm = malloc(sizeof(double*)*sycout_polimage_camim->pixels);

	sycout_polimage_camim->AlrRe[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);
	sycout_polimage_camim->AlrIm[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);
	sycout_polimage_camim->AudRe[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);
	sycout_polimage_camim->AudIm[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);

	for (i = 0; i < sycout_polimage_camim->pixels; i++) {
		if (i > 0) {
			sycout_polimage_camim->AlrRe[i] = sycout_polimage_camim->AlrRe[i-1] + sycout_polimage_camim->pixels;
			sycout_polimage_camim->AlrIm[i] = sycout_polimage_camim->AlrIm[i-1] + sycout_polimage_camim->pixels;
			sycout_polimage_camim->AudRe[i] = sycout_polimage_camim->AudRe[i-1] + sycout_polimage_camim->pixels;
			sycout_polimage_camim->AudIm[i] = sycout_polimage_camim->AudIm[i-1] + sycout_polimage_camim->pixels;
		}

		for (j = 0; j < sycout_polimage_camim->pixels; j++) {
			sycout_polimage_camim->AlrRe[i][j] = 0;
			sycout_polimage_camim->AlrIm[i][j] = 0;
			sycout_polimage_camim->AudRe[i][j] = 0;
			sycout_polimage_camim->AudIm[i][j] = 0;
		}
	}
}
void sycout_polimage_init_particle(particle *p) {}
void sycout_polimage_init_step(void) {}
void sycout_polimage_deinit_run(void) {
	#pragma omp critical
	{
		int i, j;
		for (i = 0; i < sycout_polimage_camim->pixels; i++) {
			for (j = 0; j < sycout_polimage_camim->pixels; j++) {
				sycout_polimage_result->AlrRe[i][j] += sycout_polimage_camim->AlrRe[i][j];
				sycout_polimage_result->AlrIm[i][j] += sycout_polimage_camim->AlrIm[i][j];
				sycout_polimage_result->AudRe[i][j] += sycout_polimage_camim->AudRe[i][j];
				sycout_polimage_result->AudIm[i][j] += sycout_polimage_camim->AudIm[i][j];
			}
		}
	}
}
/**
 * Gather data from one step of SOFT.
 */
void sycout_polimage_step(struct sycout_data *data) {
    int i = (int)(data->i*sycout_polimage_result->pixels);
    int j = (int)(data->j*sycout_polimage_result->pixels);
	double *polarization = sycamera_get_polarization(),
			sd = sqrt(data->differential);
	
	if (polarization == NULL) return;

	sycout_polimage_camim->AlrRe[j][i] += polarization[0] * sd;
	sycout_polimage_camim->AlrIm[j][i] += polarization[1] * sd;
	sycout_polimage_camim->AudRe[j][i] += polarization[2] * sd;
	sycout_polimage_camim->AudIm[j][i] += polarization[3] * sd;
}

/**
 * Combine the polimage in a file with the
 * given polimage (basically, add them together)
 */
void sycout_polimage_combine(sFILE *sf, camera_polimage *ci) {
	int i, j;
	sfilesize_t dims[2];
	double **AlrRe = sf->get_doubles(sf, "AlrRe", dims),
		   **AlrIm = sf->get_doubles(sf, "AlrIm", dims),
		   **AudRe = sf->get_doubles(sf, "AudRe", dims),
		   **AudIm = sf->get_doubles(sf, "AudIm", dims);
	
    for (i = 0; i < ci->pixels; i++) {
        for (j = 0; j < ci->pixels; j++) {
			ci->AlrRe[i][j] += AlrRe[i][j];
			ci->AlrIm[i][j] += AlrIm[i][j];
			ci->AudRe[i][j] += AudRe[i][j];
			ci->AudIm[i][j] += AudIm[i][j];
		}
	}
}

/**
 * Write contents of polarized image to file.
 */
void sycout_polimage_output(sFILE *sf, camera_polimage *ci) {
	double *wall[2];
	domain *d = domain_get();
	wall[0] = d->r;
	wall[1] = d->z;

	sf->write_list(sf, "detectorPosition", Rdet->val, 3);
	sf->write_list(sf, "detectorDirection", ddet->val, 3);
	sf->write_list(sf, "detectorVisang", &visang, 1);
	sf->write_image(sf, "AlrRe", ci->AlrRe, ci->pixels);
	sf->write_image(sf, "AlrIm", ci->AlrIm, ci->pixels);
	sf->write_image(sf, "AudRe", ci->AudRe, ci->pixels);
	sf->write_image(sf, "AudIm", ci->AudIm, ci->pixels);
	sf->write_array(sf, "wall", wall, 2, d->n);
}
void sycout_polimage_write(int mpi_rank, int nprocesses) {
	sFILE *sf;
	enum sfile_type ftype;

#ifdef USE_MPI
	printf("[%d] (sycout polimage) Waiting for 'output ready' from previous process.\n", mpi_rank);
	smpi_wor(SYCOUT_MPIID_POLIMAGE);
	printf("[%d] (sycout polimage) Received 'output ready' signal from previous process.\n", mpi_rank);
#endif

	ftype = sfile_get_filetype(sycout_polimage_filename);
	if (ftype == FILETYPE_UNKNOWN) {
		ftype = FILETYPE_SDT;
		if (mpi_rank == 0)
			printf("[%d] (sycout polimage) WARNING: Unable to determine filetype of output. Defaulting to SDT.\n", mpi_rank);
	}

	sf = sfile_init(ftype);

	/* First, if this is not the root process
	 * (which always writes first), read in
	 * the image that has already been written
	 * and add it to our image.
	 */
	if (mpi_rank > 0) {
		if (sf->open(sf, sycout_polimage_filename, SFILE_MODE_READ)) {
			sycout_polimage_combine(sf, sycout_polimage_result);
			sf->close(sf);
		} else {
    		fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'\n", mpi_rank, sycout_polimage_filename);
		}
	}

	if (!sf->open(sf, sycout_polimage_filename, SFILE_MODE_WRITE)) {
        fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'\n", mpi_rank, sycout_polimage_filename);
        exit(1);
	}

	sycout_polimage_output(sf, sycout_polimage_result);
	sf->close(sf);

#ifdef USE_MPI
    printf("[%d] (sycout polimage) Done, sending 'output ready' to next process.\n", mpi_rank);
    smpi_sor(SYCOUT_MPIID_POLIMAGE);

	if (mpi_rank == nprocesses-1) {
#endif
		printf("Wrote camera polimage to '%s'\n", sycout_polimage_filename);
#ifdef USE_MPI
	}
#endif
}
