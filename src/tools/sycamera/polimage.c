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
	sycout_polimage_result->StokesI = malloc(sizeof(double*)*sycout_polimage_result->pixels);
	sycout_polimage_result->StokesQ = malloc(sizeof(double*)*sycout_polimage_result->pixels);
	sycout_polimage_result->StokesU = malloc(sizeof(double*)*sycout_polimage_result->pixels);
	sycout_polimage_result->StokesV = malloc(sizeof(double*)*sycout_polimage_result->pixels);

	sycout_polimage_result->StokesI[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);
	sycout_polimage_result->StokesQ[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);
	sycout_polimage_result->StokesU[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);
	sycout_polimage_result->StokesV[0] = malloc(sizeof(double)*sycout_polimage_result->pixels*sycout_polimage_result->pixels);

	for (i = 0; i < sycout_polimage_result->pixels; i++) {
		if (i > 0) {
			sycout_polimage_result->StokesI[i] = sycout_polimage_result->StokesI[i-1] + sycout_polimage_result->pixels;
			sycout_polimage_result->StokesQ[i] = sycout_polimage_result->StokesQ[i-1] + sycout_polimage_result->pixels;
			sycout_polimage_result->StokesU[i] = sycout_polimage_result->StokesU[i-1] + sycout_polimage_result->pixels;
			sycout_polimage_result->StokesV[i] = sycout_polimage_result->StokesV[i-1] + sycout_polimage_result->pixels;
		}

		for (j = 0; j < sycout_polimage_result->pixels; j++) {
			sycout_polimage_result->StokesI[i][j] = 0;
			sycout_polimage_result->StokesQ[i][j] = 0;
			sycout_polimage_result->StokesU[i][j] = 0;
			sycout_polimage_result->StokesV[i][j] = 0;
		}
	}
}
void sycout_polimage_init_run(void) {
	/* Initialize local image */
    int i, j;
	sycout_polimage_camim = malloc(sizeof(camera_polimage));
	sycout_polimage_camim->pixels = sycout_polimage_result->pixels;
	sycout_polimage_camim->StokesI = malloc(sizeof(double*)*sycout_polimage_camim->pixels);
	sycout_polimage_camim->StokesQ = malloc(sizeof(double*)*sycout_polimage_camim->pixels);
	sycout_polimage_camim->StokesU = malloc(sizeof(double*)*sycout_polimage_camim->pixels);
	sycout_polimage_camim->StokesV = malloc(sizeof(double*)*sycout_polimage_camim->pixels);

	sycout_polimage_camim->StokesI[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);
	sycout_polimage_camim->StokesQ[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);
	sycout_polimage_camim->StokesU[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);
	sycout_polimage_camim->StokesV[0] = malloc(sizeof(double)*sycout_polimage_camim->pixels*sycout_polimage_camim->pixels);

	for (i = 0; i < sycout_polimage_camim->pixels; i++) {
		if (i > 0) {
			sycout_polimage_camim->StokesI[i] = sycout_polimage_camim->StokesI[i-1] + sycout_polimage_camim->pixels;
			sycout_polimage_camim->StokesQ[i] = sycout_polimage_camim->StokesQ[i-1] + sycout_polimage_camim->pixels;
			sycout_polimage_camim->StokesU[i] = sycout_polimage_camim->StokesU[i-1] + sycout_polimage_camim->pixels;
			sycout_polimage_camim->StokesV[i] = sycout_polimage_camim->StokesV[i-1] + sycout_polimage_camim->pixels;
		}

		for (j = 0; j < sycout_polimage_camim->pixels; j++)
			sycout_polimage_camim->StokesI[i][j] = 0.0;
		for (j = 0; j < sycout_polimage_camim->pixels; j++)
			sycout_polimage_camim->StokesQ[i][j] = 0.0;
		for (j = 0; j < sycout_polimage_camim->pixels; j++)
			sycout_polimage_camim->StokesU[i][j] = 0.0;
		for (j = 0; j < sycout_polimage_camim->pixels; j++)
			sycout_polimage_camim->StokesV[i][j] = 0.0;
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
				sycout_polimage_result->StokesI[i][j] += sycout_polimage_camim->StokesI[i][j];
				sycout_polimage_result->StokesQ[i][j] += sycout_polimage_camim->StokesQ[i][j];
				sycout_polimage_result->StokesU[i][j] += sycout_polimage_camim->StokesU[i][j];
				sycout_polimage_result->StokesV[i][j] += sycout_polimage_camim->StokesV[i][j];
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
	double *pol = sycamera_get_polarization(),
			sd = data->differential;
	
	if (pol == NULL) return;

	sycout_polimage_camim->StokesI[j][i] += pol[0] * sd;
	sycout_polimage_camim->StokesQ[j][i] += pol[1] * sd;
	sycout_polimage_camim->StokesU[j][i] += pol[2] * sd;
	sycout_polimage_camim->StokesV[j][i] += pol[3] * sd;
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
	sf->write_image(sf, "StokesI", ci->StokesI, ci->pixels);
	sf->write_image(sf, "StokesQ", ci->StokesQ, ci->pixels);
	sf->write_image(sf, "StokesU", ci->StokesU, ci->pixels);
	sf->write_image(sf, "StokesV", ci->StokesV, ci->pixels);
	sf->write_array(sf, "wall", wall, 2, d->n);
}
/**
 * "High-level" interface to outputting an image.
 * Assembles the image from all MPI processes and
 * makes the root process write the output file in
 * an appropriate format.
 */
void sycout_polimage_write(int mpi_rank, int nprocesses) {
	sFILE *sf;
	enum sfile_type ftype;

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
	/*
	if (mpi_rank > 0) {
		if (sf->open(sf, sycout_polimage_filename, SFILE_MODE_READ)) {
			sycout_polimage_combine(sf, sycout_polimage_result);
			sf->close(sf);
		} else {
    		fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'\n", mpi_rank, sycout_polimage_filename);
		}
	}
	*/

#ifdef USE_MPI
	/**
	 * First, if root process, gather images
	 * from other processes.
	 */
	int pixels2 = sycout_polimage_camim->pixels*sycout_polimage_camim->pixels;
	if (mpi_rank == 0) {
		int i, j;
		double *tmp = malloc(sizeof(double)*pixels2);
		for (i = 1; i < nprocesses; i++) {
			/* I */
			smpi_receive_matrix(tmp, pixels2, i, SYCOUT_MPIID_POLIMAGE);
			for (j = 0; j < pixels2; j++)
				sycout_polimage_result->StokesI[0][j] += tmp[j];

			/* Q */
			smpi_receive_matrix(tmp, pixels2, i, SYCOUT_MPIID_POLIMAGE);
			for (j = 0; j < pixels2; j++)
				sycout_polimage_result->StokesQ[0][j] += tmp[j];

			/* U */
			smpi_receive_matrix(tmp, pixels2, i, SYCOUT_MPIID_POLIMAGE);
			for (j = 0; j < pixels2; j++)
				sycout_polimage_result->StokesU[0][j] += tmp[j];

			/* V */
			smpi_receive_matrix(tmp, pixels2, i, SYCOUT_MPIID_POLIMAGE);
			for (j = 0; j < pixels2; j++)
				sycout_polimage_result->StokesV[0][j] += tmp[j];
		}
		free(tmp);

		/* Then write the output file */
#endif/*USE_MPI*/

		if (!sf->open(sf, sycout_polimage_filename, SFILE_MODE_WRITE)) {
			fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'\n", mpi_rank, sycout_polimage_filename);
			exit(1);
		}

		sycout_polimage_output(sf, sycout_polimage_result);
		sf->close(sf);

#ifdef USE_MPI
	} else {	/* Else, send data to root process */
		smpi_send_matrix(sycout_polimage_camim->StokesI[0], sycout_polimage_camim->pixels*sycout_polimage_camim->pixels, 0, SYCOUT_MPIID_POLIMAGE);
		smpi_send_matrix(sycout_polimage_camim->StokesQ[0], sycout_polimage_camim->pixels*sycout_polimage_camim->pixels, 0, SYCOUT_MPIID_POLIMAGE);
		smpi_send_matrix(sycout_polimage_camim->StokesU[0], sycout_polimage_camim->pixels*sycout_polimage_camim->pixels, 0, SYCOUT_MPIID_POLIMAGE);
		smpi_send_matrix(sycout_polimage_camim->StokesV[0], sycout_polimage_camim->pixels*sycout_polimage_camim->pixels, 0, SYCOUT_MPIID_POLIMAGE);
	}
#endif

	printf("Wrote camera polimage to '%s'\n", sycout_polimage_filename);
}

