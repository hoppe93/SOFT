/* Polarization spectrometer handler */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <sycamera.h>
#include "config.h"
#include "sycamera.h"

#ifdef USE_MPI
#	include <mpi.h>
#	include "smpi.h"
#endif

char *sycout_polspectrometer_filename;
double **sycout_polspectrometer_result, *sycout_polspectrometer_wavelengths,
	   *sycout_polspectrometer_lwavelengths, **sycout_polspectrometer_lresult;
int sycout_polspectrometer_nlambdas;
long long int sycout_polspectrometer_lcounts, sycout_polspectrometer_counts;

#pragma omp threadprivate(sycout_polspectrometer_filename,sycout_polspectrometer_lwavelengths, \
	sycout_polspectrometer_nlambdas,sycout_polspectrometer_lresult,sycout_polspectrometer_lcounts)

void sycout_polspectrometer_init(struct general_settings *settings) {
	int i;
	sycout_polspectrometer_result = NULL;
	sycout_polspectrometer_wavelengths = NULL;
	sycout_polspectrometer_counts = 0;

	/* Load settings */
	for (i = 0; i < settings->n; i++) {
		if (!strcmp(settings->setting[i], "name")) {
			sycout_polspectrometer_filename = settings->value[i];
		} else {
			fprintf(stderr, "WARNING: Unrecognized polspectrometer setting '%s'.\n", settings->setting[i]);
		}
	}
}

void sycout_polspectrometer_init_run(void) {
	int i;
	double *lambdas = sycamera_get_wavelengths();
	sycout_polspectrometer_nlambdas = sycamera_get_spectrum_length();
	sycout_polspectrometer_lcounts = 0;

	sycout_polspectrometer_lwavelengths = malloc(sizeof(double)*sycout_polspectrometer_nlambdas);
	sycout_polspectrometer_lresult = malloc(sizeof(double*)*4);
	sycout_polspectrometer_lresult[0] = malloc(sizeof(double)*sycout_polspectrometer_nlambdas*4);
	sycout_polspectrometer_lresult[1] = sycout_polspectrometer_lresult[0] + sycout_polspectrometer_nlambdas;
	sycout_polspectrometer_lresult[2] = sycout_polspectrometer_lresult[1] + sycout_polspectrometer_nlambdas;
	sycout_polspectrometer_lresult[3] = sycout_polspectrometer_lresult[2] + sycout_polspectrometer_nlambdas;

	for (i = 0; i < sycout_polspectrometer_nlambdas; i++) {
		sycout_polspectrometer_lwavelengths[i] = lambdas[i];
		sycout_polspectrometer_lresult[0][i] = 0.0;
		sycout_polspectrometer_lresult[1][i] = 0.0;
		sycout_polspectrometer_lresult[2][i] = 0.0;
		sycout_polspectrometer_lresult[3][i] = 0.0;
	}
}

void sycout_polspectrometer_init_particle(particle *p) {}
void sycout_polspectrometer_init_step(void) {}
void sycout_polspectrometer_deinit_run(void) {
	#pragma omp critical
	{
		int i;
		if (sycout_polspectrometer_result == NULL) {
			sycout_polspectrometer_result = malloc(sizeof(double*)*4);
			sycout_polspectrometer_result[0] = malloc(sizeof(double)*sycout_polspectrometer_nlambdas*4);
			sycout_polspectrometer_result[1] = sycout_polspectrometer_result[0] + sycout_polspectrometer_nlambdas;
			sycout_polspectrometer_result[2] = sycout_polspectrometer_result[1] + sycout_polspectrometer_nlambdas;
			sycout_polspectrometer_result[3] = sycout_polspectrometer_result[2] + sycout_polspectrometer_nlambdas;
			sycout_polspectrometer_wavelengths = malloc(sizeof(double)*sycout_polspectrometer_nlambdas);

			for (i = 0; i < sycout_polspectrometer_nlambdas; i++) {
				sycout_polspectrometer_result[0][i] = 0.0;
				sycout_polspectrometer_result[1][i] = 0.0;
				sycout_polspectrometer_result[2][i] = 0.0;
				sycout_polspectrometer_result[3][i] = 0.0;
				sycout_polspectrometer_wavelengths[i] = sycout_polspectrometer_lwavelengths[i];
			}
		}

		for (i = 0; i < sycout_polspectrometer_nlambdas; i++) {
			sycout_polspectrometer_result[0][i] += sycout_polspectrometer_lresult[0][i];
			sycout_polspectrometer_result[1][i] += sycout_polspectrometer_lresult[1][i];
			sycout_polspectrometer_result[2][i] += sycout_polspectrometer_lresult[2][i];
			sycout_polspectrometer_result[3][i] += sycout_polspectrometer_lresult[3][i];
		}
		sycout_polspectrometer_counts += sycout_polspectrometer_lcounts;
	}
}

/**
 * Accumulate spectrum
 */
void sycout_polspectrometer_step(struct sycout_data *data) {
	int i;
	double **p = sycamera_get_polarization_spectrum();
	for (i = 0; i < sycout_polspectrometer_nlambdas; i++) {
		sycout_polspectrometer_lresult[0][i] += p[0][i] * data->differential;
		sycout_polspectrometer_lresult[1][i] += p[1][i] * data->differential;
		sycout_polspectrometer_lresult[2][i] += p[2][i] * data->differential;
		sycout_polspectrometer_lresult[3][i] += p[3][i] * data->differential;
	}

	sycout_polspectrometer_lcounts++;
}

void sycout_polspectrometer_output(FILE *f, double *wavelengths, double **spectrum, int n) {
	int i;
	fprintf(f, "#wavelength,Alr2,Aud2,ARe,AIm\n");
	for (i = 0; i < n; i++) {
		fprintf(f, "%.12e,%.12e,%.12e,%.12e,%.12e\n", wavelengths[i], spectrum[0][i], spectrum[1][i], spectrum[2][i], spectrum[3][i]);
	}
}

/*
void sycout_polspectrometer_normalize(double **spectrum, int n) {
	int i;
	for (i = 0; i < n; i++) {
		if (sycout_polspectrometer_counts > 0) {
			spectrum[0][i] /= (double)sycout_polspectrometer_counts;
			spectrum[1][i] /= (double)sycout_polspectrometer_counts;
			spectrum[2][i] /= (double)sycout_polspectrometer_counts;
			spectrum[3][i] /= (double)sycout_polspectrometer_counts;
		}
	}
}
*/

/**
 * "High-level" interface for generating output.
 * This function assembles spectra from all threads
 * and makes the root process write an output file.
 */
void sycout_polspectrometer_write(int mpi_rank, int nprocesses) {
	FILE *f;

	if (sycout_polspectrometer_result == NULL) return;

#ifdef USE_MPI
	if (mpi_rank == 0) {
		/* Receive spectra from other processes */
		int i, j;
		double *tmp = malloc(sizeof(double)*sycout_polspectrometer_nlambdas*4);
		for (i = 1; i < nprocesses; i++) {
			smpi_receive_matrix(tmp, sycout_polspectrometer_nlambdas, i, SYCOUT_MPIID_POLSPECTROMETER);
			for (j = 0; j < sycout_polspectrometer_nlambdas*4; j++)
				sycout_polspectrometer_result[0][j] += tmp[j];
		}

		free(tmp);

#endif/*USE_MPI*/
		f = fopen(sycout_polspectrometer_filename, "w");
		if (f == NULL) {
			fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'.\n", mpi_rank, sycout_polspectrometer_filename);
			perror("ERROR");
			exit(1);
		}

	/*
	#ifdef USE_MPI
		if (mpi_rank == nprocesses-1) {
	#endif
		sycout_polspectrometer_normalize(sycout_polspectrometer_result, sycout_polspectrometer_nlambdas);
	#ifdef USE_MPI
		}
	#endif
	*/

		sycout_polspectrometer_output(f, sycout_polspectrometer_wavelengths, sycout_polspectrometer_result, sycout_polspectrometer_nlambdas);
		fclose(f);

#ifdef USE_MPI
	} else {
		smpi_send_matrix(sycout_polspectrometer_result[0], sycout_polspectrometer_nlambdas*4, 0, SYCOUT_MPIID_POLSPECTROMETER);
	}
#endif/*USE_MPI*/
/*
	if (mpi_rank > 0) {
		f = fopen(sycout_polspectrometer_filename, "r");
		if (f != NULL) {
			sycout_polspectrometer_combine(f, sycout_polspectrometer_result, sycout_polspectrometer_nlambdas, mpi_rank, nprocesses);
			fclose(f);
		} else {
			fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'.\n", mpi_rank, sycout_polspectrometer_filename);
			perror("WARNING");
		}
	}

*/
	printf("Wrote spectrum to '%s'\n", sycout_polspectrometer_filename);
}

