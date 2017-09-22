/* Spectrometer handler */

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

char *sycout_spectrometer_filename;
double *sycout_spectrometer_result, *sycout_spectrometer_wavelengths,
	   *sycout_spectrometer_lwavelengths, *sycout_spectrometer_lresult;
int sycout_spectrometer_nlambdas;
long long int sycout_spectrometer_lcounts, sycout_spectrometer_counts;

#pragma omp threadprivate(sycout_spectrometer_filename,sycout_spectrometer_lwavelengths, \
	sycout_spectrometer_nlambdas,sycout_spectrometer_lresult,sycout_spectrometer_lcounts)

void sycout_spectrometer_init(struct general_settings *settings) {
	int i;
	sycout_spectrometer_result = NULL;
	sycout_spectrometer_wavelengths = NULL;
	sycout_spectrometer_counts = 0;

	/* Load settings */
	for (i = 0; i < settings->n; i++) {
		if (!strcmp(settings->setting[i], "name")) {
			sycout_spectrometer_filename = settings->value[i];
		} else {
			fprintf(stderr, "WARNING: Unrecognized spectrometer setting '%s'.\n", settings->setting[i]);
		}
	}
}

void sycout_spectrometer_init_run(void) {
	int i;
	double *lambdas = sycamera_get_wavelengths();
	sycout_spectrometer_nlambdas = sycamera_get_spectrum_length();
	sycout_spectrometer_lcounts = 0;

	sycout_spectrometer_lwavelengths = malloc(sizeof(double)*sycout_spectrometer_nlambdas);
	sycout_spectrometer_lresult = malloc(sizeof(double)*sycout_spectrometer_nlambdas);

	for (i = 0; i < sycout_spectrometer_nlambdas; i++) {
		sycout_spectrometer_lwavelengths[i] = lambdas[i];
		sycout_spectrometer_lresult[i] = 0.0;
	}
}

void sycout_spectrometer_output(FILE *f, double *wavelengths, double *spectrum, int n);

void sycout_spectrometer_init_particle(particle *p) {}
void sycout_spectrometer_init_step(void) {}
void sycout_spectrometer_deinit_run(void) {
	#pragma omp critical
	{
		int i;
		if (sycout_spectrometer_result == NULL) {
			sycout_spectrometer_result = malloc(sizeof(double)*sycout_spectrometer_nlambdas);
			sycout_spectrometer_wavelengths = malloc(sizeof(double)*sycout_spectrometer_nlambdas);

			for (i = 0; i < sycout_spectrometer_nlambdas; i++) {
				sycout_spectrometer_result[i] = 0.0;
				sycout_spectrometer_wavelengths[i] = sycout_spectrometer_lwavelengths[i];
			}
		}

		for (i = 0; i < sycout_spectrometer_nlambdas; i++) {
			sycout_spectrometer_result[i] += sycout_spectrometer_lresult[i];
		}
		sycout_spectrometer_counts += sycout_spectrometer_lcounts;
	}
}

/**
 * Accumulate spectrum
 */
void sycout_spectrometer_step(struct sycout_data *data) {
	int i;
	double *p = sycamera_get_spectrum();
	for (i = 0; i < sycout_spectrometer_nlambdas; i++) {
		sycout_spectrometer_lresult[i] += p[i] * data->differential;
	}

	sycout_spectrometer_lcounts++;
}

void sycout_spectrometer_combine(FILE *f, double *spectrum, int n, int mpi_rank, int nprocesses) {
	int i;
	double tval1, tval2;
	if (fscanf(f, "%le", &tval1) != 1) {
		fprintf(stderr, "Error reading back spectrum count during sycout spectrometer output.\n");
		exit(-1);
	}

	sycout_spectrometer_counts += tval1;

	for (i = 0; i < n; i++) {
		if (fscanf(f, "%le,%le", &tval1, &tval2) != 2) {
			fprintf(stderr, "Error reading back spectrum during sycout spectrometer output.\n");
			exit(-1);
		}

		spectrum[i] += tval2;
	}
}
void sycout_spectrometer_output(FILE *f, double *wavelengths, double *spectrum, int n) {
	int i;
	for (i = 0; i < n; i++) {
		fprintf(f, "%.12e,%.12e\n", wavelengths[i], spectrum[i]);
	}
}
void sycout_spectrometer_normalize(double *spectrum, int n) {
	int i;
	for (i = 0; i < n; i++) {
		if (sycout_spectrometer_counts > 0)
			spectrum[i] /= (double)sycout_spectrometer_counts;
	}
}

void sycout_spectrometer_write(int mpi_rank, int nprocesses) {
	FILE *f;

#ifdef USE_MPI
	printf("[%d] (sycout spectrometer) Waiting for 'output ready' from previous process.\n", mpi_rank);
	smpi_wor(SYCOUT_MPIID_SPECTROMETER);
	printf("[%d] (sycout spectrometer) Received 'output ready' signal from previous process.\n", mpi_rank);
#endif

	if (sycout_spectrometer_result == NULL) return;

	if (mpi_rank > 0) {
		f = fopen(sycout_spectrometer_filename, "r");
		if (f != NULL) {
			sycout_spectrometer_combine(f, sycout_spectrometer_result, sycout_spectrometer_nlambdas, mpi_rank, nprocesses);
			fclose(f);
		} else {
			fprintf(stderr, "[%d] WARNING: Unable to open file for reading: '%s'.\n", mpi_rank, sycout_spectrometer_filename);
			perror("WARNING");
		}
	}

	f = fopen(sycout_spectrometer_filename, "w");
	if (f == NULL) {
		fprintf(stderr, "[%d] ERROR: Unable to open file for writing: '%s'.\n", mpi_rank, sycout_spectrometer_filename);
		perror("ERROR");
		exit(1);
	}

#ifdef USE_MPI
	if (mpi_rank == nprocesses-1) {
#endif
	sycout_spectrometer_normalize(sycout_spectrometer_result, sycout_spectrometer_nlambdas);
#ifdef USE_MPI
	}
#endif

	sycout_spectrometer_output(f, sycout_spectrometer_wavelengths, sycout_spectrometer_result, sycout_spectrometer_nlambdas);
	fclose(f);

#ifdef USE_MPI
    printf("[%d] (sycout spectrometer) Done, sending 'output ready' to next process.\n", mpi_rank);
    smpi_sor(SYCOUT_MPIID_IMAGE);

	if (mpi_rank == nprocesses-1) {
#endif
		printf("Wrote spectrum to '%s'\n", sycout_spectrometer_filename);
#ifdef USE_MPI
	}
#endif
}

