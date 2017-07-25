/* Compute Greens function Ihat */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "config.h"
#include "particles.h"
#include "sfile.h"
#include "sycamera.h"
#include "sycout.h"

#ifdef USE_MPI
#	include "smpi.h"
#endif

long long int sycout_green_nvel1, sycout_green_nvel2,
	sycout_green_nrad, sycout_green_nwav,
	sycout_green_pixels,
	sycout_green_ivel1, sycout_green_ivel2,
	sycout_green_irad,
	sycout_green_cvel2, sycout_green_cvel1,
	sycout_green_crad, sycout_green_cwav,
	sycout_green_subpixels,
	sycout_green_suboffseti, sycout_green_suboffsetj;

size_t sycout_green_func_sz;

#pragma omp threadprivate(sycout_green_irad,sycout_green_ivel1,sycout_green_ivel2)

char *sycout_green_output;
enum sycout_green_type sycout_green_tp;
enum sfile_type sycout_green_outformat;
double sycout_green_dlambda;
double *sycout_green_func;

void sycout_green_init(struct general_settings *settings) {
	sycout_green_tp = SYCOUT_GREEN_IMAGE;
	sycout_green_outformat = FILETYPE_SDT;
	sycout_green_output = NULL;
	sycout_green_pixels = 0;
	sycout_green_func = NULL;
	sycout_green_func_sz = 0;

	sycout_green_suboffseti = 0;
	sycout_green_suboffsetj = 0;
	sycout_green_subpixels  = 0;

	int i, outformat_set=0;
	for (i = 0; i < settings->n; i++) {
		if (!strcmp(settings->setting[i], "format")) {
			sycout_green_outformat = sfile_name2filetype(settings->value[i]);	
			outformat_set = 1;
		} else if (!strcmp(settings->setting[i], "output")) {
			sycout_green_output = settings->value[i];
		} else if (!strcmp(settings->setting[i], "pixels")) {
			sycout_green_pixels  = atoi(settings->value[i]);
		} else if (!strcmp(settings->setting[i], "suboffseti")) {
			sycout_green_suboffseti = atoi(settings->value[i]);
		} else if (!strcmp(settings->setting[i], "suboffsetj")) {
			sycout_green_suboffsetj = atoi(settings->value[i]);
		} else if (!strcmp(settings->setting[i], "subpixels")) {
			sycout_green_subpixels = atoi(settings->value[i]);
		} else if (!strcmp(settings->setting[i], "function")) {
			if (!strcmp(settings->value[i], "full")) {
				sycout_green_tp = SYCOUT_GREEN_FULL;
			} else if (!strcmp(settings->value[i], "image")) {
				sycout_green_tp = SYCOUT_GREEN_IMAGE;
			} else if (!strcmp(settings->value[i], "spectrum")) {
				sycout_green_tp = SYCOUT_GREEN_SPECTRUM;
			} else if (!strcmp(settings->value[i], "total")) {
				sycout_green_tp = SYCOUT_GREEN_TOTAL;
			} else {
				fprintf(stderr, "ERROR: sycout green: Invalid choice for option 'function': %s.\n", settings->value[i]);
				exit(-1);
			}
		} else {
			fprintf(stderr, "ERROR: Unrecognized option '%s'.\n", settings->setting[i]);
			exit(-1);
		}
	}

	/* Output set? */
	if (sycout_green_output == NULL) {
		fprintf(stderr, "ERROR: sycout green: No output filename specified.\n");
		exit(-1);
	}

	/* Green's function type */
	if (sycout_green_tp == SYCOUT_GREEN_IMAGE || sycout_green_tp == SYCOUT_GREEN_FULL) {
		if (sycout_green_pixels <= 0) {
			fprintf(stderr, "ERROR: The number of pixels for the Green's function has not been set.\n");
			exit(-1);
		}

		/* Verify pixel bounds */
		if (sycout_green_subpixels == 0)
			sycout_green_subpixels = sycout_green_pixels;
		
		if (sycout_green_suboffseti+sycout_green_subpixels > sycout_green_pixels) {
			fprintf(stderr, "ERROR: Subset image large than in i direction than actual image.\n");
			exit(-1);
		}
		if (sycout_green_suboffsetj+sycout_green_subpixels > sycout_green_pixels) {
			fprintf(stderr, "ERROR: Subset image large than in j direction than actual image.\n");
			exit(-1);
		}
	} else if (sycout_green_tp == SYCOUT_GREEN_SPECTRUM) {
	} else if (sycout_green_tp == SYCOUT_GREEN_TOTAL) {}

	/* Was the output format set? */
	if (outformat_set) {
		if (sycout_green_outformat == FILETYPE_UNKNOWN) {
			fprintf(stderr, "ERROR: sycout green: Unrecognized output file format specified.\n");
			exit(-1);
		}
	} else {
		sycout_green_outformat = sfile_get_filetype(sycout_green_output);
		if (sycout_green_outformat == FILETYPE_UNKNOWN) {
			fprintf(stderr, "ERROR: sycout green: Unrecognized output file format of file '%s'.\n", sycout_green_output);
			exit(-1);
		}
	}
}
void sycout_green_init_run(void) {
	int ijk[3];
	size_t i;

	#pragma omp critical
	{
		if (sycout_green_func == NULL) {
			particles_gridsize3(ijk);
			sycout_green_nrad  = ijk[0];
			sycout_green_nvel1 = ijk[1];
			sycout_green_nvel2 = ijk[2];

			sycout_green_nwav  = sycamera_spectrum_length();

			size_t size;
			/* Green's function type */
			if (sycout_green_tp == SYCOUT_GREEN_IMAGE) {
				size = sycout_green_subpixels*sycout_green_subpixels*sycout_green_nvel1*sycout_green_nvel2*sycout_green_nrad;

				sycout_green_cvel2 = sycout_green_subpixels*sycout_green_subpixels;
				sycout_green_cvel1 = sycout_green_nvel2 * sycout_green_cvel2;
				sycout_green_crad  = sycout_green_cvel1*sycout_green_nvel1;
			} else if (sycout_green_tp == SYCOUT_GREEN_SPECTRUM) {
				size = sycout_green_nvel1*sycout_green_nvel2*sycout_green_nrad*sycout_green_nwav;

				sycout_green_cvel2 = sycout_green_nwav;
				sycout_green_cvel1 = sycout_green_nvel2 * sycout_green_cvel2;
				sycout_green_crad  = sycout_green_cvel1*sycout_green_nvel1;
			} else if (sycout_green_tp == SYCOUT_GREEN_TOTAL) {
				size = sycout_green_nvel1*sycout_green_nvel2*sycout_green_nrad;

				sycout_green_cvel2 = 1;
				sycout_green_cvel1 = sycout_green_nvel2;
				sycout_green_crad  = sycout_green_cvel1*sycout_green_nvel1;
			} else {/* Both */
				size = sycout_green_subpixels*sycout_green_subpixels*sycout_green_nvel1*sycout_green_nvel2*sycout_green_nrad*sycout_green_nwav;

				sycout_green_cvel2 = sycout_green_subpixels*sycout_green_subpixels*sycout_green_nwav;
				sycout_green_cvel1 = sycout_green_nvel2 * sycout_green_cvel2;
				sycout_green_crad  = sycout_green_cvel1*sycout_green_nvel1;
			}
			
			/* Print note about how much memory is required */
			char suffix[] = " kMGTP";
			double rsize = size*sizeof(double);
			i = 0;
			while (rsize > 1000.0) {rsize /= 1000.0; i++;}

#ifdef USE_MPI
			printf("Green's function requires 2 x %.2f %cB\n", rsize, suffix[i]);
#else
			printf("Green's function requires %.2f %cB\n", rsize, suffix[i]);
#endif

			sycout_green_func_sz = size;
			sycout_green_func = malloc(sizeof(double)*size);
			
			/* Initialize to zero */
			//memset(sycout_green_func, 0, size*sizeof(double));
			for (i = 0; i < size; i++) {
				sycout_green_func[i] = 0.0;
			}
		}
	}
}
void sycout_green_init_particle(particle *p) {
	/* Set indices */
	sycout_green_irad  = p->ir;
	sycout_green_ivel1 = p->iv1;
	sycout_green_ivel2 = p->iv2;
}

void sycout_green_deinit_run(void) {}
void sycout_green_step(struct sycout_data *data) {
	//long long int i = (long long int)(data->i*sycout_green_subpixels),
	//	j = (long long int)(data->j*sycout_green_subpixels),
	long long int i = (long long int)(data->i*sycout_green_pixels),
		j = (long long int)(data->j*sycout_green_pixels),
		index;
	
	if (sycout_green_tp == SYCOUT_GREEN_IMAGE) {	/* IMAGE */
		/* Ignore pixels outside of subimage */
		if (i < sycout_green_suboffseti || i >= sycout_green_subpixels+sycout_green_suboffseti)
			return;
		if (j < sycout_green_suboffsetj || j >= sycout_green_subpixels+sycout_green_suboffsetj)
			return;

		i -= sycout_green_suboffseti;
		j -= sycout_green_suboffsetj;
			
		/*
		index = i * sycout_green_cpixels2 +
				j * sycout_green_cpixels  +
				sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2;
		index = index*sycout_green_nrad + sycout_green_irad;
		*/
		index = sycout_green_irad * sycout_green_crad +
				sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2 * sycout_green_cvel2 +
				i * sycout_green_subpixels + j;

		sycout_green_func[index] += data->brightness * data->RdPhi * data->Jdtdrho / particles_get_drho();
	} else if (sycout_green_tp == SYCOUT_GREEN_SPECTRUM) {	/* SPECTRUM */
		/*
		index = sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2;
		index = index*sycout_green_nrad + sycout_green_irad;
		index *= sycout_green_nwav;
		*/
		index = sycout_green_irad * sycout_green_crad +
				sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2 * sycout_green_cvel2;

		double *spectrum = sycamera_spectrum_get();
		double diffel = data->RdPhi * data->Jdtdrho;
		long long int  i;
		for (i = 0; i < sycout_green_nwav; i++) {
			sycout_green_func[index+i] += spectrum[i] * diffel;
		}
	} else if (sycout_green_tp == SYCOUT_GREEN_FULL) {	/* FULL */
		/* Ignore pixels outside of subimage */
		if (i < sycout_green_suboffseti || i >= sycout_green_subpixels+sycout_green_suboffseti)
			return;
		if (j < sycout_green_suboffsetj || j >= sycout_green_subpixels+sycout_green_suboffsetj)
			return;

		i -= sycout_green_suboffseti;
		j -= sycout_green_suboffsetj;
			
		/*
		index = i * sycout_green_cpixels2 +
				j * sycout_green_cpixels  +
				sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2;
		index = index*sycout_green_nrad + sycout_green_irad;
		index *= sycout_green_nwav;
		*/
		index = sycout_green_irad * sycout_green_crad +
				sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2 * sycout_green_cvel2 +
				i * sycout_green_subpixels + j;

		double *spectrum = sycamera_spectrum_get();
		double diffel = data->RdPhi * data->Jdtdrho / particles_get_drho();
		long long int si, gi, pixels2 = sycout_green_subpixels*sycout_green_subpixels;
		for (si = gi = 0; si < sycout_green_nwav; si++, gi += pixels2) {
			sycout_green_func[index+gi] += spectrum[si] * diffel;
		}
	} else if (sycout_green_tp == SYCOUT_GREEN_TOTAL) {	/* TOTAL */
		/*
		index = sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2;
		index = index*sycout_green_nrad + sycout_green_irad;
		*/
		index = sycout_green_irad * sycout_green_crad +
				sycout_green_ivel1 * sycout_green_cvel1 +
				sycout_green_ivel2;

		double diffel = data->RdPhi * data->Jdtdrho / particles_get_drho();
		sycout_green_func[index] += data->brightness * diffel;
	}
}
#ifdef USE_MPI
void sycout_green_mpi_receive(int nprocesses) {
	int i;
	size_t j;
	double *buffer = malloc(sycout_green_func_sz*sizeof(double));
	for (i = 0; i < nprocesses-1; i++) {
		fprintf(stderr, "[0]: Waiting for Green's function from MPI process %d.\n", i);
		smpi_receive_matrix(buffer, sycout_green_func_sz, i, SYCOUT_MPIID_GREEN);

		/* Add matrices */
		for (j = 0; j < sycout_green_func_sz; j++) {
			sycout_green_func[j] += buffer[j];
		}

		fprintf(stderr, "[0]: Received and assimilated matrix from MPI process %d.\n", i);
	}
}
#endif/*USE_MPI*/
void sycout_green_write(int mpi_rank, int nprocesses) {
	int i;
	double *b, *v, pixels;
	char *s;
#ifdef USE_MPI
	if (mpi_rank == nprocesses-1) {
		sycout_green_mpi_receive(nprocesses);
#endif
		sFILE *sf = sfile_init(sycout_green_outformat);
		sf->open(sf, sycout_green_output, SFILE_MODE_WRITE);

		/* Type of function */
		switch (sycout_green_tp) {
			case SYCOUT_GREEN_IMAGE:
				sf->write_string(sf, "type", "image", 5);
				break;
			case SYCOUT_GREEN_SPECTRUM:
				sf->write_string(sf, "type", "spectrum", 8);
				break;
			case SYCOUT_GREEN_TOTAL:
				sf->write_string(sf, "type", "total", 5);
				break;
			case SYCOUT_GREEN_FULL:
				sf->write_string(sf, "type", "both", 4);
				break;
			default:
				sf->write_string(sf, "type", "unknown", 7);
				break;
		}

		/* Number of pixels */
		if (sycout_green_tp == SYCOUT_GREEN_SPECTRUM || sycout_green_tp == SYCOUT_GREEN_TOTAL) {
			pixels = 1;
			sf->write_list(sf, "pixels", &pixels, 1);
		} else {
			pixels = sycout_green_subpixels;
			sf->write_list(sf, "pixels", &pixels, 1);
		}

		/* Type of param1 */
		s = particles_param1_name();
		sf->write_string(sf, "param1name", s, strlen(s));

		/* Type of param2 */
		s = particles_param2_name();
		sf->write_string(sf, "param2name", s, strlen(s));

		/* r vector */
		b = particles_get_bounds();
		v = malloc(sizeof(double)*b[2]);

		for (i = 0; i < b[2]; i++)
			v[i] = b[0] + ((double)i)/((double)b[2]) * (b[1]-b[0]);

		sf->write_list(sf, "r", v, b[2]);
		free(v);

		/* param1 vector */
		v = malloc(sizeof(double)*b[5]);

		for (i = 0; i < b[5]; i++)
			v[i] = b[3] + ((double)i)/((double)b[5]) * (b[4]-b[3]);

		sf->write_list(sf, "param1", v, b[5]);
		free(v);

		/* param2 vector */
		v = malloc(sizeof(double)*b[8]);

		for (i = 0; i < b[8]; i++)
			v[i] = b[6] + ((double)i)/((double)b[8]) * (b[7]-b[6]);

		sf->write_list(sf, "param2", v, b[8]);
		free(v);

		/* Wavelength vector */
		if (sycout_green_tp == SYCOUT_GREEN_IMAGE || sycout_green_tp == SYCOUT_GREEN_TOTAL) {
			double wavelengths = 0;
			sf->write_list(sf, "wavelengths", &wavelengths, 1);
		} else {
			v = sycamera_spectrum_get_wavelengths();
			sf->write_list(sf, "wavelengths", v, sycout_green_nwav);
		}

		/* Actual Green's function */
		sf->write_list(sf, "func", sycout_green_func, sycout_green_func_sz);

		sf->close(sf);

		fprintf(stdout, "[%d]: Wrote Green's function to '%s'.\n", mpi_rank, sycout_green_output);
#ifdef USE_MPI
	} else {
		/* Send Green's function to root process over MPI */
		fprintf(stdout, "[%d]: Sending Green's function to process %d.\n", mpi_rank, nprocesses-1);
		smpi_send_matrix(sycout_green_func, sycout_green_func_sz, nprocesses-1, SYCOUT_MPIID_GREEN);
		fprintf(stdout, "[%d]: Done sending Green's function to process %d.\n", mpi_rank, nprocesses-1);
	}
#endif
}

