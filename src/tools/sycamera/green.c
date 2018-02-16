/* Compute Greens function Ihat */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "config.h"
#include "particles.h"
#include "sfile.h"
#include "sycamera.h"
#include "sycout.h"
#include "util.h"

#ifdef USE_MPI
#	include "smpi.h"
#endif

size_t
	sycout_green_pixels,
	sycout_green_nvel1, sycout_green_nvel2,
	sycout_green_nrad, sycout_green_nwav,
	sycout_green_ivel1, sycout_green_ivel2,
	sycout_green_irad,
	sycout_green_subpixels,
	sycout_green_suboffseti, sycout_green_suboffsetj;

int sycout_green_weighWdf, sycout_green_haswav,
	sycout_green_hasrho, sycout_green_hasvel1, sycout_green_hasvel2,
    sycout_green_stokesparams;

size_t sycout_green_func_sz;
#define SYCOUT_GREEN_MAXDIMS 6
enum sycout_green_dimension sycout_green_format[SYCOUT_GREEN_MAXDIMS];
size_t sycout_green_factors[SYCOUT_GREEN_MAXDIMS];

#pragma omp threadprivate(sycout_green_irad,sycout_green_ivel1,sycout_green_ivel2)

char *sycout_green_output;
enum sfile_type sycout_green_outformat;
double sycout_green_dlambda;
double *sycout_green_func;      /* Green's function */

void sycout_green_init(struct general_settings *settings) {
	sycout_green_outformat = FILETYPE_SDT;
	sycout_green_output = NULL;
	sycout_green_pixels = 0;
	sycout_green_func = NULL;
	sycout_green_func_sz = 0;
	sycout_green_weighWdf = 0;
	sycout_green_haswav = 0;
	sycout_green_hasrho = 0;
	sycout_green_hasvel1 = 0;
	sycout_green_hasvel2 = 0;
    sycout_green_stokesparams = 0;

	sycout_green_suboffseti = 0;
	sycout_green_suboffsetj = 0;
	sycout_green_subpixels  = 0;

	int i,j,k, outformat_set=0;
	for (i = 0; i < SYCOUT_GREEN_MAXDIMS; i++)
		sycout_green_format[i] = SYCOUT_GREEN_NONE;

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
        } else if (!strcmp(settings->setting[i], "stokesparams")) {
            if (!strcmp(settings->value[i], "yes"))
                sycout_green_stokesparams = 1;
            else
                sycout_green_stokesparams = 0;
		} else if (!strcmp(settings->setting[i], "function")) {
			enum sycout_green_dimension dim;
			j = 0;
			while (settings->value[i][j]) {
				dim = SYCOUT_GREEN_NONE;
				switch (settings->value[i][j]) {
					case '1': dim = SYCOUT_GREEN_VEL1; sycout_green_weighWdf++; sycout_green_hasvel1 = 1; break;
					case '2': dim = SYCOUT_GREEN_VEL2; sycout_green_weighWdf++; sycout_green_hasvel2 = 1; break;
					case 'i': dim = SYCOUT_GREEN_IMAGEI; break;
					case 'j': dim = SYCOUT_GREEN_IMAGEJ; break;
					case 'r': dim = SYCOUT_GREEN_RADIUS; sycout_green_weighWdf++; sycout_green_hasrho = 1; break;
					case 'w': dim = SYCOUT_GREEN_SPECTRUM; sycout_green_haswav = 1; break;
					default:
						fprintf(stderr, "ERROR: Unrecognized Green's function dimension: %c. Allowed values: 12irw.\n", settings->value[i][j]);
						exit(EXIT_FAILURE);
				}

				for (k = 0; k < SYCOUT_GREEN_MAXDIMS; k++) {
					if (sycout_green_format[k] == SYCOUT_GREEN_NONE) {
						sycout_green_format[k] = dim;
						break;
					} else if (sycout_green_format[k] == dim) {
						fprintf(stderr, "ERROR: Green's function dimensions may not be repeated in the format.\n");
						exit(EXIT_FAILURE);
					}
				}

				j++;
			}

			/* If all of r, v1 and v2 are free parameters (i.e.
			 * dimensions in the GF) we shouldn't weigh with
			 * any distribution function. */
			if (sycout_green_weighWdf == 3) sycout_green_weighWdf = 0;
			else sycout_green_weighWdf = 1;
		} else {
			fprintf(stderr, "ERROR: Unrecognized option '%s'.\n", settings->setting[i]);
			exit(-1);
		}
	}

	/* Output set? */
	if (sycout_green_output == NULL) {
		fprintf(stderr, "ERROR: (sycout green): No output filename specified.\n");
		exit(EXIT_FAILURE);
	}

	/* Check Green's function format */
	if (sycout_green_format[0] == SYCOUT_GREEN_NONE) {
		fprintf(stderr, "ERROR: (sycout green): The Green's function is 0-dimensional.\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < SYCOUT_GREEN_MAXDIMS && sycout_green_format[i] != SYCOUT_GREEN_NONE; i++) {
		if (sycout_green_format[i] == SYCOUT_GREEN_IMAGEI) {
			if (sycout_green_pixels == 0) {
				fprintf(stderr, "ERROR: (sycout green): Invalid number of pixels set: %zu.\n", sycout_green_pixels);
				exit(EXIT_FAILURE);
			}

			if (sycout_green_subpixels == 0)
				sycout_green_subpixels = sycout_green_pixels;

			if (sycout_green_suboffseti+sycout_green_subpixels > sycout_green_pixels) {
				fprintf(stderr, "ERROR: (sycout green): Subset image larger in i direction than actual image.\n");
				exit(EXIT_FAILURE);
			}
		} else if (sycout_green_format[i] == SYCOUT_GREEN_IMAGEJ) {
			if (sycout_green_pixels == 0) {
				fprintf(stderr, "ERROR: (sycout green): Invalid number of pixels set: %zu.\n", sycout_green_pixels);
				exit(EXIT_FAILURE);
			}

			if (sycout_green_subpixels == 0)
				sycout_green_subpixels = sycout_green_pixels;
			
			if (sycout_green_suboffsetj+sycout_green_subpixels > sycout_green_pixels) {
				fprintf(stderr, "ERROR: (sycout green): Subset image larger in j direction than actual image.\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	/* Was the output format set? */
	if (outformat_set) {
		if (sycout_green_outformat == FILETYPE_UNKNOWN) {
			fprintf(stderr, "ERROR: (sycout green): Unrecognized output file format specified.\n");
			exit(-1);
		}
	} else {
		sycout_green_outformat = sfile_get_filetype(sycout_green_output);
		if (sycout_green_outformat == FILETYPE_UNKNOWN) {
			fprintf(stderr, "ERROR: (sycout green): Unrecognized output file format of file '%s'.\n", sycout_green_output);
			exit(-1);
		}
	}
}
void sycout_green_init_run(void) {
	int ijk[3], j;
	size_t i;

	#pragma omp critical
	{
		if (sycout_green_func == NULL) {
			particles_gridsize3(ijk);
			sycout_green_nrad  = ijk[0];
			sycout_green_nvel1 = ijk[1];
			sycout_green_nvel2 = ijk[2];
			sycout_green_nwav  = sycamera_get_spectrum_length();

			printf("Green's function format: ");
			char *s, S[20];
			for (i = 0; i < SYCOUT_GREEN_MAXDIMS && sycout_green_format[i] != SYCOUT_GREEN_NONE; i++) {
				if (i > 0) printf(" x ");
				switch (sycout_green_format[i]) {
					case SYCOUT_GREEN_RADIUS: printf("RADIUS"); break;
					case SYCOUT_GREEN_SPECTRUM: printf("SPECTRUM"); break;
					case SYCOUT_GREEN_IMAGEI: printf("PIXELS-I"); break;
					case SYCOUT_GREEN_IMAGEJ: printf("PIXELS-J"); break;
					case SYCOUT_GREEN_VEL1:
						s = particles_param1_name();
						strcpy(S, s);
						for (j=0; S[j]; j++) S[j] = chrupr(S[j]);
						printf("%s", S);
						break;
					case SYCOUT_GREEN_VEL2:
						s = particles_param2_name();
						strcpy(S, s);
						for (j=0; S[j]; j++) S[j] = chrupr(S[j]);
						printf("%s", S);
						break;
					default: break;
				}
			}

            if (sycout_green_stokesparams)
                printf(" x STOKES-PARAMETERS");

			putc('\n', stdout);

			/* Compute size of Green's function and build list of factors */
			size_t size = 1;
			for (i = 0; i < SYCOUT_GREEN_MAXDIMS; i++) {
				sycout_green_factors[i] = 1;
            }

			for (j = SYCOUT_GREEN_MAXDIMS-1; j >= 0; j--) {
				if (sycout_green_format[j] == SYCOUT_GREEN_NONE) continue;

				switch (sycout_green_format[j]) {
					case SYCOUT_GREEN_VEL1:
						size *= sycout_green_nvel1;
						if (j-1 >= 0) sycout_green_factors[j-1] = sycout_green_nvel1;
						break;
					case SYCOUT_GREEN_VEL2:
						size *= sycout_green_nvel2;
						if (j-1 >= 0) sycout_green_factors[j-1] = sycout_green_nvel2;
						break;
					case SYCOUT_GREEN_IMAGEI:
					case SYCOUT_GREEN_IMAGEJ:
						size *= sycout_green_subpixels;
						if (j-1 >= 0) sycout_green_factors[j-1] = sycout_green_subpixels;
						break;
					case SYCOUT_GREEN_RADIUS:
						size *= sycout_green_nrad;
						if (j-1 >= 0) sycout_green_factors[j-1] = sycout_green_nrad;
						break;
					case SYCOUT_GREEN_SPECTRUM:
						size *= sycout_green_nwav;
						if (j-1 >= 0) sycout_green_factors[j-1] = sycout_green_nwav;
						break;
					default:
						fprintf(
							stderr,
							"WARNING: (sycout green): Unrecognized Green's function dimension: %d. Ignoring when allocating memory.\n",
							sycout_green_format[j]
						);
						break;
				}
			}

            /* If we are to store Stokes parameters,
             * rather than just plain intensities, we
             * multiply the size by the number of Stokes
             * parameters (4, IQUV) */
            if (sycout_green_stokesparams)
                size *= 4;

			/* Generate factors array (eval cumulative product) */
			for (j = SYCOUT_GREEN_MAXDIMS-2; j >= 0; j--) {
				if (sycout_green_format[j] == SYCOUT_GREEN_NONE) continue;
				else {
                    sycout_green_factors[j] *= sycout_green_factors[j+1];
                }
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

/**
 * Compute the index of the Green's function
 * corresponding to the current integration-space
 * point.
 *
 * data:     Step data (containing information about the integrand)
 * index:    Contains the Green's function index on return
 * wavindex: Contains the Green's function index of the first wavelength on return.
 */
void sycout_green_get_index(struct sycout_data *data, size_t *index, size_t *wavindex) {
	size_t I = (size_t)(data->i*sycout_green_pixels),
		J = (size_t)(data->j*sycout_green_pixels);

	if (sycout_green_pixels > 0) {
		if (I < sycout_green_suboffseti || I >= sycout_green_subpixels+sycout_green_suboffseti)
			return;
		if (J < sycout_green_suboffsetj || J >= sycout_green_subpixels+sycout_green_suboffsetj)
			return;

		I -= sycout_green_suboffseti;
		J -= sycout_green_suboffsetj;
	}
	
	/* Compute index */
	size_t i;
	for (i = 0; i < SYCOUT_GREEN_MAXDIMS && sycout_green_format[i] != SYCOUT_GREEN_NONE; i++) {
		switch (sycout_green_format[i]) {
			case SYCOUT_GREEN_VEL1: *index += sycout_green_ivel1*sycout_green_factors[i]; break;
			case SYCOUT_GREEN_VEL2: *index += sycout_green_ivel2*sycout_green_factors[i]; break;
			case SYCOUT_GREEN_IMAGEI: *index += I*sycout_green_factors[i]; break;
			case SYCOUT_GREEN_IMAGEJ: *index += J*sycout_green_factors[i]; break;
			case SYCOUT_GREEN_RADIUS: *index += sycout_green_irad*sycout_green_factors[i]; break;
			case SYCOUT_GREEN_SPECTRUM: *wavindex = i; break;
			default: fprintf(
					stderr,
					"ERROR: (sycout green): Unrecognized GF dimension: %d. Ignoring when computing index.\n",
					sycout_green_format[i]
				);
				break;
		}
	}
}

void sycout_green_step(struct sycout_data *data) {
    size_t index = 0, wavindex = 0, i;
    sycout_green_get_index(data, &index, &wavindex);

	/* Store Stokes params? There are four
	 * parameters per parameter-space point,
	 * so multiply index by 4. */
	if (sycout_green_stokesparams)
		index *= 4;

	/* Compute differential element */
	double diffel = data->RdPhi * data->Jdtdrho * data->Jp;

	if (sycout_green_weighWdf) diffel *= data->distribution_function;
	if (!sycout_green_hasvel1 && particles_get_dvel1() != 0) diffel *= fabs(particles_get_dvel1());
	if (!sycout_green_hasvel2 && particles_get_dvel2() != 0) diffel *= fabs(particles_get_dvel2());
	if (sycout_green_hasrho && particles_get_drho() != 0) diffel /= fabs(particles_get_drho());
	/* Also, if we have both of d^2p we must multiply by the momentum space Jacobian */
	if (!sycout_green_hasvel1 && !sycout_green_hasvel2) {
		diffel *= particles_get_differential_factor_current();
	}

	/* Set value of Green's function */
	if (sycout_green_haswav) {
        double *spectrum=NULL, **polspec=NULL;
        if (sycout_green_stokesparams)
            polspec = sycamera_get_polarization_spectrum();
        else
            spectrum = sycamera_get_spectrum();

		/* If we have to weigh with the distribution function,
		 * then any or all of rho, vel1 and vel2 are NOT part
		 * of the Green's function. That implies that two threads
		 * can't work on the same element of the GF at a time,
		 * which requires 'omp critical' pragmas. */
		if (sycout_green_weighWdf) {
			#pragma omp critical
			{
                if (sycout_green_stokesparams) {
                    for (i = 0; i < sycout_green_nwav; i++) {
                        sycout_green_func[index+4*i+0] += polspec[0][i] * diffel;
                        sycout_green_func[index+4*i+1] += polspec[1][i] * diffel;
                        sycout_green_func[index+4*i+2] += polspec[2][i] * diffel;
                        sycout_green_func[index+4*i+3] += polspec[3][i] * diffel;
                    }
                } else {
                    for (i = 0; i < sycout_green_nwav; i++) {
                        sycout_green_func[index+i] += spectrum[i] * diffel;
                        index += sycout_green_factors[wavindex];
                    }
                }
			}
		} else {
            if (sycout_green_stokesparams) {
                for (i = 0; i < sycout_green_nwav; i++) {
                    sycout_green_func[index+4*i+0] += polspec[0][i] * diffel;
                    sycout_green_func[index+4*i+1] += polspec[1][i] * diffel;
                    sycout_green_func[index+4*i+2] += polspec[2][i] * diffel;
                    sycout_green_func[index+4*i+3] += polspec[3][i] * diffel;
                }
            } else {
                for (i = 0; i < sycout_green_nwav; i++) {
                    sycout_green_func[index+i] += spectrum[i] * diffel;
                    index += sycout_green_factors[wavindex];
                }
            }
		}
	} else {
        double *stokp=NULL;
        if (sycout_green_stokesparams)
            stokp = sycamera_get_polarization();

		/* If we have to weigh with the distribution function,
		 * then any or all of rho, vel1 and vel2 are NOT part
		 * of the Green's function. That implies that two threads
		 * can't work on the same element of the GF at a time,
		 * which requires 'omp critical' pragmas. */
		if (sycout_green_weighWdf) {
			#pragma omp critical
			{
                if (sycout_green_stokesparams) {
                    sycout_green_func[index+0] += stokp[0];
                    sycout_green_func[index+1] += stokp[1];
                    sycout_green_func[index+2] += stokp[2];
                    sycout_green_func[index+3] += stokp[3];
                } else
                    sycout_green_func[index] += data->brightness * diffel;
			}
		} else {
            if (sycout_green_stokesparams) {
                sycout_green_func[index+0] += stokp[0];
                sycout_green_func[index+1] += stokp[1];
                sycout_green_func[index+2] += stokp[2];
                sycout_green_func[index+3] += stokp[3];
            } else {
                sycout_green_func[index] += data->brightness * diffel;
            }
		}
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

		/* Format of function */
		char formats[SYCOUT_GREEN_MAXDIMS+1];
		int haspixels = 0;
		for (i = 0; i < SYCOUT_GREEN_MAXDIMS && sycout_green_format[i] != SYCOUT_GREEN_NONE; i++) {
			switch (sycout_green_format[i]) {
				case SYCOUT_GREEN_VEL1:
					formats[i] = '1';
					break;
				case SYCOUT_GREEN_VEL2:
					formats[i] = '2';
					break;
				case SYCOUT_GREEN_RADIUS:
					formats[i] = 'r';
					break;
				case SYCOUT_GREEN_IMAGEI:
					haspixels = 1;
					formats[i] = 'i';
					break;
				case SYCOUT_GREEN_IMAGEJ:
					haspixels = 1;
					formats[i] = 'j';
					break;
				case SYCOUT_GREEN_SPECTRUM:
					formats[i] = 'w';
					break;
				default:
					fprintf(
						stderr,
						"WARNING: (sycout green): Unrecognized GF dimension: %d. Ignoring during output.\n",
						sycout_green_format[i]
					);
					break;
			}
		}

		formats[i] = 0;
		sf->write_string(sf, "format", formats, i);

		/* Number of pixels */
		if (!haspixels) {
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
		if (sycout_green_hasrho) {
			v = malloc(sizeof(double)*b[2]);

			for (i = 0; i < b[2]; i++)
				v[i] = b[0] + ((double)i)/((double)b[2]) * (b[1]-b[0]);

			sf->write_list(sf, "r", v, b[2]);
			free(v);
		} else {
			double d = 0;
			sf->write_list(sf, "r", &d, 1);
		}

		/* param1 vector */
		if (sycout_green_hasvel1) {
			v = malloc(sizeof(double)*b[5]);

			for (i = 0; i < b[5]; i++)
				v[i] = b[3] + ((double)i)/((double)b[5]) * (b[4]-b[3]);

			sf->write_list(sf, "param1", v, b[5]);
			free(v);
		} else {
			double d = 0;
			sf->write_list(sf, "param1", &d, 1);
		}

		/* param2 vector */
		if (sycout_green_hasvel2) {
			v = malloc(sizeof(double)*b[8]);

			for (i = 0; i < b[8]; i++)
				v[i] = b[6] + ((double)i)/((double)b[8]) * (b[7]-b[6]);

			sf->write_list(sf, "param2", v, b[8]);
			free(v);
		} else {
			double d = 0;
			sf->write_list(sf, "param2", &d, 1);
		}

		/* Wavelength vector */
		if (!sycout_green_haswav) {
			double wavelengths = 0;
			sf->write_list(sf, "wavelengths", &wavelengths, 1);
		} else {
			v = sycamera_get_wavelengths();
			sf->write_list(sf, "wavelengths", v, sycout_green_nwav);
		}

		/* Actual Green's function */
		sf->write_list(sf, "func", sycout_green_func, sycout_green_func_sz);

        /* Write Stokes parameter flag */
        sf->write_attribute_scalar(sf, "func", "stokesparameters", sycout_green_stokesparams);

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

