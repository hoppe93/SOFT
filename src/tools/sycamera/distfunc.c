/* Distribution function management */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "distfunc.h"
#include "sfile.h"

distfunc *distfunc_function;

void distfunc_init_run(void) {
    df_interp_init_run();
}
double distfunc_eval(double r, double xi, double p) {
    return df_interp_eval(r,xi,p);
}
void distfunc_load(const char *filename) {
	sFILE *sf;
	enum sfile_type st;
	sfilesize_t dims[2];
	double **temp;
	char *punits;
	size_t i;
	
	st = sfile_get_filetype(filename);
	if (st == FILETYPE_UNKNOWN) {
		fprintf(stderr, "ERROR: (distfunc): Unrecognized type of distribution function file: %s.\n", filename);
		exit(EXIT_FAILURE);
	}

	sf = sfile_init(st);
	if (!sf->open(sf, filename, SFILE_MODE_READ)) {
		fprintf(stderr, "ERROR: (distfunc): Unable to open distribution function: %s.\n", filename);
		exit(EXIT_FAILURE);
	}

    distfunc_function = malloc(sizeof(distfunc));

	/* Load name, description and units of p */
	distfunc_function->name = sf->get_string(sf, "name");
	distfunc_function->desc = sf->get_string(sf, "description");
	punits = sf->get_string(sf, "punits");
	if (distfunc_function->name == NULL) { fprintf(stderr, "ERROR: (distfunc): The distribution function does not have any name.\n"); exit(EXIT_FAILURE); }
	if (distfunc_function->desc == NULL) { fprintf(stderr, "ERROR: (distfunc): The distribution function does not contain any description.\n"); exit(EXIT_FAILURE); }
	if (punits == NULL) { fprintf(stderr, "ERROR: (distfunc): No unit has been specified for the momentum.\n"); exit(EXIT_FAILURE); }

	/* Load coordinates */
	temp = sf->get_doubles(sf, "r", dims);
	if (temp == NULL) { fprintf(stderr, "ERROR: (distfunc): The distribution function does not contain an 'r' vector.\n"); exit(EXIT_FAILURE); }
	distfunc_function->nr   = dims[1];
	distfunc_function->r = *temp; free(temp);

	temp = sf->get_doubles(sf, "xi", dims);
	if (temp == NULL) { fprintf(stderr, "ERROR: (distfunc): The distribution function does not contain a 'xi' vector.\n"); exit(EXIT_FAILURE); }
	distfunc_function->nxi = dims[1];
	distfunc_function->xi = *temp; free(temp);

	temp = sf->get_doubles(sf, "p", dims);
	if (temp == NULL) { fprintf(stderr, "ERROR: (distfunc): The distribution function does not contain a 'p' vector.\n"); exit(EXIT_FAILURE); }
	distfunc_function->np   = dims[1];
	distfunc_function->p = *temp; free(temp);

	distfunc_function->rmin 	= distfunc_function->r[0];
	distfunc_function->rmax 	= distfunc_function->r[distfunc_function->nr-1];
	distfunc_function->ximin 	= distfunc_function->xi[0];
	distfunc_function->ximax 	= distfunc_function->xi[distfunc_function->nxi-1];
	distfunc_function->pmin 	= distfunc_function->p[0];
	distfunc_function->pmax 	= distfunc_function->p[distfunc_function->np-1];

	if (!strcmp(punits, "ev")) {
		distfunc_function->pmin *= MOMENTUM;
		distfunc_function->pmax *= MOMENTUM;
		for (i = 0; i < distfunc_function->np; i++) {
			distfunc_function->p[i] *= MOMENTUM;
		}
	} else if (!strcmp(punits, "normalized")) {
		distfunc_function->pmin *= ELECTRONMASS*LIGHTSPEED;
		distfunc_function->pmax *= ELECTRONMASS*LIGHTSPEED;
		for (i = 0; i < distfunc_function->np; i++) {
			distfunc_function->p[i] *= ELECTRONMASS*LIGHTSPEED;
		}
	} else if (!strcmp(punits, "si")) {
		/* This is what SOFT uses internally,
		 * so no need to do anything... */
	} else {
		fprintf(stderr, "ERROR: (distfunc): Unrecognized units of momentum: '%s'.\n", punits);
		exit(EXIT_FAILURE);
	}

	/* Initialize first parts of distribution function */
	distfunc_function->value = sf->get_doubles(sf, "f", dims);
	if (distfunc_function->value == NULL) { fprintf(stderr, "ERROR: (distfunc): The distribution function does not contain an 'f' matrix.\n"); exit(EXIT_FAILURE); }
	if (dims[0] != distfunc_function->nr || dims[1] != distfunc_function->nxi*distfunc_function->np) {
		fprintf(
			stderr,
			"ERROR: (distfunc): Dimensions of grid (r x xi x p = %zu x %zu x %zu) and distribution function (f) do not match.\n",
			distfunc_function->nr, distfunc_function->nxi, distfunc_function->np
		);
		exit(EXIT_FAILURE);
	}

	sf->close(sf);

    df_interp_init(distfunc_function);
    printf("Loaded distribution function from '%s'.\n", filename);
    df_readfile_unload();
}

void distfunc_test(void) {
    distfunc_load("../../softdist/analytical/slim.func");

    /* Check r */
    size_t i;
    for (i = 0; i < distfunc_function->nr-1; i++) {
        if (distfunc_function->r[i] >= distfunc_function->r[i+1]) {
            printf("Error: Radius must be strictly increasing! Error found in index %zu of radius vector (%e >= %e).\n", i,
                distfunc_function->r[i], distfunc_function->r[i+1]);
            exit(-1);
        }
    }

    for (i = 0; i < distfunc_function->nxi-1; i++) {
        if (distfunc_function->xi[i] >= distfunc_function->xi[i+1]) {
            printf("Error: cos(theta) must be strictly increasing! Error found in index %zu of cos(theta) vector (%e >= %e).\n", i,
                distfunc_function->xi[i], distfunc_function->xi[i+1]);
            exit(-1);
        }
    }
    for (i = 0; i < distfunc_function->np-1; i++) {
        if (distfunc_function->p[i] >= distfunc_function->p[i+1]) {
            printf("Error: Momentum must be strictly increasing! Error found in index %zu of momentum vector (%e >= %e).\n", i,
                distfunc_function->p[i], distfunc_function->p[i+1]);
            exit(-1);
        }
    }

    /* Check values */
	/*
    int j,k;
    for (i = 0; i < distfunc_function->nr; i++) {
        for (j = 0; j < distfunc_function->nxi*distfunc_function->np; j++) {
			for (k = 0; k < distfunc_function->np; k++) {
				if (distfunc_function->value[i][j*distfunc_function->nxi+k] != 1.) {
					printf("Warning: %d:%d:%d: Value does not equal 1\n", i, j, k);
				}
			}
        }
    }
	*/

	i = distfunc_function->nr-1;
	size_t j=distfunc_function->nxi-5;
	size_t k=distfunc_function->np-3;
	printf("f at index %zu:%zu:%zu = %e\n", i, j, k, distfunc_function->value[i][j*distfunc_function->np+k]);

	double r = distfunc_function->r[i];
	double c = distfunc_function->xi[j];
	double p = distfunc_function->p[k];
	double dc = distfunc_function->xi[j]-distfunc_function->xi[j-1];
	double nc = c-dc;
	printf("%zu:%zu:%zu <=> r = %e, cos(theta) = %e, p = %e\n", i, j, k, r, c, p);

    df_interp_init_run();

	printf("f(r,c,p) = %e\n", df_interp_eval(r, c, p));
	printf("f(r,nc,p) = %e\n", df_interp_eval(r, nc, p));
	printf("new c = %e\n", nc);

    //printf("Distribution function passed all tests.\n");

/*
    double pmin=1., pmax=50., dp=2., p,
        r = .72, xi = 0.99;
    for (p = pmin; p <= pmax; p += dp) {
        printf("(r = %e, p = %e, cos(theta) = %e) = %e\n", r, p, xi, distfunc_eval(r, xi, p));
    }
*/
}
