/* Distribution function management */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "distfunc.h"
#include "particles.h"
#include "sfile.h"

distfunc *distfunc_function;
int distfunc_include_drifts=0;

void distfunc_init_run(int include_drifts) {
    distfunc_include_drifts = include_drifts;
    df_interp_init_run();
}
double distfunc_eval(double r, double xi, double p) {
    if (distfunc_include_drifts) {
        double dr = particles_get_orbit_drift_shift();
        return df_interp_eval(r-dr, xi, p);
    } else
        return df_interp_eval(r, xi, p);
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
    distfunc_function->logarithmic = 0;

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

    df_interp_init(distfunc_function, 1);
    printf("Loaded distribution function from '%s'.\n", filename);
    df_readfile_unload();
}

void distfunc_test(void) {
	size_t i, j, pi, xii;
    distfunc_load("/mnt/HDD/runaway/mathias/softruns/softarticle/alexdist/alexdist-rad.mat");
    //distfunc_load("/home/hoppe/Skrivbord/alexdist-rad.mat");

	pi = 17;
	xii = 7;

	/* Print p */
	printf("p = [");
	for (i = 0; i < distfunc_function->np && i < 5; i++) {
		printf("%.4e  ", distfunc_function->p[i]);
	}
	printf("...  %.4e]\n", distfunc_function->pmax);

	/* Print radial distribution */
	printf("f(r) = ");
	for (i = 0; i < distfunc_function->nr && i < 5; i++) {
		printf("%e  ", distfunc_function->value[i][0]);
	}
	printf("\n");

	/* Print first momentum distribution */
	printf("f(xi) = ");
	for (i = xii; i < distfunc_function->nxi && i < xii+5; i++) {
		printf("%e  ", distfunc_function->value[0][pi+i*distfunc_function->np]);
	}
	printf("\n");

	/* Print first momentum distribution */
	printf("f(p) = ");
	for (i = pi; i < distfunc_function->np && i < pi+5; i++) {
		printf("%e  ", distfunc_function->value[0][i+xii*distfunc_function->np]);
	}
	printf("\n");

	df_interp_init_run();

	double r = distfunc_function->r[0],
		   xi = distfunc_function->xi[xii],
		   p = distfunc_function->p[pi],
		   v = df_interp_eval(r, xi, p);
	printf("f(r = %.3f, p = %.3e, xi = %.3e) = %e\n", r, p, xi, v);

	/* Create part of matrix using interpolation */
	for (j = xii; j < xii+4; j++) {
		for (i = pi; i < pi+4; i++) {
			xi = distfunc_function->xi[j];
			p = distfunc_function->p[i];
			printf("%.3e  ", df_interp_eval(r, xi, p));
		}
		printf("\n");
	}
}
