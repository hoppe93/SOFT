/* Distribution function management */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "distfunc.h"

distfunc *distfunc_function;

void distfunc_init_run(void) {
    df_interp_init_run();
}
double distfunc_eval(double r, double costheta, double p) {
    return df_interp_eval(r,costheta,p);
}
void distfunc_load(const char *filename) {
    distfunc_function = malloc(sizeof(distfunc));

    df_readfile_load(filename);
    df_readfile_reset();

	/* First, we read of the array length specifiers */
	df_readfile_loaddbl(&(distfunc_function->rmin)); /* rmin */
	df_readfile_loaddbl(&(distfunc_function->rmax)); /* rmax */
	df_readfile_loadint(&(distfunc_function->nr)); /* nr */
	if (!df_readfile_stopped_at_newline()) { fprintf(stderr, "ERROR: Expected newline after nr!\n"); exit(EXIT_FAILURE); }

	df_readfile_loaddbl(&(distfunc_function->costmin)); /* costmin */
	df_readfile_loaddbl(&(distfunc_function->costmax)); /* costmax */
	df_readfile_loadint(&(distfunc_function->ncostheta)); /* ncostheta */
	if (!df_readfile_stopped_at_newline()) { fprintf(stderr, "ERROR: Expected newline after ntheta!\n"); exit(EXIT_FAILURE); }

	df_readfile_loaddbl(&(distfunc_function->pmin)); /* pmin */
	df_readfile_loaddbl(&(distfunc_function->pmax)); /* pmax */
	df_readfile_loadint(&(distfunc_function->np)); /* np */
	if (!df_readfile_stopped_at_newline()) { fprintf(stderr, "ERROR: Expected newline after np!\n"); exit(EXIT_FAILURE); }

	/* Allocate memory for coordinate holders */
	distfunc_function->r 	= malloc(sizeof(double)*distfunc_function->nr);
	distfunc_function->costheta = malloc(sizeof(double)*distfunc_function->ncostheta);
	distfunc_function->p 	= malloc(sizeof(double)*distfunc_function->np);
	//distfunc_function->value = malloc(sizeof(double*)*distfunc_function->nr);

    /* Load r values */
	int i=0;
    do {
        df_readfile_loaddbl(distfunc_function->r+i); i++;
    } while (!df_readfile_eof() && !df_readfile_stopped_at_newline());

    /* Load cos(theta) values */
    i=0;
    do {
        df_readfile_loaddbl(distfunc_function->costheta+i); i++;
    } while (!df_readfile_eof() && !df_readfile_stopped_at_newline());

	/* Load p values */
    i=0;
    do {
		df_readfile_loaddbl(distfunc_function->p+i); i++;
	} while (!df_readfile_eof() && !df_readfile_stopped_at_newline());

    /*
	double dr=0., dp=0., dcost=0.;

	if (distfunc_function->nr > 1)
		dr = (distfunc_function->rmax - distfunc_function->rmin)/(distfunc_function->nr-1);
	if (distfunc_function->ncostheta > 1)
		dcost = (distfunc_function->costmax - distfunc_function->costmin)/(distfunc_function->ncostheta-1);
	if (distfunc_function->np > 1)
		dp = (distfunc_function->pmax - distfunc_function->pmin)/(distfunc_function->np-1);
    */

	/* Load r, theta and f values */
	int r=0, t=0, p=0;

	/* Initialize first parts of distribution function */
	distfunc_function->value = malloc(sizeof(double*)*distfunc_function->nr);
	distfunc_function->value[0] = malloc(sizeof(double)*distfunc_function->ncostheta*distfunc_function->np);

	while (!df_readfile_eof()) {
		/* First, load r value (if the value is not repeated) */
		//distfunc_function->costheta[t] = distfunc_function->costmin + dcost*t;

		/* Load function values until end-of-line */
		p = 0;
		do {
			//distfunc_function->p[p] = distfunc_function->pmin + dp*p;
			df_readfile_loaddbl(distfunc_function->value[r]+t*distfunc_function->np+p);
			p++;
		} while (!df_readfile_stopped_at_newline());

		i++, t++;
		if (t >= distfunc_function->ncostheta) {
			t = 0;
			//distfunc_function->r[r] = distfunc_function->rmin + dr*r;
			r++;
			if (r == distfunc_function->nr) break;

			distfunc_function->value[r] = malloc(sizeof(double)*distfunc_function->ncostheta*distfunc_function->np);
		}
	}

    df_interp_init(distfunc_function);
    printf("Loaded distribution function from '%s'.\n", filename);
    df_readfile_unload();
}

void distfunc_test(void) {
    distfunc_load("../../softdist/analytical/slim.func");

    /* Check r */
    int i;
    for (i = 0; i < distfunc_function->nr-1; i++) {
        if (distfunc_function->r[i] >= distfunc_function->r[i+1]) {
            printf("Error: Radius must be strictly increasing! Error found in index %d of radius vector (%e >= %e).\n", i,
                distfunc_function->r[i], distfunc_function->r[i+1]);
            exit(-1);
        }
    }

    for (i = 0; i < distfunc_function->ncostheta-1; i++) {
        if (distfunc_function->costheta[i] >= distfunc_function->costheta[i+1]) {
            printf("Error: cos(theta) must be strictly increasing! Error found in index %d of cos(theta) vector (%e >= %e).\n", i,
                distfunc_function->costheta[i], distfunc_function->costheta[i+1]);
            exit(-1);
        }
    }
    for (i = 0; i < distfunc_function->np-1; i++) {
        if (distfunc_function->p[i] >= distfunc_function->p[i+1]) {
            printf("Error: Momentum must be strictly increasing! Error found in index %d of momentum vector (%e >= %e).\n", i,
                distfunc_function->p[i], distfunc_function->p[i+1]);
            exit(-1);
        }
    }

    /* Check values */
	/*
    int j,k;
    for (i = 0; i < distfunc_function->nr; i++) {
        for (j = 0; j < distfunc_function->ncostheta*distfunc_function->np; j++) {
			for (k = 0; k < distfunc_function->np; k++) {
				if (distfunc_function->value[i][j*distfunc_function->ncostheta+k] != 1.) {
					printf("Warning: %d:%d:%d: Value does not equal 1\n", i, j, k);
				}
			}
        }
    }
	*/

	i = distfunc_function->nr-1;
	int j=distfunc_function->ncostheta-5;
	int k=distfunc_function->np-3;
	printf("f at index %d:%d:%d = %e\n", i, j, k, distfunc_function->value[i][j*distfunc_function->np+k]);

	double r = distfunc_function->r[i];
	double c = distfunc_function->costheta[j];
	double p = distfunc_function->p[k];
	double dc = distfunc_function->costheta[j]-distfunc_function->costheta[j-1];
	double nc = c-dc;
	printf("%d:%d:%d <=> r = %e, cos(theta) = %e, p = %e\n", i, j, k, r, c, p);

    df_interp_init_run();

	printf("f(r,c,p) = %e\n", df_interp_eval(r, c, p));
	printf("f(r,nc,p) = %e\n", df_interp_eval(r, nc, p));
	printf("new c = %e\n", nc);

    //printf("Distribution function passed all tests.\n");

/*
    double pmin=1., pmax=50., dp=2., p,
        r = .72, costheta = 0.99;
    for (p = pmin; p <= pmax; p += dp) {
        printf("(r = %e, p = %e, cos(theta) = %e) = %e\n", r, p, costheta, distfunc_eval(r, costheta, p));
    }
*/
}
