/* Calculation of the synchrotron power emitted in a given
   wavelength range. The function Pcyl was originally given
   by [Bekefi, 1966] and implemented in MATLAB by A. Stahl.
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <omp.h>
#include "sycamera.h"

double *sycamera_pcyl_lambdas, sycamera_pcyl_dl, *sycamera_pcyl_comp;
int sycamera_pcyl_lambda_resolution, sycamera_pcyl_haswarned;
gsl_interp_accel *sycamera_pcyl_acc;
gsl_spline *sycamera_pcyl_spline;

#pragma omp threadprivate(sycamera_pcyl_haswarned,sycamera_pcyl_acc,sycamera_pcyl_spline,sycamera_pcyl_comp)

void sycamera_pcyl_init(double lambda0, double lambda1, int res) {
    sycamera_pcyl_haswarned = 0;

    sycamera_pcyl_lambdas = malloc(sizeof(double)*res);
    sycamera_pcyl_lambda_resolution = res;
    int i;
    for (i = 0; i < res; i++) {
        /* We reverse this spectrum for simpler code below */
        sycamera_pcyl_lambdas[res-i-1] = lambda0 + i * (lambda1-lambda0) / (res-1);
    }

    /* Compute step length */
    sycamera_pcyl_dl = fabs(sycamera_pcyl_lambdas[1]-sycamera_pcyl_lambdas[0]);
}
void sycamera_pcyl_init_run(void) {
	int i;
	/* Initialize spectrum storage */
	sycamera_pcyl_comp = malloc(sizeof(double)*sycamera_pcyl_lambda_resolution);	/* Storage of spectrum components */
	for (i = 0; i < sycamera_pcyl_lambda_resolution; i++) {
		sycamera_pcyl_comp[i] = 0;
	}

	/* Initialize interpolation */
    /* Initialize power spline */
    sycamera_pcyl_acc = gsl_interp_accel_alloc();

    /* Linear interpolation (somewhat faster, less accurate) */
    //sycamera_pcyl_spline = gsl_spline_alloc(gsl_interp_linear, sycamera_pcyl_lookup_count);
    /* Cubic spline interpolation */
    sycamera_pcyl_spline = gsl_spline_alloc(gsl_interp_cspline, sycamera_pcyl_lookup_count);

    gsl_spline_init(sycamera_pcyl_spline, sycamera_pcyl_lookup_lambda, sycamera_pcyl_lookup_int, sycamera_pcyl_lookup_count);
}

/**
 * Computes the integral of P_cyl from lambda0 to lambda1
 *
 * lambda0: The lower wavelength to evaluate the integral at
 * lambda1: The upper wavelength to evaluate the integral at
 * p: Particle momentum
 * xi: Cosine of the pitch angle
 * Bmag: Magnetic field strength in the point of evaluation
 */
double sycamera_pcyl_int(double ppar2, double pperp2, double Bmag, double mass, double fraction) {
    double ic = 1/((LIGHTSPEED*LIGHTSPEED)*mass*mass);
    double ppar2l = ppar2*ic;
    double pperp2l = pperp2*ic;
    double gamma2 = 1+ppar2l+pperp2l;
	double gammaParallel2 = gamma2 / (1+pperp2l);
    double gammaParallel = sqrt(gamma2 / (1+pperp2l));
	double betaparpar = ppar2l / sqrt(gamma2*gamma2 - gamma2);	// = beta*cospitch*cospitch
    double lowerBound = SYCAMERA_PCYL_LOWERBOUND * gammaParallel * mass / (Bmag*gamma2);
    double prefactor = SYCAMERA_PCYL_CONSTANT * gammaParallel2  * (1 - betaparpar) / gamma2;

    int i;
    double sum=0., ilambda, lb, pf;
    for (i = 0; i < sycamera_pcyl_lambda_resolution; i++) {
        ilambda = 1 / sycamera_pcyl_lambdas[i];
        lb = lowerBound * ilambda;
        pf = prefactor*ilambda*ilambda*ilambda;

        /* Lambdas outside the bounds are assumed negligible */
        if (lb >= sycamera_pcyl_lookup_lambda[sycamera_pcyl_lookup_count-1]) {
			sycamera_pcyl_comp[i] = 0;
            continue;
        } else if (lb <= sycamera_pcyl_lookup_lambda[0]) {
            if (!sycamera_pcyl_haswarned) {
                fprintf(stderr, "WARNING: The detector spectrum used is outside the supported bounds. Please adjust the lower bound to > %e.\n", lowerBound/(sycamera_pcyl_lookup_lambda[0]));
                sycamera_pcyl_haswarned = 1;
            }
			sycamera_pcyl_comp[i] = 0;
            continue;
        }

		double e = gsl_spline_eval(sycamera_pcyl_spline, lb, sycamera_pcyl_acc);
		double c = pf * e * fraction;
		sycamera_pcyl_comp[i] = c;
        sum += c;
    }

    return sum * sycamera_pcyl_dl;
}

double *sycamera_pcyl_get_wavelengths(void) { return sycamera_pcyl_lambdas; }
double *sycamera_pcyl_get_spectrum(void) { return sycamera_pcyl_comp; }
int sycamera_pcyl_get_spectrum_length(void) { return sycamera_pcyl_lambda_resolution; }

void sycamera_pcyl_test(void) {
/*
    sycamera_pcyl_init(5e-7,1e-6,50);
    //double p2[8] = {1,25,100,400,625,900,1600,2500};
    double sin2theta = 0.009900990099010;
    double Bmag = 4.923295e+00;
    int branch = 0;

    if (branch == 0) {
        double p2=1.603286e-20*1.603286e-20;
        double sum = sycamera_pcyl_int(p2*(1-sin2theta), p2*sin2theta, Bmag, 9.10938e-31);
        printf("p2 = %e   ==>   %e\n", p2, sum);
    } else {
        double p2=1;
        unsigned int i;
        for (i = 0; i < 12345678; i++) {
            p2 += 0.0001;
            sycamera_pcyl_int(p2*(LIGHTSPEED*LIGHTSPEED), sin2theta, Bmag, 9.10938e-31);
            //printf("p2 = %e   ==>   %e\n", p2, sum);
        }
    }
}
void sycamera_pdist_test(void) {
*/
	double B = 5.0;
	double lambda0 = 1e-6,
		   lambda1 = 4e-7;
	double ppar = 4.95e7 * 5.36e-28;
	double pperp = 6.48e6 * 5.36e-28;
	double mass = 9.10938356e-31;

	double ppar2 = ppar*ppar, pperp2 = pperp*pperp;

	FILE *f;

	sycamera_pcyl_init_run();

	f = fopen("pcyl.out", "w");
	if (!f) {
		perror("ERROR");
		fprintf(stderr, "Unable to generate output file. Exiting...\n");
		exit(-1);
	}

	int i, N = 1;
	double mu = 0.0, dmu = 0.3 / (N-1),/* s2mu = 0.0,*/
		intensity = 0.0;
	for (i = 0; i < N; i++, mu += dmu) {
		sycamera_pcyl_init(lambda1,lambda0,50);

		//intensity = sycamera_pdist_int(1 / gamma2, sin2p, beta2*cos2p, beta2, B, s2mu);
        intensity = sycamera_pcyl_int(ppar2, pperp2, B, mass, 1.0);

		fprintf(f, "%.12e \t %.12e\n", mu, intensity);
	}

	fclose(f);
}

