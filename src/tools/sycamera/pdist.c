/**
 * Calculation of the radiation emitted in a given
 * frequency interval at a given viewing angle.
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "sycamera.h"

double *sycamera_pdist_omegas, sycamera_pdist_dl,
	sycamera_pdist_prefactor, sycamera_pdist_omega_B_factor,
	sycamera_pdist_omega_low, sycamera_pdist_omega_up,
	sycamera_pdist_ximax, sycamera_pdist_ximin;
int sycamera_pdist_resolution, sycamera_pdist_haswarned;
gsl_interp_accel *sycamera_pdist_acc1, *sycamera_pdist_acc2;
gsl_spline *sycamera_pdist_spline1, *sycamera_pdist_spline2;
enum sycamera_polarization_type sycamera_pdist_polt;

#pragma omp threadprivate(sycamera_pdist_haswarned,sycamera_pdist_acc1,sycamera_pdist_acc2,sycamera_pdist_spline1,sycamera_pdist_spline2)

void sycamera_pdist_init(double omega0, double omega1, enum sycamera_polarization_type polt) {
	sycamera_pdist_omega_low = omega0;
	sycamera_pdist_omega_up = omega1;
	sycamera_pdist_polt = polt;
}
void sycamera_pdist_init_run(void) {
	sycamera_pdist_acc1 = gsl_interp_accel_alloc();
	sycamera_pdist_acc2 = gsl_interp_accel_alloc();
	sycamera_pdist_spline1 = gsl_spline_alloc(gsl_interp_cspline, sycamera_pdist_lookup_count);
	sycamera_pdist_spline2 = gsl_spline_alloc(gsl_interp_cspline, sycamera_pdist_lookup_count);

	gsl_spline_init(sycamera_pdist_spline1, sycamera_pdist_lookup_omega, sycamera_pdist_lookup_int1, sycamera_pdist_lookup_count);
	gsl_spline_init(sycamera_pdist_spline2, sycamera_pdist_lookup_omega, sycamera_pdist_lookup_int2, sycamera_pdist_lookup_count);

	sycamera_pdist_ximax = sycamera_pdist_lookup_omega[sycamera_pdist_lookup_count-1];
	sycamera_pdist_ximin = sycamera_pdist_lookup_omega[0];
}
void sycamera_pdist_init_particle(double mass) {
	sycamera_pdist_prefactor = 9.0 * CHARGE*CHARGE*CHARGE*CHARGE / (256.0 * PI*PI*PI * EPS0 * LIGHTSPEED * mass*mass);
	sycamera_pdist_omega_B_factor = CHARGE/mass;
	//sycamera_pdist_omega_c_factor = 3.0 * CHARGE / (2.0*mass);
}

double sycamera_pdist_int(
	double gammai2, double gamma3, double gammapar2, double beta,
	double beta2, double betapar2, double Bmag, double sinmu,
	double cosmu, double sinp, double cosp
) {
	double
		omegaB = sycamera_pdist_omega_B_factor*Bmag*sqrt(gammai2),
		cospsi = cosmu*cosp + sinmu*sinp,
		bcospsi = beta*cospsi,
		bcospsi2 = bcospsi/2.0,
		mcospsi = 1-bcospsi,
		mcospsi2 = mcospsi*mcospsi,
		gpar2mcospsi = gammapar2*mcospsi;
	
	double pf = sycamera_pdist_prefactor * Bmag*Bmag * beta2 *
		(1 - beta*cosp*cosmu) * gammai2 / (sqrt(gpar2mcospsi*bcospsi2)*mcospsi2);

	double cf = 2.0 / (3.0*omegaB) * sqrt(gpar2mcospsi*mcospsi2/bcospsi2);
	double lower = sycamera_pdist_omega_low * cf;
	double upper = sycamera_pdist_omega_up * cf;

	if (lower < sycamera_pdist_ximin)
		lower = sycamera_pdist_ximin;
	else if (lower > sycamera_pdist_ximax)
		lower = sycamera_pdist_ximax;

	if (upper < sycamera_pdist_ximin)
		upper = sycamera_pdist_ximin;
	else if (upper > sycamera_pdist_ximax)
		upper = sycamera_pdist_ximax;
	
	/* Compute integrals */
	double I13l, I13u, I23l, I23u, I13, I23;

	switch (sycamera_pdist_polt) {
		case SYCAMERA_POLARIZATION_PARALLEL:
			I23l = gsl_spline_eval(sycamera_pdist_spline2, lower, sycamera_pdist_acc2);
			I23u = gsl_spline_eval(sycamera_pdist_spline2, upper, sycamera_pdist_acc2);
			I23 = I23l - I23u;

			return pf * I23;
		case SYCAMERA_POLARIZATION_PERPENDICULAR:
			I13l = gsl_spline_eval(sycamera_pdist_spline1, lower, sycamera_pdist_acc1);
			I13u = gsl_spline_eval(sycamera_pdist_spline1, upper, sycamera_pdist_acc1);
			I13 = I13l - I13u;

			return pf * bcospsi2/mcospsi * (1-cospsi*cospsi) * I13;
		case SYCAMERA_POLARIZATION_BOTH:
		default:
			I13l = gsl_spline_eval(sycamera_pdist_spline1, lower, sycamera_pdist_acc1);
			I13u = gsl_spline_eval(sycamera_pdist_spline1, upper, sycamera_pdist_acc1);
			I23l = gsl_spline_eval(sycamera_pdist_spline2, lower, sycamera_pdist_acc2);
			I23u = gsl_spline_eval(sycamera_pdist_spline2, upper, sycamera_pdist_acc2);
		
			I13 = I13l - I13u;
			I23 = I23l - I23u;
		
			return pf * (I23 + bcospsi2/mcospsi * (1-cospsi*cospsi) * I13);
	}
}

void sycamera_pdist_test(void) {
	double B = 5.0;
	double omega0 = 2.0*PI*LIGHTSPEED / 1e-6,
		   omega1 = 2.0*PI*LIGHTSPEED / 4e-7;
	double ppar = 1e8 * 5.36e-28;
	double pperp = 1e6 * 5.36e-28;
	double mass = 9.10938356e-31;

	double ppar2 = ppar*ppar, pperp2 = pperp*pperp;
	double p2 = ppar2 + pperp2;
	double gamma2 = 1 + p2 / (mass*mass*LIGHTSPEED*LIGHTSPEED);
	double beta2 = 1 - 1 / gamma2;
	double sin2p = pperp2 / p2, cos2p = ppar2 / p2;

	FILE *f;

	sycamera_pdist_init(omega0, omega1, SYCAMERA_POLARIZATION_BOTH);
	sycamera_pdist_init_run();
	sycamera_pdist_init_particle(mass);

	f = fopen("pdist.out", "w");
	if (!f) {
		perror("ERROR");
		fprintf(stderr, "Unable to generate output file. Exiting...\n");
		exit(-1);
	}

	int i, N = 1000;
	double mu = 0.0, dmu = 0.1 / (N-1), smu=0.0, cmu=0.0,
		intensity = 0.0;
	for (i = 0; i < N; i++, mu += dmu) {
		smu = sin(mu);
		cmu = cos(mu);

		intensity = sycamera_pdist_int(
			1/gamma2, gamma2*sqrt(gamma2), 1/(1-beta2*cos2p), sqrt(beta2), beta2, beta2*cos2p,
			B, smu, cmu, sqrt(sin2p), sqrt(cos2p)
		);

		fprintf(f, "%.12e \t %.12e\n", mu, intensity);
	}

	fclose(f);
}

