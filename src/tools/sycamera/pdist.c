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
	sycamera_pdist_ximax, sycamera_pdist_ximin, *sycamera_pdist_polarization,
	*sycamera_pdist_spectrum, *sycamera_pdist_wavelengths,
	sycamera_pdist_spec_prefactor, **sycamera_pdist_polarization_spectrum,
	sycamera_pdist_dlambda;
int sycamera_pdist_spectrum_resolution, sycamera_pdist_haswarned;
gsl_interp_accel *sycamera_pdist_acc1, *sycamera_pdist_acc2,
				 *sycamera_pdist_spec_acc1, *sycamera_pdist_spec_acc2;
gsl_spline *sycamera_pdist_spline1, *sycamera_pdist_spline2,
		   * sycamera_pdist_spec_spline1, *sycamera_pdist_spec_spline2;

#pragma omp threadprivate(sycamera_pdist_haswarned,sycamera_pdist_acc1,sycamera_pdist_acc2, \
						  sycamera_pdist_spec_acc1,sycamera_pdist_spec_acc2,sycamera_pdist_spline1, \
						  sycamera_pdist_spline2,sycamera_pdist_spec_spline1,sycamera_pdist_spec_spline2, \
						  sycamera_pdist_polarization,sycamera_pdist_spectrum,sycamera_pdist_polarization_spectrum)

void sycamera_pdist_init(double omega0, double omega1, int spectrum_resolution) {
	sycamera_pdist_omega_low = omega0;
	sycamera_pdist_omega_up = omega1;
	sycamera_pdist_spectrum_resolution = spectrum_resolution;

	/* Generate wavelengths */
	double dlambda = 2*PI*LIGHTSPEED * (1.0/omega0 - 1.0/omega1) / ((double)spectrum_resolution);
	sycamera_pdist_wavelengths = malloc(sizeof(double)*sycamera_pdist_spectrum_resolution);
	sycamera_pdist_wavelengths[0] = 2*PI*LIGHTSPEED / omega1;

	int i;
	for (i = 1; i < sycamera_pdist_spectrum_resolution; i++) {
		sycamera_pdist_wavelengths[i] = sycamera_pdist_wavelengths[i-1] + dlambda;
	}

	sycamera_pdist_dlambda = dlambda;
}
void sycamera_pdist_init_run(void) {
	sycamera_pdist_acc1 = gsl_interp_accel_alloc();
	sycamera_pdist_acc2 = gsl_interp_accel_alloc();
	sycamera_pdist_spec_acc1 = gsl_interp_accel_alloc();
	sycamera_pdist_spec_acc2 = gsl_interp_accel_alloc();

	sycamera_pdist_spline1 = gsl_spline_alloc(gsl_interp_cspline, sycamera_pdist_lookup_count);
	sycamera_pdist_spline2 = gsl_spline_alloc(gsl_interp_cspline, sycamera_pdist_lookup_count);
	sycamera_pdist_spec_spline1 = gsl_spline_alloc(gsl_interp_cspline, sycamera_pdist_spec_lookup_count);
	sycamera_pdist_spec_spline2 = gsl_spline_alloc(gsl_interp_cspline, sycamera_pdist_spec_lookup_count);

	gsl_spline_init(sycamera_pdist_spline1, sycamera_pdist_lookup_omega, sycamera_pdist_lookup_int1, sycamera_pdist_lookup_count);
	gsl_spline_init(sycamera_pdist_spline2, sycamera_pdist_lookup_omega, sycamera_pdist_lookup_int2, sycamera_pdist_lookup_count);
	gsl_spline_init(sycamera_pdist_spec_spline1, sycamera_pdist_spec_lookup_xi, sycamera_pdist_spec_lookup_f1, sycamera_pdist_spec_lookup_count);
	gsl_spline_init(sycamera_pdist_spec_spline2, sycamera_pdist_spec_lookup_xi, sycamera_pdist_spec_lookup_f2, sycamera_pdist_spec_lookup_count);

	sycamera_pdist_ximax = sycamera_pdist_lookup_omega[sycamera_pdist_lookup_count-1];
	sycamera_pdist_ximin = sycamera_pdist_lookup_omega[0];

	sycamera_pdist_spectrum = malloc(sizeof(double)*sycamera_pdist_spectrum_resolution);
	sycamera_pdist_polarization = malloc(sizeof(double)*4);
	sycamera_pdist_polarization_spectrum = malloc(sizeof(double*)*4);

	sycamera_pdist_polarization_spectrum[0] = malloc(sizeof(double)*sycamera_pdist_spectrum_resolution);
	sycamera_pdist_polarization_spectrum[1] = malloc(sizeof(double)*sycamera_pdist_spectrum_resolution);
	sycamera_pdist_polarization_spectrum[2] = malloc(sizeof(double)*sycamera_pdist_spectrum_resolution);
	sycamera_pdist_polarization_spectrum[3] = malloc(sizeof(double)*sycamera_pdist_spectrum_resolution);
}
void sycamera_pdist_init_particle(double mass) {
	sycamera_pdist_prefactor = 9.0 * CHARGE*CHARGE*CHARGE*CHARGE / (256.0 * PI*PI*PI * EPS0 * LIGHTSPEED * mass*mass);
	sycamera_pdist_spec_prefactor = CHARGE*CHARGE*CHARGE / (16.0 * PI * EPS0 * mass);
	sycamera_pdist_omega_B_factor = CHARGE/mass;
}

double sycamera_pdist_int(
	double gammai2, double gamma3, double gammapar2, double beta,
	double beta2, double betapar2, double Bmag, double sinmu,
	double cosmu, double sinp, double cosp, double polcosa, double polsina
) {
	double
		omegaB = sycamera_pdist_omega_B_factor*Bmag*sqrt(gammai2),
		cospsi = cosmu*cosp + sinmu*sinp,
		sinpsi2 = 1-cospsi*cospsi,
		bcospsi = beta*cospsi,
		bcospsi2 = bcospsi/2.0,
		mcospsi = 1-bcospsi,
		mcospsi2 = mcospsi*mcospsi,
		gpar2mcospsi = gammapar2*mcospsi;
	
	double pf = sycamera_pdist_prefactor * Bmag*Bmag * beta2 *
		(1 - beta*cosp*cosmu) * gammai2 / (sqrt(gpar2mcospsi*bcospsi2)*mcospsi2);
	double pf_spec = sycamera_pdist_spec_prefactor * Bmag * beta2 / (bcospsi*mcospsi);
	double mf_spec = (bcospsi*sinpsi2)/mcospsi;

	double cf = (2.0/3.0)/omegaB * sqrt(gpar2mcospsi*mcospsi2/bcospsi2),
		   xicf = cf * 2*PI*LIGHTSPEED;
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
	double I13l, I13u, I23l, I23u, I13, I23,
		   xi2K13, xi2K23, xi, lambda, Apar2, Aperp2,
		   pT, pol0, pol1, pol2, pol3,
		   polsin2a=polsina*polsina, polcos2a=polcosa*polcosa;

	sycamera_pdist_polarization[0] = 0.0;
	sycamera_pdist_polarization[1] = 0.0;
	sycamera_pdist_polarization[2] = 0.0;
	sycamera_pdist_polarization[3] = 0.0;

	/* Compute spectrum */
	int i;
	for (i = 0; i < sycamera_pdist_spectrum_resolution; i++) {
		lambda = sycamera_pdist_wavelengths[i];
		xi = xicf / lambda;
		xi2K13 = gsl_spline_eval(sycamera_pdist_spec_spline1, xi, sycamera_pdist_spec_acc1);
		xi2K23 = gsl_spline_eval(sycamera_pdist_spec_spline2, xi, sycamera_pdist_spec_acc2);

		pT = pf_spec / (lambda*lambda);
		Apar2 = pT*xi2K23;
		Aperp2 = pT*mf_spec*xi2K13;

		/* Compute polarization stuff */
		pol0 = Aperp2*polsin2a + Apar2*polcos2a;	/* |A_{left-right}|^2 */
		pol1 = Apar2*polcos2a + Apar2*polsin2a;		/* |A_{up-down}|^2 */
		pol2 = (Aperp2-Apar2)*polsina*polcosa;		/* Re(A_{up-down} x A_{left-right}*) */

		/* Im(A_{up-down} x A_{left-right}*) */
		if (Apar2*Aperp2 <= 0) pol3 = 0.0;
		else pol3 = sqrt(Apar2*Aperp2);					

		//sycamera_pdist_polarization[0] += pol0 * sycamera_pdist_dlambda;
		sycamera_pdist_polarization[0] += (Apar2 + Aperp2)*sycamera_pdist_dlambda;
		sycamera_pdist_polarization[1] += pol1 * sycamera_pdist_dlambda;
		sycamera_pdist_polarization[2] += pol2 * sycamera_pdist_dlambda;
		sycamera_pdist_polarization[3] += pol3 * sycamera_pdist_dlambda;

		sycamera_pdist_polarization_spectrum[0][i] = pol0;
		sycamera_pdist_polarization_spectrum[1][i] = pol1;
		sycamera_pdist_polarization_spectrum[2][i] = pol2;
		sycamera_pdist_polarization_spectrum[3][i] = pol3;

		sycamera_pdist_spectrum[i] = Apar2 + Aperp2;
	}

	I13l = gsl_spline_eval(sycamera_pdist_spline1, lower, sycamera_pdist_acc1);
	I13u = gsl_spline_eval(sycamera_pdist_spline1, upper, sycamera_pdist_acc1);
	I23l = gsl_spline_eval(sycamera_pdist_spline2, lower, sycamera_pdist_acc2);
	I23u = gsl_spline_eval(sycamera_pdist_spline2, upper, sycamera_pdist_acc2);

	I13 = I13l - I13u;
	I23 = I23l - I23u;

	return pf * (I23 + bcospsi2/mcospsi * (1-cospsi*cospsi) * I13);
}

double *sycamera_pdist_get_wavelengths(void) {
	return sycamera_pdist_wavelengths;
}
double *sycamera_pdist_get_spectrum(void) {
	return sycamera_pdist_spectrum;
}
int sycamera_pdist_get_spectrum_length(void) {
	return sycamera_pdist_spectrum_resolution;
}

double *sycamera_pdist_get_polarization(void) {
	return sycamera_pdist_polarization;
}
double **sycamera_pdist_get_polarization_spectrum(void) {
	return sycamera_pdist_polarization_spectrum;
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

	sycamera_pdist_init(omega0, omega1, 100);
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
			B, smu, cmu, sqrt(sin2p), sqrt(cos2p), 1., 0.
		);

		fprintf(f, "%.12e \t %.12e\n", mu, intensity);
	}

	fclose(f);
}

