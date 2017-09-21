/**
 * Computation of the Bremsstrahlung spectrum.
 */

#include <math.h>
#include <omp.h>
#include "constants.h"
#include "sycamera.h"

double *sycamera_bsspec_photen, sycamera_bsspec_dk, *sycamera_bsspec_comp;
double sycamera_bsspec_prefac;
int sycamera_bsspec_lambda_resolution;

#pragma omp threadprivate(sycamera_bsspec_comp)

void sycamera_bsspec_init(double lambda0, double lambda1, int res) {
	int i;
	double r0, alpha;

	sycamera_bsspec_photen = malloc(sizeof(double)*res);
	sycamera_bsspec_lambda_resolution = res;
	sycamera_bsspec_dk = (lambda1-lambda0) / (res-1);

	for (i = 0; i < res; i++) {
		sycamera_bsspec_photen[res-i-1] = lambda0 + i * sycamera_bsspec_dk;
	}

	sycamera_bsspec_dk = fabs(sycamera_bsspec_dk);

	r0 = CHARGE*CHARGE/(ELECTRONMASS*LIGHTSPEED*LIGHTSPEED);
	alpha = CHARGE*CHARGE / (2.0*EPS0*PLANCKH*LIGHTSPEED);
	sycamera_bsspec_prefac = 2.0 * sycamera_zeff * r0*r0 * alpha;
}

void sycamera_bsspec_init_run(void) {
	int i;
	/* Initialize spectrum storage */
	sycamera_bsspec_comp = malloc(sizeof(double)*sycamera_bsspec_lambda_resolution);
	for (i = 0; i < sycamera_bsspec_lambda_resolution; i++) {
		sycamera_bsspec_comp[i] = 0;
	}
}

/**
 * Computes the differential cross-section
 * in photon energy for bremsstrahlung.
 *
 * E0: Incident electron energy
 * p0: Incident electron momentum
 * p02: Incident electron momentum squared
 * eps0 = 2*ln((E0 + p0) / (E0-p0)) ~ 2*ln(E0+p0)
 * k: Emitted photon energy
 */
double sycamera_bsspec_dsigma(double E0, double p0, double p02, double eps0, double k) {
	double
		E = E0-k,
		p2 = E*E-1,
		p = sqrt(p2),
		p3 = p2*p,
		p03 = p02*p0,
		E0E = E0*E,
		p0p = p0*p,
		p0p2 = p0p*p0p,
		p0p3 = p0p2*p0p,
		pf = p / (p0*k),
		eps = 2*log(E+p),
		L = 2*log((E0E + p0p - 1) / k),
		T1 = (4.0/3.0) - 2*E0E*(p02 + p2) / p0p2 +
			eps0*E/p03 + eps*E0/p3 - eps*eps0/p0p,
		T2 = L * (8*E0E/(3*p0p) + k*k/p0p3 * (E0E*E0E + p0p2)),
		T3 = L * k/(2*p0p) * ((E0E+p02)*eps0/p03 -
			(E0E+p2)*eps/p3 + 2*k*E0E/p0p2);
	
	printf("%e\n", pf * (T1+T2+T3));
	return pf * (T1+T2+T3);
}
/**
 * Computes the integral of the bremsstrahlung spectrum from
 * lambda0 to lambda1.
 */
double sycamera_bsspec_int(double ppar2, double pperp2, double mass, double fraction) {
#define mc (9.10938356e-31*LIGHTSPEED)
#define hc_mc2 ((PLANCKH*LIGHTSPEED) / (mc*LIGHTSPEED))

	double p02 = (ppar2+pperp2)/(mc*mc), E0 = sqrt(1+p02), p0 = sqrt(p02), k,
		eps0 = 2*log(E0+p0), sum=0, c, v;
	
	int i;
	for (i = 0; i < sycamera_bsspec_lambda_resolution; i++) {
		k = sycamera_bsspec_photen[i];

		if (k < E0-1) {
			c = fraction * sycamera_bsspec_dsigma(E0, p0, p02, eps0, k);
			sycamera_bsspec_comp[i] = c;
			sum += c;
		} else {
			sycamera_bsspec_comp[i] = 0;
		}
	}


	v = p0 / (E0*mass);	/* Electron speed */
	return sycamera_bsspec_prefac * sum * sycamera_bsspec_dk * v;
}

double *sycamera_bsspec_get_wavelengths(void) { return sycamera_bsspec_photen; }
double *sycamera_bsspec_get_spectrum(void) { return sycamera_bsspec_comp; }
int sycamera_bsspec_get_spectrum_length(void) { return sycamera_bsspec_lambda_resolution; }
void sycamera_bsspec_test(void) {
	double p0 = 12.0, k = 5.0;
	double p02 = p0*p0, E0 = sqrt(1+p02),
		   eps0 = 2*log(E0+p0);

	double r = sycamera_bsspec_dsigma(E0, p0, p02, eps0, k);
	printf("sigma(p0=%.2f, k=%.2f) = %e\n", p0, k, r);
}

