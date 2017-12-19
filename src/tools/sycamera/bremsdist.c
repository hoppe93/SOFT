/**
 * Angular and spectral distribution of
 * bremsstrahlung. */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "constants.h"
#include "global.h"
#include "sycamera.h"

double sycamera_bremsdist_k0, sycamera_bremsdist_k1,
	sycamera_bremsdist_dk, sycamera_bremsdist_prefactor,
	*sycamera_bremsdist_spectrum, *sycamera_bremsdist_wavelengths;
int sycamera_bremsdist_spectrum_resolution;

#pragma omp threadprivate(sycamera_bremsdist_wavelengths,sycamera_bremsdist_spectrum,sycamera_bremsdist_prefactor)

void sycamera_bremsdist_init(double k0, double k1, int spectrum_resolution) {
	int i;
	double dk;

	sycamera_bremsdist_k0 = k0;
	sycamera_bremsdist_k1 = k1;
	sycamera_bremsdist_spectrum_resolution = spectrum_resolution;

	/* Generate photon energies */
	dk = (k1-k0) / ((double)spectrum_resolution);
	if (dk == 0) {
		sycamera_bremsdist_wavelengths = malloc(sizeof(double));
		sycamera_bremsdist_spectrum_resolution = 1;
		sycamera_bremsdist_dk = 1.0;
	} else {
		sycamera_bremsdist_wavelengths = malloc(sizeof(double)*sycamera_bremsdist_spectrum_resolution);
		sycamera_bremsdist_dk = dk;
	}

	sycamera_bremsdist_wavelengths[0] = k0;

	for (i = 0; i < sycamera_bremsdist_spectrum_resolution; i++) {
		sycamera_bremsdist_wavelengths[i] = sycamera_bremsdist_wavelengths[i-1] + dk;
	}
}

void sycamera_bremsdist_init_run(void) {
	int i;
	sycamera_bremsdist_spectrum = malloc(sizeof(double)*sycamera_bremsdist_spectrum_resolution);
	for (i = 0; i < sycamera_bremsdist_spectrum_resolution; i++)
		sycamera_bremsdist_spectrum[i] = 0;
}

/**
 * mass: Mass of particle
 * p0: Momentum of particle, normalized to m*c
 */
void sycamera_bremsdist_init_particle(double mass, double p0) {
	double r0 = CHARGE*CHARGE / (4.0*PI*EPS0*mass*LIGHTSPEED*LIGHTSPEED);	/* Electron radius */
	double alpha = CHARGE*CHARGE / (4.0*PI*EPS0*PLANCKHBAR*LIGHTSPEED);

	sycamera_bremsdist_prefactor =
		sycamera_zeff*sycamera_zeff*r0*r0*alpha / (8.0*PI*p0);
}
/**
 * Calculates the integrand for the angular & spectral
 * distribution of bremsstrahlung.
 *
 * E0: Energy of the incoming electron, normalized to mc^2.
 * E02: Squared energy of the incoming electron, normalized to m^2c^4.
 * p0: Momentum of the incoming electron, normalized to mc.
 * sinmu: Sine of angle between line-of-sight and magnetic field.
 * cosmu: Cosine of angle between line-of-sight and magnetic field.
 * sinp: Sine of pitch angle.
 * cosp: Cosine of pitch angle.
 */
double sycamera_bremsdist_int(
	/*double gammai2, double gamma3, double gammapar2, double beta,
	double beta2, double betapar2, double Bmag, double sinmu,
	double cosmu, double sinp, double cosp, double polcosa, double polsina*/
	double E0, double E02, double beta,
	double sinmu, double cosmu, double sinp, double cosp
) {
	int i;
	double k;

	double
		cospsi = cosmu*cosp + sinmu*sinp,
		sinpsi = sinmu*cosp - cosmu*sinp,
		sinpsi2= sinpsi*sinpsi,
		p0 = E0*beta,
		p02=p0*p0,
		delta0 = E0-p0*cospsi,
		delta02 = delta0*delta0,
		delta04 = delta02*delta02;

	double E, E2, EE0, k2, T, T1, T2, T3, L, Q, Q2, p, pp0, epsilon, epsQ, s=0;

	for (i = 0, k=sycamera_bremsdist_k0; i < sycamera_bremsdist_spectrum_resolution; i++, k += sycamera_bremsdist_dk) {
		E = E0-k;
		if (E <= 1) {
			sycamera_bremsdist_spectrum[i] = 0;
			continue;
		}

		E2 = E*E; k2 = k*k;
		EE0 = E*E0;
		p = sqrt(E2 - 1.0); pp0 = p*p0;
		Q2 = p02 + k2 - 2.0*p0*k*cospsi;
		Q = sqrt(Q2);
		
		/* Awful logarithms... */
		L = 2.0 * log(((EE0-1.0) + pp0) / k);
		epsilon = 2.0*log(E+p);
		epsQ = log((Q+p)*(Q+p)/(2.0*k*delta0));

		T1 = 8.0*sinpsi2*(2.0*E02+1.0) / (p02*delta04);
		T1 -= 2.0*(5.0*E02 + 2.0*E*E0 + 3.0) / (p02 * delta02);
		T1 -= 2.0*(p02-k2) / (Q2 * delta02);
		T1 += 4.0*E / (p02*delta0);

		T2 = 4.0*E0*sinpsi2*(3.0*k - p02*E) / (p02*delta04);
		T2 += 4.0*E02*(E02 + E2) / (p02*delta02);
		T2 += (2.0 - 2.0*(7.0*E02 - 3.0*E*E0 + E2)) / (p02*delta02);
		T2 += 2.0*k*(E02 + E*E0 - 1.0) / (p02*delta0);

		T3 = 4.0 / delta02;
		T3 -= 6.0*k/delta0;
		T3 -= 2.0*k*(p02-k2) / (Q2*delta0);

		T = sycamera_bremsdist_prefactor * (p/k) * LIGHTSPEED*beta *
			(T1 + (L/(p*p0))*T2 + epsQ/(p*Q)*T3 - 4.0*epsilon/(p*delta0));

		sycamera_bremsdist_spectrum[i] = T;
		s += T;
	}

	return s*sycamera_bremsdist_dk;
}

double *sycamera_bremsdist_get_wavelengths(void) {
	return sycamera_bremsdist_wavelengths;
}
double *sycamera_bremsdist_get_spectrum(void) {
	return sycamera_bremsdist_spectrum;
}
int sycamera_bremsdist_get_spectrum_length(void) {
	return sycamera_bremsdist_spectrum_resolution;
}

double *sycamera_bremsdist_get_polarization(void) { return NULL; }
double **sycamera_bremsdist_get_polarization_spectrum(void) { return NULL; }

