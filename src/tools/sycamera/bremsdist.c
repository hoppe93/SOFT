/**
 * Angular and spectral distribution of
 * bremsstrahlung. */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "constants.h"
#include "global.h"
#include "sycamera.h"

double sycamera_bremsdist_k0, sycamera_bremsdist_k1,
	sycamera_bremsdist_dk, sycamera_bremsdist_prefactor,
	*sycamera_bremsdist_spectrum, *sycamera_bremsdist_wavelengths;
int sycamera_bremsdist_spectrum_resolution;

/* Interpolation objects */
gsl_interp_accel *sycamera_bremsdist_Zeta11_acc_x, *sycamera_bremsdist_Zeta11_acc_y,
	*sycamera_bremsdist_Zeta12_acc_x, *sycamera_bremsdist_Zeta12_acc_y,
	*sycamera_bremsdist_Zeta21_acc_x, *sycamera_bremsdist_Zeta21_acc_y;
gsl_spline2d *sycamera_bremsdist_Zeta11,
	*sycamera_bremsdist_Zeta12,
	*sycamera_bremsdist_Zeta21;

#pragma omp threadprivate(sycamera_bremsdist_wavelengths,sycamera_bremsdist_spectrum,sycamera_bremsdist_prefactor,\
	sycamera_bremsdist_Zeta11_acc_x,sycamera_bremsdist_Zeta11_acc_y,sycamera_bremsdist_Zeta12_acc_x,sycamera_bremsdist_Zeta12_acc_y,\
	sycamera_bremsdist_Zeta21_acc_x,sycamera_bremsdist_Zeta21_acc_y,sycamera_bremsdist_Zeta11,sycamera_bremsdist_Zeta12,\
	sycamera_bremsdist_Zeta21)

void sycamera_bremsdist_init(double k0, double k1, int spectrum_resolution) {
	int i;
	double dk;

	sycamera_bremsdist_k0 = k0;
	sycamera_bremsdist_k1 = k1;
	sycamera_bremsdist_spectrum_resolution = spectrum_resolution;

	if (k0 < sycamera_Zeta11_lookup_kmin) {
		fprintf(stderr, "ERROR: Lowest permitted photon energy is %e.\n", sycamera_Zeta11_lookup_kmin);
		exit(-1);
	}

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
	
	sycamera_bremsdist_Zeta11_acc_x = gsl_interp_accel_alloc();
	sycamera_bremsdist_Zeta11_acc_y = gsl_interp_accel_alloc();
	sycamera_bremsdist_Zeta12_acc_x = gsl_interp_accel_alloc();
	sycamera_bremsdist_Zeta12_acc_y = gsl_interp_accel_alloc();
	sycamera_bremsdist_Zeta21_acc_x = gsl_interp_accel_alloc();
	sycamera_bremsdist_Zeta21_acc_y = gsl_interp_accel_alloc();

	sycamera_bremsdist_Zeta11 = gsl_spline2d_alloc(gsl_interp2d_bicubic, sycamera_Zeta11_lookup_ny, sycamera_Zeta11_lookup_nx);
	sycamera_bremsdist_Zeta12 = gsl_spline2d_alloc(gsl_interp2d_bicubic, sycamera_Zeta12_lookup_ny, sycamera_Zeta12_lookup_nx);
	sycamera_bremsdist_Zeta21 = gsl_spline2d_alloc(gsl_interp2d_bicubic, sycamera_Zeta21_lookup_ny, sycamera_Zeta21_lookup_nx);

	gsl_spline2d_init(sycamera_bremsdist_Zeta11, sycamera_Zeta11_lookup_y, sycamera_Zeta11_lookup_x, sycamera_Zeta11_lookup_Z, sycamera_Zeta11_lookup_ny, sycamera_Zeta11_lookup_nx);
	gsl_spline2d_init(sycamera_bremsdist_Zeta12, sycamera_Zeta12_lookup_y, sycamera_Zeta12_lookup_x, sycamera_Zeta12_lookup_Z, sycamera_Zeta12_lookup_ny, sycamera_Zeta12_lookup_nx);
	gsl_spline2d_init(sycamera_bremsdist_Zeta21, sycamera_Zeta21_lookup_y, sycamera_Zeta21_lookup_x, sycamera_Zeta21_lookup_Z, sycamera_Zeta21_lookup_ny, sycamera_Zeta21_lookup_nx);
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

void sycamera_bremsdist_evalZ(double x, double y) {
	double Zeta11 = gsl_spline2d_eval(sycamera_bremsdist_Zeta11, y, x, sycamera_bremsdist_Zeta11_acc_y, sycamera_bremsdist_Zeta11_acc_x);
	printf("Z11(%.3f, %.3f) = %e\n", x, y, Zeta11);
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
	double E0, double E02, double beta,
	double sinmu, double cosmu, double sinp, double cosp
) {
	int i;
	double k;

	double
		cosmucosp = cosmu*cosp,
		sinmusinp = sinmu*sinp,
		p0 = E0*beta,
		p02=p0*p0;

	double E, E2, EE0, k2, Ttot, L, p, p2, p4, pp0, epsilon, s=0, T[12],
		c0, c1, c2, P1, P2, P3, xi, xi2, xi3, kappa, eta, xi_kappa,
		kkappa2, xi_kappa2, x, y, Zeta11, Zeta12, Zeta21, TT1, TT2, TT3;

	kappa = E0 - p0*cosmucosp;
	eta = p0*sinmusinp / kappa;
	xi = 1 / sqrt(1-eta*eta); xi2 = xi*xi; xi3 = xi*xi2;
	xi_kappa = xi / kappa;
	xi_kappa2 = xi_kappa*xi_kappa;
	c0 = 1 - cosmucosp*cosmucosp;
	c1 = 2*cosmucosp*sinmusinp;
	c2 = sinmusinp*sinmusinp;

	P1 = xi;
	P2 = 0.5 * (3*xi2 - 1);
	P3 = 0.5 * (5*xi3 - 3*xi);

	double
		avDelta0 = xi_kappa,
		avDelta02 = xi_kappa2 * P1,
		avQ2Delta02, avSin2aDelta04;
	
	if (sinmu > 0) {
		avSin2aDelta04 =
			c0 * xi_kappa2*xi_kappa2 * P3 -
			c1 * xi_kappa2*xi_kappa/(kappa*eta) * (xi*P3 - P2) -
			c2 * xi_kappa2/(kappa*kappa*eta*eta) * (xi2*P3 - 2*xi*P2 + P1);
	} else {
		avSin2aDelta04 = sinp*sinp / (kappa*kappa*kappa*kappa);
	}
	
	for (i = 0, k=sycamera_bremsdist_k0; i < sycamera_bremsdist_spectrum_resolution; i++, k += sycamera_bremsdist_dk) {
		E = E0-k;
		if (E <= 1) {
			sycamera_bremsdist_spectrum[i] = 0;
			continue;
		}

		k2 = k*k;
		kkappa2 = 2*k*kappa;

		E2 = E*E;
		EE0 = E*E0;
		p2 = E2 - 1;
		p4 = p2*p2;
		p = sqrt(p2);

		pp0 = p*p0;
		L = 2*log((EE0-1+pp0)/k);
		epsilon = 2*log(E0+p);

		x = eta;
		y = kkappa2 / (p2 + kkappa2);

		avQ2Delta02 = 4*k2 / (p4*sqrt((p2 + kkappa2)*(p2 + kkappa2) - (kkappa2*eta)*(kkappa2*eta))) - 2*k*xi_kappa / p4 + xi_kappa2 / p2 * P1,
		Zeta11 = gsl_spline2d_eval(sycamera_bremsdist_Zeta11, y, x, sycamera_bremsdist_Zeta11_acc_y, sycamera_bremsdist_Zeta11_acc_x);
		Zeta12 = gsl_spline2d_eval(sycamera_bremsdist_Zeta12, y, x, sycamera_bremsdist_Zeta12_acc_y, sycamera_bremsdist_Zeta12_acc_x);
		Zeta21 = gsl_spline2d_eval(sycamera_bremsdist_Zeta21, y, x, sycamera_bremsdist_Zeta21_acc_y, sycamera_bremsdist_Zeta21_acc_x);

		T[0]  = 8*(2*E02+1)/p02 * avSin2aDelta04;
		T[1]  =-2*(5*E02+2*EE0+3)/p02*avDelta02;
		T[2]  =-2*(p02-k2)*avQ2Delta02;
		T[3]  = 4*E/p02*avDelta0;
		TT1   = T[0]+T[1]+T[2]+T[3];

		T[4]  = 4*E0*(3*k - p02*E)/p02*avSin2aDelta04;
		T[5]  = 4*E02*(E02+E2)/p02*avDelta02;
		T[6]  = (2-2*(7*E02-3*EE0+E2))/p02*avDelta02;
		T[7]  = 2*k*(E02+EE0-1)/p02*avDelta0;
		TT2   = T[4]+T[5]+T[6]+T[7];

		T[8]  =-4*epsilon / p * avDelta0;

		T[9]  = (4/(p2*kappa*kappa)) * Zeta12;
		T[10] =-(6*k/(p2*kappa)) * Zeta11;
		T[11] =-2*k*(p02-k2) / (p4*kappa) * Zeta21;
		TT3   = T[9]+T[10]+T[11];

		Ttot = TT1 + TT2 + TT3 + T[8];

		Ttot = sycamera_bremsdist_prefactor * (p/k) * LIGHTSPEED*beta *
			(TT1 + (L/(p*p0))*TT2 + TT3 + T[8]);

		sycamera_bremsdist_spectrum[i] = Ttot;
		s += Ttot;
	}

	return s;
}

/* The below function is a literal implementation of Koch & Motz 2BN.
 * It has not been gyro-averaged.
 */
/*
double sycamera_bremsdist_int(
	double E0, double E02, double beta,
	double sinmu, double cosmu, double sinp, double cosp
) {
	int i, j;
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

	double E, E2, EE0, k2, Ttot, TT1, TT2, TT3, TT4, L, Q, Q2, p, pp0, epsilon, epsQ, s=0,
		T[12];

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
		
		// Awful logarithms...
		L = 2.0 * log(((EE0-1.0) + pp0) / k);
		epsilon = 2.0*log(E+p);
		epsQ = log((Q+p)*(Q+p)/(2.0*k*delta0));

		T[0] = 8.0*sinpsi2*(2.0*E02+1.0) / (p02*delta04);
		T[1] =-2.0*(5.0*E02 + 2.0*E*E0 + 3.0) / (p02 * delta02);
		T[2] =-2.0*(p02-k2) / (Q2 * delta02);
		T[3] = 4.0*E / (p02*delta0);
		TT1 = T[0] + T[1] + T[2] + T[3];

		T[4] = 4.0*E0*sinpsi2*(3.0*k - p02*E) / (p02*delta04);
		T[5] = 4.0*E02*(E02 + E2) / (p02*delta02);
		T[6] = (2.0 - 2.0*(7.0*E02 - 3.0*E*E0 + E2)) / (p02*delta02);
		T[7] = 2.0*k*(E02 + E*E0 - 1.0) / (p02*delta0);
		TT2 = T[4] + T[5] + T[6] + T[7];

		T[8] = 4.0*epsilon/(p*delta0);
		TT4 = T[8];

		T[9] = 4.0 / delta02;
		T[10] =-6.0*k/delta0;
		T[11] =2.0*k*(p02-k2) / (Q2*delta0);
		TT3 = T[9] + T[10] + T[11];

		Ttot = sycamera_bremsdist_prefactor * (p/k) * LIGHTSPEED*beta *
			(TT1 + (L/(p*p0))*TT2 + epsQ/(p*Q)*TT3 - TT4);

		sycamera_bremsdist_spectrum[i] = Ttot;
		s += Ttot;
	}

	return s*sycamera_bremsdist_dk;
}
*/

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

