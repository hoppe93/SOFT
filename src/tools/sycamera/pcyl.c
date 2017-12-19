/* Calculation of the synchrotron power emitted in a given
   wavelength range. The function Pcyl was originally given
   by [Bekefi, 1966] and implemented in MATLAB by A. Stahl.
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <omp.h>
#include "magnetic_field.h"
#include "sycamera.h"
#include "vector.h"

double *sycamera_pcyl_lambdas, sycamera_pcyl_dl, *sycamera_pcyl_comp;
int sycamera_pcyl_lambda_resolution, sycamera_pcyl_haswarned;
gsl_interp_accel *sycamera_pcyl_acc, *sycamera_pcyl_pol_acc;
gsl_spline *sycamera_pcyl_spline, *sycamera_pcyl_pol_spline;
/* Stokes parameters (polariation components):
 * [0] = I
 * [1] = Q
 * [2] = U
 * [3] = V
 */
double *sycamera_pcyl_polarization, **sycamera_pcyl_polarization_spectrum;

#pragma omp threadprivate(sycamera_pcyl_haswarned,sycamera_pcyl_acc,\
	sycamera_pcyl_pol_acc, sycamera_pcyl_pol_spline,\
	sycamera_pcyl_spline,sycamera_pcyl_comp,\
	sycamera_pcyl_polarization,sycamera_pcyl_polarization_spectrum)

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
	
	/* Initialize polarization */
	sycamera_pcyl_polarization = malloc(sizeof(double)*4);
	sycamera_pcyl_polarization_spectrum = malloc(sizeof(double*)*4);
	sycamera_pcyl_polarization_spectrum[0] = malloc(sizeof(double)*sycamera_pcyl_lambda_resolution);
	sycamera_pcyl_polarization_spectrum[1] = malloc(sizeof(double)*sycamera_pcyl_lambda_resolution);
	sycamera_pcyl_polarization_spectrum[2] = malloc(sizeof(double)*sycamera_pcyl_lambda_resolution);
	sycamera_pcyl_polarization_spectrum[3] = malloc(sizeof(double)*sycamera_pcyl_lambda_resolution);

	/* Initialize interpolation */
    /* Initialize power spline */
    sycamera_pcyl_acc = gsl_interp_accel_alloc();
    /* Linear interpolation (somewhat faster, less accurate) */
    //sycamera_pcyl_spline = gsl_spline_alloc(gsl_interp_linear, sycamera_pcyl_lookup_count);
    /* Cubic spline interpolation */
    sycamera_pcyl_spline = gsl_spline_alloc(gsl_interp_cspline, sycamera_pcyl_lookup_count);
    gsl_spline_init(sycamera_pcyl_spline, sycamera_pcyl_lookup_lambda, sycamera_pcyl_lookup_int, sycamera_pcyl_lookup_count);

	/* Initialize polarization splines */
	sycamera_pcyl_pol_acc = gsl_interp_accel_alloc();
	sycamera_pcyl_pol_spline = gsl_spline_alloc(gsl_interp_cspline, sycamera_pcyl_lookup_pol_count);
	gsl_spline_init(sycamera_pcyl_pol_spline, sycamera_pcyl_lookup_pol_lambda, sycamera_pcyl_lookup_pol_int, sycamera_pcyl_lookup_pol_count);
}
void sycamera_pcyl_init_step(step_data *sd) { }

/**
 * Computes the integral of P_cyl from lambda0 to lambda1
 *
 * lambda0: The lower wavelength to evaluate the integral at
 * lambda1: The upper wavelength to evaluate the integral at
 * p: Particle momentum
 * xi: Cosine of the pitch angle
 * Bmag: Magnetic field strength in the point of evaluation
 */
double sycamera_pcyl_int(double ppar2, double pperp2, double Bmag, double mass, double fraction, vector *rcp, vector *vhat) {
    double ic = 1/((LIGHTSPEED*LIGHTSPEED)*mass*mass);
    double ppar2l = ppar2*ic;
    double pperp2l = pperp2*ic;
    double gamma2 = 1+ppar2l+pperp2l;
	double gammaParallel2 = gamma2 / (1+pperp2l);
    double gammaParallel = sqrt(gamma2 / (1+pperp2l));
	double betaparpar = ppar2l / sqrt(gamma2*gamma2 - gamma2);	// = beta*cospitch*cospitch
    double lambdac = SYCAMERA_PCYL_LAMBDAC_CONST * gammaParallel * mass / (Bmag*gamma2);
    double prefactor = SYCAMERA_PCYL_CONSTANT * gammaParallel2  * (1 - betaparpar) / gamma2;
	double prefactorpol = SYCAMERA_PCYL_POL_CONSTANT * gamma2*gamma2 * Bmag*Bmag*Bmag * gammaParallel2 * (1-betaparpar) / (gammaParallel*gammaParallel2);
	double Apar2, Aperp2;
	
	/* Compute polarization parameters */
	double ir = 1/vnorm3(rcp),
		   ne = vdot3(rcp, e2)*ir,
		   nb = vdot3(rcp, vhat)*ir,
		   ve = vdot3(vhat, e2),
		   ne2 = ne*ne,
		   nb2 = nb*nb,
		   divfac = 1.0 / sqrt((1-nb2)*(1-ne2)),
		   cosb = (ve - nb*ne) * divfac,
		   sinb = (
			   e2->val[0] * (vhat->val[1]*rcp->val[2] - vhat->val[2]*rcp->val[1]) +
			   e2->val[1] * (vhat->val[2]*rcp->val[0] - vhat->val[0]*rcp->val[2]) +
			   e2->val[2] * (vhat->val[0]*rcp->val[1] - vhat->val[1]*rcp->val[0])
		   )*divfac*ir;
	
	sycamera_pcyl_polarization[0] = 0.0;
	sycamera_pcyl_polarization[1] = 0.0;
	sycamera_pcyl_polarization[2] = 0.0;
	sycamera_pcyl_polarization[3] = 0.0;

    int i;
    double sum=0., ilambda, lb, pf,
		   pol1=0., pol2=0., pol3=0., pol4=0.;
    for (i = 0; i < sycamera_pcyl_lambda_resolution; i++) {
        ilambda = 1 / sycamera_pcyl_lambdas[i];
        lb = lambdac * ilambda;
        pf = prefactor*ilambda*ilambda*ilambda;

        /* Lambdas outside the bounds are assumed negligible */
        if (lb >= sycamera_pcyl_lookup_lambda[sycamera_pcyl_lookup_count-1]) {
			sycamera_pcyl_comp[i] = 0;
            continue;
        } else if (lb <= sycamera_pcyl_lookup_lambda[0]) {
            if (!sycamera_pcyl_haswarned) {
                fprintf(stderr, "WARNING: The detector spectrum used is outside the supported bounds. Please adjust the lower bound to > %e.\n", lambdac/(sycamera_pcyl_lookup_lambda[0]));
                sycamera_pcyl_haswarned = 1;
            }
			sycamera_pcyl_comp[i] = 0;
            continue;
        }

		/* Compute spectrum */
		double e = gsl_spline_eval(sycamera_pcyl_spline, lb, sycamera_pcyl_acc);
		double c = pf * e * fraction;
		sycamera_pcyl_comp[i] = c;
        sum += c;

		/* Compute polarization */
		double k = gsl_spline_eval(sycamera_pcyl_pol_spline, lb, sycamera_pcyl_pol_acc);

		Apar2  = 0.5*(c + prefactorpol*k*fraction);
		Aperp2 = 0.5*(c - prefactorpol*k*fraction);

		pol1 = Apar2 + Aperp2;
		pol2 = (Aperp2 - Apar2) * (sinb*sinb - cosb*cosb);
		pol3 = 2.0 * (Aperp2 - Apar2) * sinb*cosb;
		pol4 = 0;

		sycamera_pcyl_polarization_spectrum[0][i] = pol1;
		sycamera_pcyl_polarization_spectrum[1][i] = pol2;
		sycamera_pcyl_polarization_spectrum[2][i] = pol3;
		sycamera_pcyl_polarization_spectrum[3][i] = pol4;

		sycamera_pcyl_polarization[0] += pol1;
		sycamera_pcyl_polarization[1] += pol2;
		sycamera_pcyl_polarization[2] += pol3;
		sycamera_pcyl_polarization[3] += pol4;
    }

	sycamera_pcyl_polarization[0] *= sycamera_pcyl_dl;
	sycamera_pcyl_polarization[1] *= sycamera_pcyl_dl;
	sycamera_pcyl_polarization[2] *= sycamera_pcyl_dl;
	sycamera_pcyl_polarization[3] *= sycamera_pcyl_dl;

    return sum * sycamera_pcyl_dl;
}

double *sycamera_pcyl_get_wavelengths(void) { return sycamera_pcyl_lambdas; }
double *sycamera_pcyl_get_spectrum(void) { return sycamera_pcyl_comp; }
int sycamera_pcyl_get_spectrum_length(void) { return sycamera_pcyl_lambda_resolution; }
double *sycamera_pcyl_get_polarization(void) { return sycamera_pcyl_polarization; }
double **sycamera_pcyl_get_polarization_spectrum(void) { return sycamera_pcyl_polarization_spectrum; }

