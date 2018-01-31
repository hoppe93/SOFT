/**
 * CONE WITH ANGULAR DISTRIBUTION
 *
 * This file implements the model in which
 * the angular distribution of radiation is
 * taken into account. This model does not
 * take the frequency into account.
 */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "global.h"
#include "sycamera.h"
#include "tools.h"

#ifdef _TESTMODE
vector *ddet, *e1, *e2, *Rdet;
double rdet;
#endif

enum sycamera_radiation_type cone_dist_radt;
int cone_dist_nint=3, cone_dist_nwavelengths;
double cone_dist_prefactor, cone_dist_preprefactor,
       cone_dist_charge, cone_dist_mass, cone_dist_speed, cone_dist_beta,
	   cone_dist_beta2, cone_dist_gammai2, cone_dist_gamma3, cone_dist_costheta,
	   cone_dist_sintheta, cone_dist_B, cone_dist_betapar2, cone_dist_betaperp2,
	   cone_dist_gammapar2, cone_dist_betapar, cone_dist_betaperp,
	   cone_dist_gamma, cone_dist_gamma2,
	   /* The polarization has four components:
		* [0] = |A_{left-right}|^2
		* [1] = |A_{up-down}|^2
		* [2] = Real component of A_{left-right} x A_{up-down}* = (Aperp^2 - Apar^2)*sina*cosa
		* [3] = Imaginary component A_{left-right} x A_{up-down}* = Aperp*Apar
		*/
	   *cone_dist_polarization, cone_dist_polcosa, cone_dist_polsina,
	   *cone_dist_wavelengths, *cone_dist_spectrum,
	   **cone_dist_polarization_spectrum, cone_dist_dwavelength;

#define NPOLARIZATION_COMPONENTS 4

double (*Ihat)(double,double,double)=NULL;

#pragma omp threadprivate(cone_dist_prefactor,cone_dist_preprefactor,\
    cone_dist_speed,cone_dist_costheta,cone_dist_sintheta,\
	cone_dist_beta, cone_dist_beta2, cone_dist_gammai2, cone_dist_gamma3,\
	cone_dist_betapar2, cone_dist_betaperp2, cone_dist_gammapar2, \
	cone_dist_betapar, cone_dist_betaperp, cone_dist_polarization, \
	cone_dist_polcosa, cone_dist_polsina, cone_dist_spectrum, \
	cone_dist_gamma, cone_dist_gamma2, cone_dist_polarization_spectrum)

double cone_dist_Ihat_benchmark(double sinmu, double cosmu, double sinmu2) {
#define CONEWIDTH 0.036
	return (fabs(sinmu-cone_dist_sintheta) <= CONEWIDTH);
}

void cone_dist_init(
	enum sycamera_radiation_type radt,
	double *wavelengths, int spectrum_resolution, int integral_resolution
) {
	cone_dist_radt = radt;

    if (radt == SYCAMERA_RADIATION_SYNCHROTRON) {
		Ihat = cone_dist_Ihat;
	} else if (radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) {
		Ihat = cone_dist_Ihat_spec;

		if (wavelengths == NULL) {
			fprintf(stderr, "ERROR: No spectral range set for the camera.\n");
			exit(-1);
		}
		if (wavelengths[0] >= wavelengths[1]) {
			fprintf(stderr, "ERROR: setting spectrum: The lower wavelength boundary must be given first.\n");
			exit(-1);
		}

		double omega0 = 2.0*PI*LIGHTSPEED / wavelengths[1];
		double omega1 = 2.0*PI*LIGHTSPEED / wavelengths[0];
		printf("omega = [%e, %e]\n", omega0, omega1);
		sycamera_pdist_init(omega0, omega1, spectrum_resolution);

		cone_dist_nwavelengths = spectrum_resolution;
		cone_dist_dwavelength = fabs((wavelengths[1]-wavelengths[0])/spectrum_resolution);
	} else if (radt == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM) {
		Ihat = cone_dist_Ihat_bremsspec;

		if (wavelengths == NULL) {
			fprintf(stderr, "ERROR: No spectral range set for the detector.\n");
			exit(-1);
		}
		if (wavelengths[0] > wavelengths[1]) {
			fprintf(stderr, "ERROR: The lower wavelength boundary must be equal or less than the upper.\n");
			exit(-1);
		}

		sycamera_bremsdist_init(wavelengths[0], wavelengths[1], spectrum_resolution);

		cone_dist_nwavelengths = spectrum_resolution;
		cone_dist_dwavelength = fabs((wavelengths[1]-wavelengths[0])/spectrum_resolution);
	} else {
        fprintf(stderr, "WARNING: This cone model only supports synchrotron radiation (w & w/o spectrum) and bremsstrahlung.\n");
    }

	cone_dist_nint = integral_resolution;
}
void cone_dist_init_run(void) {
	if (cone_dist_radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM)
		sycamera_pdist_init_run();
	else if (cone_dist_radt == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM)
		sycamera_bremsdist_init_run();

	cone_dist_spectrum = malloc(sizeof(double)*cone_dist_nwavelengths);
	cone_dist_polarization = malloc(sizeof(double)*4);
	cone_dist_polarization_spectrum = malloc(sizeof(double*)*4);
	cone_dist_polarization_spectrum[0] = malloc(sizeof(double)*cone_dist_nwavelengths);
	cone_dist_polarization_spectrum[1] = malloc(sizeof(double)*cone_dist_nwavelengths);
	cone_dist_polarization_spectrum[2] = malloc(sizeof(double)*cone_dist_nwavelengths);
	cone_dist_polarization_spectrum[3] = malloc(sizeof(double)*cone_dist_nwavelengths);
}

/**
 * Called to prepare the model for a new
 * particle being simulated.
 *
 * p: Initial parameters of the particle
 */
void cone_dist_init_particle(particle *p) {
	cone_dist_charge = p->charge;
	cone_dist_mass = p->mass;

    double e4 = cone_dist_charge*cone_dist_charge;
    e4 *= e4;

    cone_dist_preprefactor = e4 / (8.0*PI*EPS0*LIGHTSPEED*cone_dist_mass*cone_dist_mass);

	if (cone_dist_radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM)
		sycamera_pdist_init_particle(p->mass);
	else if (cone_dist_radt == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM) {
		double beta2 = (p->vpar*p->vpar + p->vperp*p->vperp) / (LIGHTSPEED*LIGHTSPEED);
		double p0 = p->mass * sqrt(beta2 / (1 - beta2));
		sycamera_bremsdist_init_particle(p->mass, p0);
	}
}
/**
 * Called to prepare the model for processing
 * of a new timestep (as taken by the ODE solver).
 *
 * sd: Information about the state parameters
 */
void cone_dist_init_step(step_data *sd) {
    cone_dist_prefactor = cone_dist_preprefactor * sd->vperp * sd->B*sd->B;
    cone_dist_speed = hypot(sd->vpar, sd->vperp);
    cone_dist_beta = cone_dist_speed / LIGHTSPEED;
	cone_dist_beta2 = cone_dist_beta*cone_dist_beta;
	cone_dist_gammai2 = 1 - cone_dist_beta2;
	cone_dist_gamma2 = 1/cone_dist_gammai2;
	cone_dist_gamma = sqrt(cone_dist_gamma2);
	cone_dist_gamma3 = cone_dist_gamma*cone_dist_gamma2;
    cone_dist_costheta = fabs(sd->vpar / cone_dist_speed);
    cone_dist_sintheta = sd->vperp / cone_dist_speed;
	cone_dist_betapar = cone_dist_beta*cone_dist_costheta;
	cone_dist_betapar2 = cone_dist_betapar*cone_dist_betapar;
	cone_dist_betaperp = cone_dist_beta*cone_dist_sintheta;
	cone_dist_betaperp2 = cone_dist_betaperp*cone_dist_betaperp;
	cone_dist_gammapar2 = 1 / (1 - cone_dist_betapar2);
	cone_dist_B = sd->B;
}

/**
 * Check whether radiation can at all hit the detector.
 * This function uses the simpler "delta cone" model.
 */
int cone_dist_can_radiation_hit(step_data *sd, vector *temp1, vector *temp2, vector *temp3) {
	vector *vhat=temp1, *rcp=temp2, *xhitv=temp3;

	vhat->val[0] = sd->vx;
	vhat->val[1] = sd->vy;
	vhat->val[2] = sd->vz;

	rcp->val[0] = sd->x - Rdet->val[0];
	rcp->val[1] = sd->y - Rdet->val[1];
	rcp->val[2] = sd->z - Rdet->val[2];

	double rcp_ddet  = vdot3(rcp, ddet);

	/* Various angles */
	double vmag = hypot(sd->vpar, sd->vperp);
	double sinThetap = sd->vperp / vmag;
	double sin2Thetap = sinThetap*sinThetap;
	double cosThetap = fabs(sd->vpar / vmag);
	double cos2Thetap= cosThetap*cosThetap;
	double cosphi = vdot3(vhat, ddet);
	double sin2phi = 1-cosphi*cosphi;
	double sinphi = sqrt(sin2phi);

	/* Rotations */
	double cosxi = vdot3(vhat, e1) / sinphi;
	double sinxi = vdot3(vhat, e2) / sinphi;
	double tanxi = sinxi / cosxi;

	/* Length of vector from particle to velocity intersection with camera plane */
	double X = fabs(rcp_ddet/cosphi);// + rg/tanThetap;

	/* Compute semi-axes */
	double cms = cos2Thetap - sin2phi;
	double a = (cosThetap*sinThetap*X*cosphi)/cms;
	double b = (sinThetap*X*cosphi)/sqrt(fabs(cms));
	double ola=-sin2Thetap/cms*X*sinphi;
	double x0=-ola*cosxi;
	xhitv->val[0] = X*vhat->val[0] + rcp->val[0];
	xhitv->val[1] = X*vhat->val[1] + rcp->val[1];
	xhitv->val[2] = X*vhat->val[2] + rcp->val[2];

	double xhit = xhitv->val[0]*ddet->val[1] + xhitv->val[1]*e1->val[1] + xhitv->val[2]*e2->val[1];
	//double yhit = xhitv->val[0]*ddet->val[2] + xhitv->val[1]*e1->val[2] + xhitv->val[2]*e2->val[2];

	/* Check if we're on the right solution */
	double t = -sinThetap*sinphi / cms * (cosThetap*cosphi + sinThetap*sinphi);
	if (cosphi < 0) {
		if (t > 1.0) return 0;
	} else if (t < 1.0 || cms > 0) {
		return 0;
	}

	double x1;

	if (cms > 0.) {	/* Ellipse */
		double tant = b/a*tanxi, tant2 = tant*tant;
		double cost = 1/sqrt(tant2+1);
		double sint = tant*cost;

		x1 = a*cost*cosxi + b*sint*sinxi;
	} else {	/* Hyperbola */
		return 1;
		/*
		double tanht = -b/a*tanxi, tanht2 = tanht*tanht;
		double cosht = 1/sqrt(1-tanht2);
		double sinht = tanht*cosht;

		if (tanht2 > 1)
			printf("a = %e, b = %e, tanht = %e, tanxi = %e\n", a, b, tanht, tanxi);

		x1 = a*cosht*cosxi + b*sinht*sinxi;
		*/
	}

	double s1 = xhit+x0+x1, s2 = xhit+x0-x1;
	double sgn1 = s1<0, sgn2 = s2<0;

	/* We make these conditions a bit more liberal...
	return ((fabs(s1)<rdet/2 ||
			 fabs(s2)<rdet/2) ||
			 (sgn1 != sgn2));
	*/
	return ((fabs(s1)<rdet ||
			 fabs(s2)<rdet) ||
			 (sgn1 != sgn2));
}
int cone_dist_can_radiation_hit2(step_data *sd, vector *temp1, vector *temp2) {
	vector *vhat=temp1, *rcp=temp2;

	vhat->val[0] = sd->vx;
	vhat->val[1] = sd->vy;
	vhat->val[2] = sd->vz;

	rcp->val[0] = sd->x - Rdet->val[0];
	rcp->val[1] = sd->y - Rdet->val[1];
	rcp->val[2] = sd->z - Rdet->val[2];

	double vdotr = vdot3(vhat, rcp);
	double absrcp2 = rcp->val[0]*rcp->val[0] + rcp->val[1]*rcp->val[1] + rcp->val[2]*rcp->val[2];
	double absrcp = sqrt(absrcp2);
	
	double lhs = cone_dist_costheta*vdotr + cone_dist_sintheta*sqrt(absrcp2 - vdotr*vdotr);

	/* Stricter condition (cone width = 1 / gamma) */
	//double rhs = absrcp * cone_dist_beta2;
	/* Less strict condition (cone_width = 2 / gamma) */
	double rhs = absrcp * (2*cone_dist_beta2 - 1);

	return (lhs > rhs);
}

double cone_dist_Ihat(double sinmu, double cosmu, double sinmu2) {
	double cmu_ctheta = cosmu * cone_dist_costheta,
		   smu_stheta = sinmu * cone_dist_sintheta,
		   kappa = 1 - cone_dist_beta*cmu_ctheta, ki = 1/kappa, ki2 = ki*ki,
		   eta = cone_dist_beta * smu_stheta * ki,
		   factor = cone_dist_preprefactor * (ki*ki2) *
		   		cone_dist_betaperp2 * cone_dist_B*cone_dist_B * cone_dist_gammai2 *
				cone_dist_gammapar2 * (1 - cone_dist_betapar*cosmu);
				//cone_dist_gammapar2 * (1+cone_dist_beta*cosmu);
	
	if (eta == 0.0) return factor;

	double xi = 1/sqrt(1-eta*eta),
		   xi2 = xi*xi,
		   xi3 = xi2*xi,
		   xi5 = xi3*xi2,
		   xi7 = xi5*xi2,
		   p1, p2, p2f;
	
	p1 = 3.0/2.0 * xi5 - 1.0/2.0 * xi3;
	p2 = 5.0/8.0 * xi7 - 1.0/8.0 * xi5;

	p2f = sinmu2 * cone_dist_gammai2 * ki2;

	return factor * (p1 - p2f * p2);
}
double cone_dist_Ihat_spec(double sinmu, double cosmu, double sinmu2) {
	return sycamera_pdist_int(
		cone_dist_gammai2, cone_dist_gamma3, cone_dist_gammapar2,
		cone_dist_beta, cone_dist_beta2, cone_dist_betapar2,
		cone_dist_B, sinmu, cosmu, cone_dist_sintheta, cone_dist_costheta,
		cone_dist_polcosa, cone_dist_polsina
	);
}
double cone_dist_Ihat_bremsspec(double sinmu, double cosmu, double sinmu2) {
	return sycamera_bremsdist_int(
		cone_dist_gamma, cone_dist_gamma2, cone_dist_beta, sinmu, cosmu, cone_dist_sintheta, cone_dist_costheta
	);
}

/**
 * Computes sines and cosines of the angle mu
 *
 * dX: (X - X0). Offset from detector center in e1 direction
 * dY: (Y - Y0). Offset from detector center in e2 direction
 * rcp: Vector from camera to particle
 * vhat: Unit vector in direction of guiding-center motion
 * cosmu: Contains cos(mu) on return
 * sinmu: Contains sin(mu) on return
 */
void cone_dist_get_angles(
    double dX, double dY,
    vector *rcp, vector *vhat,
    double *cosmu, double *sinmu,
	double *sinmu2
) {
    double rcpx = rcp->val[0] + dX*e1->val[0] + dY*e2->val[0],
           rcpy = rcp->val[1] + dX*e1->val[1] + dY*e2->val[1],
           rcpz = rcp->val[2] + dX*e1->val[2] + dY*e2->val[2];
    double r = sqrt(rcpx*rcpx + rcpy*rcpy + rcpz*rcpz);

	double
		bhat_ehat = vdot3(vhat, e2),
		bhat_rcp  = (vhat->val[0]*rcpx + vhat->val[1]*rcpy + vhat->val[2]*rcpz) / r,
		ehat_rcp  = (e2->val[0]*rcpx   + e2->val[1]*rcpy   + e2->val[2]*rcpz) / r,
		bhatXrcpe = (e2->val[0] * (vhat->val[1]*rcp->val[2] - vhat->val[2]*rcp->val[1]) +
					 e2->val[1] * (vhat->val[2]*rcp->val[0] - vhat->val[0]*rcp->val[2]) +
					 e2->val[2] * (vhat->val[0]*rcp->val[1] - vhat->val[1]*rcp->val[0])) / r,
		divfac = 1.0/sqrt((1-bhat_rcp*bhat_rcp)*(1-ehat_rcp*ehat_rcp));

	cone_dist_polcosa = (bhat_ehat - bhat_rcp*ehat_rcp) * divfac;
	cone_dist_polsina = bhatXrcpe * divfac;

    *cosmu = -(rcpx*vhat->val[0] + rcpy*vhat->val[1] + rcpz*vhat->val[2]) / r;
    *sinmu2= 1 - (*cosmu)*(*cosmu);
    *sinmu = sqrt(*sinmu2);
}

/**
 * Carry out the integral over the Y coordinate.
 *
 * X: Is really (X-X0), the distance from the center
 *    of the detector!
 * dY: Stepsize in Y
 */
double cone_dist_integrateY(
    double X, double dY,
    vector *rcp, vector *vhat,
	double factor	/* Numerical integration factor * geometric factor */
) {
    int i;
    double Y, cosmu, sinmu, sinmu2, dY2, s;

	factor *= cone_dist_dwavelength;
    dY2 = dY+dY;

    /* Simpson's rule */
    cone_dist_get_angles(X, -rdet, rcp, vhat, &cosmu, &sinmu, &sinmu2);
    s = (*Ihat)(sinmu, cosmu, sinmu2);
	cone_dist_add_spectrum_and_polarization(factor);

    cone_dist_get_angles(X, +rdet, rcp, vhat, &cosmu, &sinmu, &sinmu2);
    s += (*Ihat)(sinmu, cosmu, sinmu2);
	cone_dist_add_spectrum_and_polarization(factor);

    for (i = 1, Y = -rdet + dY; i < cone_dist_nint; i += 2, Y += dY2) {
        cone_dist_get_angles(X, Y, rcp, vhat, &cosmu, &sinmu, &sinmu2);

        s += 4.0 * (*Ihat)(sinmu, cosmu, sinmu2);
		cone_dist_add_spectrum_and_polarization(factor);
    }
    for (i = 2, Y = -rdet + dY2; i < cone_dist_nint-1; i += 2, Y += dY2) {
        cone_dist_get_angles(X, Y, rcp, vhat, &cosmu, &sinmu, &sinmu2);

        s += 2.0 * (*Ihat)(sinmu, cosmu, sinmu2);
		cone_dist_add_spectrum_and_polarization(factor);
    }

    return s;
}

/**
 * Computes the intensity from one particle
 * in this model.
 */
double cone_dist_get_intensity(
	step_data *sd, vector *rcp, vector *vhat, vector *empty1,
	vector *empty2, vector *empty3
) {
    int i;
    double X, dX, dX2, s,
		r2 = rcp->val[0]*rcp->val[0] + rcp->val[1]*rcp->val[1] + rcp->val[2]*rcp->val[2],
		weight = vdot3(rcp, ddet) / sqrt(r2), fac;

    //if (!cone_dist_can_radiation_hit2(sd, empty1, empty2))
	if (!cone_dist_can_radiation_hit(sd, empty1, empty2, empty3))
		return 0.0;
	
    dX = 2.0 * rdet / ((double)(cone_dist_nint - 1));
    dX2 = dX+dX;
	fac = dX*dX*weight / (9.0*r2);

	cone_dist_reset_spectrum();
	cone_dist_reset_polarization();

    /* Simpson's rule */
    s = cone_dist_integrateY(-rdet, dX, rcp, vhat, fac) +
        cone_dist_integrateY(rdet, dX, rcp, vhat, fac);

    for (i = 1, X = -rdet + dX; i < cone_dist_nint; i += 2, X += dX2) {
        s += 4.0 * cone_dist_integrateY(X, dX, rcp, vhat, 4.0*fac);
    }
    for (i = 2, X = -rdet + dX2; i < cone_dist_nint-1; i += 2, X += dX2) {
        s += 2.0 * cone_dist_integrateY(X, dX, rcp, vhat, 2.0*fac);
    }

    /* We multiply twice by the step and divide
       twice by 3 to compensate for not doing it
       in the Y integral. */
    return s * fac;
}

/**
 * Spectrum and polarization acquisition routines
 */
double *cone_dist_get_wavelengths(void) {
	if (cone_dist_radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) {
		return sycamera_pdist_get_wavelengths();
	} else if (cone_dist_radt == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM) {
		return sycamera_bremsdist_get_wavelengths();
	} else return NULL;
}
double *cone_dist_get_spectrum(void) {
	if (cone_dist_radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM ||
		cone_dist_radt == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM) {
		return cone_dist_spectrum;
	} else return NULL;
}
double *cone_dist_get_polarization(void) {
	if (cone_dist_radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) {
		return cone_dist_polarization;
	} else return NULL;
}
double **cone_dist_get_polarization_spectrum(void) {
	if (cone_dist_radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) {
		return cone_dist_polarization_spectrum;
	} else return NULL;
}
int cone_dist_get_spectrum_length(void) {
	if (cone_dist_radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) {
		return cone_dist_nwavelengths;
	} else return 0;
}

/**
 * Spectrum and polarization addition/resets.
 */
void cone_dist_add_spectrum_and_polarization(double factor) {
	if (cone_dist_radt != SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) return;

	double *spectrum, *polarization, **polspec;
	spectrum = sycamera_pdist_get_spectrum();
	polarization = sycamera_pdist_get_polarization();
	polspec = sycamera_pdist_get_polarization_spectrum();

	if (spectrum == NULL) return;

	cone_dist_add_spectrum(spectrum, factor);
	cone_dist_add_polarization(polarization, factor);
	cone_dist_add_polarization_spectrum(polspec, factor);
}
void cone_dist_add_spectrum(double *spec, double factor) {
	int i;
	for (i = 0; i < cone_dist_nwavelengths; i++)
		cone_dist_spectrum[i] += factor * spec[i];
}
void cone_dist_add_polarization(double *polarization, double factor) {
	cone_dist_polarization[0] += polarization[0] * factor;
	cone_dist_polarization[1] += polarization[1] * factor;
	cone_dist_polarization[2] += polarization[2] * factor;
	cone_dist_polarization[3] += polarization[3] * factor;
}
void cone_dist_add_polarization_spectrum(double **spectrum, double factor) {
	int i;
	for (i = 0; i < cone_dist_nwavelengths; i++) {
		cone_dist_polarization_spectrum[0][i] += spectrum[0][i] * factor;
		cone_dist_polarization_spectrum[1][i] += spectrum[1][i] * factor;
		cone_dist_polarization_spectrum[2][i] += spectrum[2][i] * factor;
		cone_dist_polarization_spectrum[3][i] += spectrum[3][i] * factor;
	}
}

void cone_dist_reset_spectrum(void) {
	int i;
	for (i = 0; i < cone_dist_nwavelengths; i++)
		cone_dist_spectrum[i] = 0.0;
}
void cone_dist_reset_polarization(void) {
	cone_dist_polarization[0] = 0.0;
	cone_dist_polarization[1] = 0.0;
	cone_dist_polarization[2] = 0.0;
	cone_dist_polarization[3] = 0.0;
}
void cone_dist_reset_polarization_spectrum(void) {
	int i;
	for (i = 0; i < cone_dist_nwavelengths; i++) {
		cone_dist_polarization_spectrum[0][i] = 0.0;
		cone_dist_polarization_spectrum[1][i] = 0.0;
		cone_dist_polarization_spectrum[2][i] = 0.0;
		cone_dist_polarization_spectrum[3][i] = 0.0;
	}
}

void cone_dist_test(void) {
	rdet = 0.003;
	ddet = vinit(3, 0, 1, 0);
	e1 = vinit(3, 1, 0, 0);
	e2 = vinit(3, 0, 0, 1);
	Rdet = vinit(3, 0, -1, 0);
#define eV2kgms 5.36e-28
	double PPAR = 3e7 * eV2kgms;
	double PPERP = 5e6 * eV2kgms;
	
	particle p = {0,0,9.10938356e-31,1.60217662e-19,0,0,NULL,NULL,0.0,0,0,0,0,0.0};
	step_data sd;
	sd.ppar2 = PPAR*PPAR;
	sd.pperp2 = PPERP*PPERP;

	double emass = 9.10938356e-31;
	double gamma2 = 1 + (sd.ppar2 + sd.pperp2) / (emass*emass*LIGHTSPEED*LIGHTSPEED);

	sd.B = 4;
	sd.x = sd.y = sd.z = 0;
	sd.vx = sd.vy = sd.vz = 0;
	sd.vpar = PPAR / (sqrt(gamma2)*emass);
	sd.vperp = PPERP / (sqrt(gamma2)*emass);

	cone_dist_init_particle(&p);
	cone_dist_init_step(&sd);

	FILE *f = fopen("angdist.out", "w");
	if (!f) {
		perror("ERROR");
		exit(-1);
	}

	double I, mu, dmu;
	int i, n = 4000;
	dmu = PI/2.0 / (n-1);
	for (i = 0, mu = 0; i < n; i++, mu+=dmu) {
		I = cone_dist_Ihat(sin(mu), cos(mu), sin(mu)*sin(mu));
		fprintf(f, "%.12e \t %.12e\n", mu, I);
	}

	fclose(f);
}
