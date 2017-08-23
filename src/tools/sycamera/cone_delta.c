/**
 * DELTA CONE MODEL
 *
 * This file implements the "delta cone model", whereby
 * the cone is modeled as having an infinitesimally thin
 * wall and opening angle = the pitch angle. The necessary
 * equations have been derived fully for the case of one
 * particle, and so the computation is rather direct.
 */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <gsl/gsl_sf_dilog.h>
#include "diag.h"
#include "global.h"
#include "sycamera.h"
#include "tools.h"

/* Projection operator */
double Pn[3][3], Pnx[3], Pne1r[3], Pne2r[3], absPne1r, absPne2r;
double cone_delta_charge, cone_delta_mass;
enum sycamera_radiation_type cone_delta_radiation_type=SYCAMERA_RADIATION_SYNCHROTRON;

#pragma omp threadprivate(Pn,Pnx,Pne1r,Pne2r,absPne1r,absPne2r)

void cone_delta_init(
	enum sycamera_radiation_type radt, enum sycamera_polarization_type polt,
	double *lambdas, int spectrum_resolution, int integral_resolution
) {
	cone_delta_radiation_type = radt;

	/* Initialize spectrum recognition */
	if (radt == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) {
		if (lambdas == NULL) {
			fprintf(stderr, "ERROR: No spectral range set for the camera.\n");
			exit(-1);
		}
		if (lambdas[0] >= lambdas[1]) {
			fprintf(stderr, "ERROR: setting spectrum: The lower wavelength boundary must be given first.\n");
			exit(-1);
		}
		sycamera_spectrum_init(lambdas[0], lambdas[1], spectrum_resolution);
	} else if (radt == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM) {
		if (lambdas == NULL) {
			fprintf(stderr, "ERROR: No spectral range set for the detector.\n");
			exit(-1);
		}
		if (lambdas[0] >= lambdas[1]) {
			fprintf(stderr, "ERROR: setting spectrum: The lower photon energy must be given first.\n");
			exit(-1);
		}
		sycamera_bsspec_init(lambdas[0], lambdas[1], spectrum_resolution);
	}

	if (polt != SYCAMERA_POLARIZATION_BOTH) {
		fprintf(stderr, "ERROR: The cone model has no support for polarized radiation.\n");
		exit(-1);
	}
}
void cone_delta_init_run(void) {
	if (cone_delta_radiation_type == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM)
		sycamera_spectrum_init_run();
	else if (cone_delta_radiation_type == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM)
		sycamera_bsspec_init_run();
}
void cone_delta_init_particle(particle *p) {
	cone_delta_charge = p->charge;
	cone_delta_mass = p->mass;
}
void cone_delta_init_step(step_data *sd) { }

/**
 * Computes the intensity detected from one
 * guiding-center in the given position.
 */
double cone_delta_get_intensity(
	step_data *sd, vector *rcp, vector *vhat, vector *empty1,
	vector *empty2, vector *empty3
) {
    double fraction = cone_delta_radiation_hits(sd, rcp, vhat, empty1, empty2, empty3);
    double totalP = 0.0;

    if (fraction > 0.0) {
    	if (cone_delta_radiation_type == SYCAMERA_RADIATION_BREMSSTRAHLUNG) {
			double p02 = (sd->pperp2+sd->ppar2) / (cone_delta_mass*cone_delta_mass*LIGHTSPEED*LIGHTSPEED);
			double p0 = sqrt(p02);
    		double E02 = 1+p02;
			double E0 = sqrt(E02);
			double E0p0 = E0*p0;
			double lnE0p0 = log(E0+p0);
			double Fx = gsl_sf_dilog(2*p0*(E0+p0));

			totalP = (12*E02 + 4) / (3*E0p0) * lnE0p0;
			totalP-= (8*E0 + 6*p0) / (3*E0p0*p0) * lnE0p0*lnE0p0 - 4.0/3.0;
			totalP+= 2/E0p0 * Fx;

    		//totalP = gamma*(log(2*gamma) - 1./3.);
			totalP *= fraction;
		} else if (cone_delta_radiation_type == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM) {
			totalP = sycamera_bsspec_int(sd->ppar2, sd->pperp2, cone_delta_mass, fraction);
    	} else if (cone_delta_radiation_type == SYCAMERA_RADIATION_CONSTANT) {
    		totalP = fraction;
    	} else if (cone_delta_radiation_type == SYCAMERA_RADIATION_SYNCHROTRON) {
    		double q2 = cone_delta_charge*cone_delta_charge;
			double betapar = sd->vpar / LIGHTSPEED;
			double cospitch = sd->vpar / hypot(sd->vpar, sd->vperp);
			//double betaperp = sd->vperp / LIGHTSPEED;
			double gammapar2 = 1/(1-betapar*betapar);
			double c3 = LIGHTSPEED*LIGHTSPEED*LIGHTSPEED;
			double m2 = cone_delta_mass*cone_delta_mass, m4 = m2*m2;
    		totalP = q2*q2*sd->B*sd->B*sd->pperp2*gammapar2*(1-betapar*cospitch)/(6*PI*EPS0*c3*m4);

			totalP *= fraction;
    	} else if (cone_delta_radiation_type == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM) {
    		totalP = sycamera_spectrum_weight(sd, cone_delta_mass, fraction);
    	} else {
    		fprintf(stderr, "ERROR: Unrecognized radiation type selected: %d\n", cone_delta_radiation_type);
    		exit(1);
    	}
    }

    return totalP;
}

/**
 * "Quick" check to see if radiation can at all hit
 * the camera.
 */
int cone_delta_can_radiation_hit(step_data *sd, vector *temps) {
	vector *vhat=temps, *rcp=temps+1, *xhitv=temps+2;

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
	double cosThetap = sd->vpar / vmag;
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
	double a = fabs((cosThetap*sinThetap*X*cosphi)/cms);
	double b = (sinThetap*X*cosphi)/sqrt(fabs(cms));
	double ola=-sin2Thetap/cms*X*sinphi;
	double x0=-ola*cosxi;
	xhitv->val[0] = X*vhat->val[0] + rcp->val[0];
	xhitv->val[1] = X*vhat->val[1] + rcp->val[1];
	xhitv->val[2] = X*vhat->val[2] + rcp->val[2];

	double xhit = xhitv->val[0]*ddet->val[1] + xhitv->val[1]*e1->val[1] + xhitv->val[2]*e2->val[1];
	//double yhit = xhitv->val[0]*ddet->val[2] + xhitv->val[1]*e1->val[2] + xhitv->val[2]*e2->val[2];
	
	/* Check if we're on the right solution */
	if (cosphi < 0) {
		if ((a-ola)*sinphi > X) return 0.0;
	} else if ((a-ola)*sinphi < X || cms > 0) return 0.0;

	double x1;

	if (cms > 0.) {	/* Ellipse */
		double tant = b/a*tanxi, tant2 = tant*tant;
		double cost = 1/sqrt(tant2+1);
		double sint = tant*cost;

		x1 = a*cost*cosxi + b*sint*sinxi;
	} else {	/* Hyperbola */
		double tanht = -b/a*tanxi, tanht2 = tanht*tanht;
		double cosht = 1/sqrt(1-tanht2);
		double sinht = tanht*cosht;

		x1 = a*cosht*cosxi + b*sinht*sinxi;
	}

	double s1 = xhit+x0+x1, s2 = xhit+x0-x1;
	double sgn1 = s1<0, sgn2 = s2<0;

	return ((fabs(s1)<rdet/2 ||
			 fabs(s2)<rdet/2) ||
			 (sgn1 != sgn2));
}

/* Calculates the intersection points of the
 * radiation ellipse.
 *
 * A & B are picked from the expression
 *   x (y) = x_0 (y_0) +(-) A*cos(t) + B*sin(t)
 * origin is x_0 (or y_0, depending on equation to solve)
 * t1 and t2 are the intersection points, which are returned
 * by this function
 * ok1 and ok2 are flags which are set by this function. A
 * 'true' value means the intersection point corresponding
 * to the flag is valid.
 *
 * The y-vals corresponding to t1 and t2 are returned
 * in the pointers 'y1', 'y2' respectively
 */
void cone_delta_ellipse_intersections(
	double Ax, double Bx, double Ay, double By, double d,
	double *t1, double *t2, double *y1, double *y2,
	int *ok1, int *ok2,
	double origin, double half_side_length
) {
	double dif = Ax-d;
	if (dif == 0.) {
		*ok1 = 0; *ok2 = 0;
		return;
	}

	/* Find intersection points */
	double fac = Ax*Ax+Bx*Bx-d*d;
	if (fac < 0.) {
		*ok1 = 0; *ok2 = 0;
		return;
	}

	double a=sqrt(fac)/dif, b=Bx/dif;
	double fac1 = b+a, fac12=fac1*fac1;
	double fac2 = b-a, fac22=fac2*fac2;
	double sqr1 = 1/(fac12+1);
	double sqr2 = 1/(fac22+1);
	*t1 = 2*atan(fac1);
	*t2 = 2*atan(fac2);

	double cost1 = (1-fac12)*sqr1, cost2 = (1-fac22)*sqr2;
	double sint1 = 2*fac1*sqr1, sint2 = 2*fac2*sqr2;

	/* Make sure they also satisfy the second condition */
	*y1 = origin - Ay*cost1 + By*sint1;
	if (-half_side_length <= *y1 && *y1 <= half_side_length)
		*ok1 = 1;
	else *ok1 = 0;

	*y2 = origin - Ay*cost2 + By*sint2;
	if (-half_side_length <= *y2 && *y2 <= half_side_length)
		*ok2 = 1;
	else *ok2 = 0;
}

/**
 * This function is similar to the one above. All arguments
 * are the hyperbolic equivalents, and the return values are
 * the intersections of the hyperbola with the camera.
 */
void cone_delta_hyperbola_intersections(
	double Ax, double Bx, double Ay, double By, double d,
	double *t1, double *t2, double *y1, double *y2,
	int *ok1, int *ok2,
	double origin, double half_side_length
) {
	if (Ax+Bx == 0.0) { *ok1 = 0; *ok2 = 0; return; }

	/* Calculate intersections */
	double u = d*d - Ax*Ax + Bx*Bx;
	double b = 1/(Ax+Bx);
	/* If u < 0, the solutions are imaginary <=> no solutions */
	if (u < 0) { *ok1 = 0; *ok2 = 0; return; }
	double sqr = sqrt(u);
	double arg1 = (sqr-d)*b;
	double arg2 = (-sqr-d)*b;

	//if (arg1 <= 0.0) arg1 = -arg1;
	/* Only allow real solutions */
	if (arg1 > 0) {
		*t1 = log(arg1);

		double cosht1 = cosh(*t1), sinht1 = sinh(*t1);
		/* Make sure it also satisfies the second condition */
		*y1 = origin - Ay*cosht1 + By*sinht1;
		if (-half_side_length <= *y1 && *y1 <= half_side_length)
			*ok1 = 1;
		else *ok1 = 0;
	} else *ok1 = 0;

	//if (arg2 <= 0.0) arg2 = -arg2;
	/* Only allow real solutions */
	if (arg2 > 0) {
		*t2 = log(arg2);

		double cosht2 = cosh(*t2), sinht2 = sinh(*t2);
		/* Make sure they also satisfy the second condition */
		*y2 = origin - Ay*cosht2 + By*sinht2;
		if (-half_side_length <= *y2 && *y2 <= half_side_length)
			*ok2 = 1;
		else *ok2 = 0;
	} else *ok2 = 0;
}

/**
 * Find intersections between projected curve and camera.
 *
 * t: Array to store intersection points (8 elements)
 * xval: Array to store corresponding x values (8 elements)
 * yval: Array to store corresponding y values (8 elements)
 * tindex: Returned as number of intersection points found
 * x0: Curve offset in x from center of camera
 * y0: Curve offset in y from center of camera
 * a: Major semi-axis
 * b: Minor semi-axis
 * rdet2: Half detector-side (rdet/2)
 * intersect: Proper intersection function to use
 */
void cone_delta_find_intersections(
	double *t, double *xval, double *yval, int *tindex,
	double x0, double y0, double a, double b, double rdet2,
	double cosxi, double sinxi,
	void (*intersect)(
		double, double, double, double, double,
		double*, double*, double*, double*,
		int*, int*, double, double
	)
) {
	double t1,t2,x1,x2,y1,y2;
	int ok1=0,ok2=0;

	/************************************************/
	/*********** FIND INTERSECTION POINTS ***********/
	/************************************************/
	/*         x = rd/2, -rd/2 < y < +rd/2          */
	x1 = x2 = rdet2;
	(*intersect)(
		a*cosxi, b*sinxi, a*sinxi, b*cosxi, x0-rdet2,
		&t1, &t2, &y1, &y2,
		&ok1, &ok2, y0, rdet2
	);
	if(ok2){t[(*tindex)]=t2;xval[(*tindex)]=x2;yval[(*tindex)]=y2;(*tindex)++;}
	if(ok1){t[(*tindex)]=t1;xval[(*tindex)]=x1;yval[(*tindex)]=y1;(*tindex)++;}

	/************************************************/
	/*         x = -rd/2, -rd/2 < y < +rd/2          */
	x1 = x2 = -rdet2;
	(*intersect)(
		a*cosxi, b*sinxi, a*sinxi, b*cosxi, x0+rdet2,
		&t1, &t2, &y1, &y2,
		&ok1, &ok2, y0, rdet2
	);
	if(ok2){t[(*tindex)]=t2;xval[(*tindex)]=x2;yval[(*tindex)]=y2;(*tindex)++;}
	if(ok1){t[(*tindex)]=t1;xval[(*tindex)]=x1;yval[(*tindex)]=y1;(*tindex)++;}

	/************************************************/
	/*         y = rd/2, -rd/2 < x < +rd/2          */
	y1 = y2 = rdet2;
	(*intersect)(
		-a*sinxi, b*cosxi, -a*cosxi, b*sinxi, y0-rdet2,
		&t1, &t2, &x1, &x2,
		&ok1, &ok2, x0, rdet2
	);
	if(ok2){t[(*tindex)]=t2;xval[(*tindex)]=x2;yval[(*tindex)]=y2;(*tindex)++;}
	if(ok1){t[(*tindex)]=t1;xval[(*tindex)]=x1;yval[(*tindex)]=y1;(*tindex)++;}

	/************************************************/
	/*         y = -rd/2, -rd/2 < x < +rd/2         */
	y1 = y2 = -rdet2;
	(*intersect)(
		-a*sinxi, b*cosxi, -a*cosxi, b*sinxi, y0+rdet2,
		&t1, &t2, &x1, &x2,
		&ok1, &ok2, x0, rdet2
	);
	if(ok2){t[(*tindex)]=t2;xval[(*tindex)]=x2;yval[(*tindex)]=y2;(*tindex)++;}
	if(ok1){t[(*tindex)]=t1;xval[(*tindex)]=x1;yval[(*tindex)]=y1;(*tindex)++;}
}

/**
 * Prepare for computing the ellipse overlap
 */
void cone_delta_oncone_init(vector *empty, vector *X, double Xmag, double sinxi, double cosxi) {
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			if (i == j) Pn[i][j] = 1;
			else Pn[i][j] = 0;

			Pn[i][j] -= X->val[i]*X->val[j]/(Xmag*Xmag);
		}
	}

	vector *e1r = empty;
	e1r->val[0] = e1->val[0]*cosxi + e2->val[0]*sinxi;
	e1r->val[1] = e1->val[1]*cosxi + e2->val[1]*sinxi;
	e1r->val[2] = e1->val[2]*cosxi + e2->val[2]*sinxi;

	Pne1r[0] = Pn[0][0]*e1r->val[0] + Pn[0][1]*e1r->val[1] + Pn[0][2]*e1r->val[2];
	Pne1r[1] = Pn[1][0]*e1r->val[0] + Pn[1][1]*e1r->val[1] + Pn[1][2]*e1r->val[2];
	Pne1r[2] = Pn[2][0]*e1r->val[0] + Pn[2][1]*e1r->val[1] + Pn[2][2]*e1r->val[2];

	vector *e2r = empty;
	e2r->val[0] = -e1->val[0]*sinxi + e2->val[0]*cosxi;
	e2r->val[1] = -e1->val[1]*sinxi + e2->val[1]*cosxi;
	e2r->val[2] = -e1->val[2]*sinxi + e2->val[2]*cosxi;

	Pne2r[0] = Pn[0][0]*e2r->val[0] + Pn[0][1]*e2r->val[1] + Pn[0][2]*e2r->val[2];
	Pne2r[1] = Pn[1][0]*e2r->val[0] + Pn[1][1]*e2r->val[1] + Pn[1][2]*e2r->val[2];
	Pne2r[2] = Pn[2][0]*e2r->val[0] + Pn[2][1]*e2r->val[1] + Pn[2][2]*e2r->val[2];

	absPne1r = sqrt(Pne1r[0]*Pne1r[0] + Pne1r[1]*Pne1r[1] + Pne1r[2]*Pne1r[2]);
	absPne2r = sqrt(Pne2r[0]*Pne2r[0] + Pne2r[1]*Pne2r[1] + Pne2r[2]*Pne2r[2]);
}

/**
 * Convert a point on the ellipse/hyperbola to a point on the unit cone
 *
 * empty: An empty vector which this function can use freely
 * X:     Vector from particle to point in plane which is intersected by particle's velocity vector
 * Xmag:  Magnitude of vector X
 * x:     X coordinate (in the system which has inclination along the x-axis only)
 * y:     Y coordinate (in the system which has inclination along the x-axis only)
 * sinxi:  Sine of the rotation angle between the system described two & three lines above and n-e1-e2
 * cosxi:  Cosine of the rotation angle between the system described on the two lines above and n-e1-e2
 */
double cone_delta_oncone(
	vector *empty, double x, double y
) {
	vector *r = empty;
	r->val[0] = x*e1->val[0] + y*e2->val[0];
	r->val[1] = x*e1->val[1] + y*e2->val[1];
	r->val[2] = x*e1->val[2] + y*e2->val[2];

	Pnx[0] = Pn[0][0]*r->val[0] + Pn[0][1]*r->val[1] + Pn[0][2]*r->val[2];
	Pnx[1] = Pn[1][0]*r->val[0] + Pn[1][1]*r->val[1] + Pn[1][2]*r->val[2];
	Pnx[2] = Pn[2][0]*r->val[0] + Pn[2][1]*r->val[1] + Pn[2][2]*r->val[2];

	double absPnx = sqrt(Pnx[0]*Pnx[0] + Pnx[1]*Pnx[1] + Pnx[2]*Pnx[2]);

	double cosDelta = (Pnx[0]*Pne1r[0] + Pnx[1]*Pne1r[1] + Pnx[2]*Pne1r[2])/(absPnx*absPne1r);
	double sinDelta = (Pnx[0]*Pne2r[0] + Pnx[1]*Pne2r[1] + Pnx[2]*Pne2r[2])/(absPnx*absPne2r);

	if (sinDelta > 0) {
		return acos(cosDelta);
	} else {
		return 2*PI-acos(cosDelta);
	}
}

/**
 * Checks whether the emitted radiation hits the detector, and if
 * so how much of it that hits.
 *
 * sd: Step data structure
 * rcp: Vector from camera to particle
 * vhat: Vector pointing in the direction of guiding-center motion
 * empty1, 2 & 3: Empty pre-allocated vectors
 * rcp_ddet: The dot product between rcp and camera normal
 *
 * Returns the fraction of emitted power that hits the detector
 * (a value between 0 and 1).
 */
double cone_delta_radiation_hits(
	step_data *sd, vector *rcp, vector *vhat, vector *empty1,
	vector *empty2, vector *empty3
) {
	/* Calculate some helper numbers (useful when doing projections) */
	double vmag = hypot(sd->vpar, sd->vperp);
	double sinThetap = fabs(sd->vperp) / vmag;
	double sin2Thetap = sinThetap*sinThetap;
	double cosThetap = fabs(sd->vpar) / vmag;
	double cos2Thetap= cosThetap*cosThetap;
	double cosphi = vhat->val[0]*ddet->val[0] + vhat->val[1]*ddet->val[1] + vhat->val[2]*ddet->val[2];
	double sin2phi = 1-cosphi*cosphi;
	double sinphi = sqrt(sin2phi);		/* sinphi should always be > 0 by definition */

	/* Rotate parametrization */
	double cosxi=1, sinxi=0;
	if (sinphi != 0.0) {
		cosxi = vdot3(vhat, e1)/sinphi;
		sinxi = vdot3(vhat, e2)/sinphi;
	}

	/* Calculate length of cone focus point vector, projected in the camera normal */
	double X = fabs(vdot3(rcp, ddet)/cosphi);
	vector *Xv = empty1;
	Xv->val[0] = X*vhat->val[0];
	Xv->val[1] = X*vhat->val[1];
	Xv->val[2] = X*vhat->val[2];

	double cms = cos2Thetap - sin2phi;
	//double a = (cosThetap*sinThetap*X)/cms;
	double a = fabs((cosThetap*sinThetap*X*cosphi)/cms);
	double b = fabs(sinThetap*X*cosphi)/sqrt(fabs(cms));

	/* Now we have two cases. For cms > 0, we get an ellipse,
	 * while for cms < 0 we get a hyperbola. The way we find the
	 * intersections and calculate the total overlap is analogous
	 * for both objects, but with a slightly different formula
	 * for finding the intersection points.
	 */
	void (*intersect)(
		double, double, double, double, double,
		double*, double*, double*, double*,
		int*, int*, double, double
	);

	if(cms>0) {	/* Ellipse */
		intersect = &cone_delta_ellipse_intersections;
	} else { 	/* Hyperbola */
		intersect = &cone_delta_hyperbola_intersections;
	}

	/* Calculate offset from center of detector */
	double ola=-sin2Thetap/cms*X*sinphi;
	/* These are (x0*cosxi + x_hit) etc. NOT just x0 as written in the docs (ola = x0 of the docs) */
	/*
	double x0=ola*cosxi - (rd->val[0]*ddet->val[1] + rd->val[1]*e1->val[1] + rd->val[2]*e2->val[1]);
	double y0=-ola*sinxi - (rd->val[0]*ddet->val[2] + rd->val[1]*e1->val[2] + rd->val[2]*e2->val[2]);
	*/

	vector *xhit = empty3;
	xhit->val[0] = Xv->val[0] + rcp->val[0];
	xhit->val[1] = Xv->val[1] + rcp->val[1];
	xhit->val[2] = Xv->val[2] + rcp->val[2];
	double xh = xhit->val[0]*ddet->val[1] + xhit->val[1]*e1->val[1] + xhit->val[2]*e2->val[1];
	double yh = xhit->val[0]*ddet->val[2] + xhit->val[1]*e1->val[2] + xhit->val[2]*e2->val[2];

	double x0 = -ola*cosxi + xh;
	double y0 = ola*sinxi + yh;

	/* Check if we're on the right solution (forward cone) */
	if (cosphi < 0) {
		if ((a-ola)*sinphi > X) return 0.0;
	} else if ((a-ola)*sinphi < X || cms > 0) return 0.0;

	double t[8], xval[8], yval[8], rdet2=rdet/2;
	int tindex=0;

	/* If the pitch angle is 0, all radiation hits in
	 * a single point. */
	if (sinThetap == 0.) {
		if (fabs(x0)<rdet2 && fabs(y0)<rdet2) {
			/* Full hit! */
			return 1.0;
		} else return 0.0;	/* No hit */
	}

	/* Find intersection points */
	cone_delta_find_intersections(t, xval, yval, &tindex, x0, y0, a, b, rdet2, cosxi, sinxi, intersect);

	/* No solutions found! This means that the radiation
	 * either hits entirely inside the camera, or entirely
	 * outside. */
	if (tindex < 2) {
		if (cms > 0) {	/* If ellipse */
			/* If any point on the ellipse is inside the square,
			 * then so must all others be. Same for the opposite case */
			double x = x0 + a*cosxi;
			double y = y0 - a*sinxi;
			if (fabs(x)<rdet2 && fabs(y)<rdet2) {
				/* Full hit! */
				return 1.0;
			} else return 0.0;	/* No hit */
		} else return 0.0;	/* No solutions for hyperbola => no hit */
	}

	/***********************************************/
	/*     COMPUTE CORRESPONDING POINT ON CONE     */
	/***********************************************/
	int i, swapped=1, tlen=tindex;
	double swp;
	/* Bubble sort the intersection points in order smallest to largest */
	/* (array will be like [0 -> 2pi] or [-inf -> inf]) */
	while (swapped) {
		swapped = 0;
		for (i = 1; i < tlen; i++) {
			if (t[i-1] > t[i]) {
				swp = t[i];
				t[i] = t[i-1];
				t[i-1] = swp;

				swp = xval[i];
				xval[i] = xval[i-1];
				xval[i-1] = swp;

				swp = yval[i];
				yval[i] = yval[i-1];
				yval[i-1] = swp;

				swapped = 1;
			}
		}
		tlen--;
	}

	/* A point is only either a 'to'-point or a 'from'-point.
	 * 'to'-points bring positive contributions to the total,
	 * while 'from'-points bring negative contributions.
	 * Determine which is which
	 */
	double sgn = -1;
	double p = (t[1]+t[0])/2, lastpoint=0.0;
	double fraction = 0.0, contrib;
	int count = tindex;

	/* Prepare for computing angles */
	cone_delta_oncone_init(empty2, Xv, X, sinxi, cosxi);

	/* QUICK REMINDER
	 * 'From'-points are points you INTEGRATE from.
	 * 'To'-points are points you INTEGRATE to.
	 */

	/* If an ellipse... */
	if (cms > 0) {
		/* If the first point is a 'to'-point, change sign */
		if (fabs(x0+a*cosxi*cos(p)+b*sinxi*sin(p)) >= rdet2||
			fabs(y0-a*sinxi*cos(p)+b*cosxi*sin(p)) >= rdet2) {
			sgn = 1;

			/* Get its corresponding 'from'-point */
			//double xp = a*cosxi*cos(t[tindex-1]) + b*sinxi*sin(t[tindex-1]);
			//double yp =-a*sinxi*cos(t[tindex-1]) + b*cosxi*sin(t[tindex-1]);
			double xp = xval[tindex-1] - x0;
			double yp = yval[tindex-1] - y0;
			lastpoint = cone_delta_oncone(empty2, xp, yp);
			/* No need to handle the final point again */
			count--;
			fraction += -lastpoint;
		}
	} /* Note: the first point of a hyperbola MUST be a 'from'-point,
	     since its domain is (-inf, inf) */

	for (i = 0; i < count; i++) {
		/* Calculate contribution */
		//double xp = a*cosxi*cos(t[i]) + b*sinxi*sin(t[i]);
		//double yp =-a*sinxi*cos(t[i]) + b*cosxi*sin(t[i]);
		double xp = xval[i] - x0;
		double yp = yval[i] - y0;
		contrib = cone_delta_oncone(empty2, xp, yp);

		fraction += sgn*contrib;
		/* If this is a "to"-point and this point is smaller
		 * than the previous point, the line segment connecting
		 * them passes through 0. Therefore, add a term of 2pi */
		if (sgn > 0 && contrib < lastpoint) {
			fraction += 2*PI;
		}

		lastpoint = contrib;
		sgn = -1*sgn;
	}

	fraction /= 2*PI;
	return fraction;
}

double *cone_delta_get_wavelengths(void) {
	if (cone_delta_radiation_type == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM)
		return sycamera_spectrum_get_wavelengths();
	else if (cone_delta_radiation_type == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM)
		return sycamera_bsspec_get_wavelengths();
	else return NULL;
}
double *cone_delta_get_spectrum(void) {
	if (cone_delta_radiation_type == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM)
		return sycamera_spectrum_get();
	else if (cone_delta_radiation_type == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM)
		return sycamera_bsspec_get_spectrum();
	else return NULL;
}
int cone_delta_get_spectrum_length(void) {
	if (cone_delta_radiation_type == SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM)
		return sycamera_spectrum_length();
	else if (cone_delta_radiation_type == SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM)
		return sycamera_bsspec_get_spectrum_length();
	else return 0;
}

