/* Interpolation of magnetic_field object
 * input given in cylindrical coordinates
 */
#include <omp.h>
#include "magfield.h"
#include "magnetic_field.h"
#include "magnetic_num.h"
#include "vector.h"

/* GNU Scientific Library */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_version.h>

/* Since GSL 2.0 'interp2d' has been merged
 * with GSL. For older versions (which I actually
 * use) we must use the independent release of
 * 'interp2d' for bicubic interpolation.
 */
#if GSL_MAJOR_VERSION < 2
#	include <interp2d.h>
#	include <interp2d_spline.h>
#	define gsl_spline2d interp2d_spline
#	define gsl_spline2d_init interp2d_spline_init
#	define gsl_spline2d_alloc interp2d_spline_alloc
#	define gsl_interp2d_bicubic interp2d_bicubic
#	define gsl_spline2d_eval interp2d_spline_eval
#	define gsl_spline2d_eval_deriv_x interp2d_spline_eval_deriv_x
#	define gsl_spline2d_eval_deriv_y interp2d_spline_eval_deriv_y
#else
#	include <gsl/gsl_interp2d.h>
#	include <gsl/gsl_spline2d.h>
#endif

#include <math.h>

/* This value is used in case the magnetic field
 * is probed outside its area of definition to guarantee
 * a value.
 */
#define EPS 1e-7

/* GSL accelerators that help speed up interpolations */
gsl_interp_accel *ra, *za;
/* interp2d interpolation objects */
gsl_spline2d *Br;
gsl_spline2d *Bphi;
gsl_spline2d *Bz;
vector *B_interp;

double rmin, rmax, zmin, zmax;
double **jacobian;

#pragma omp threadprivate(rmin,rmax,zmin,zmax,jacobian,B_interp,ra,za,Br,Bphi,Bz,B_interp)

/**
 * Initializatlize magnetic_field for interpolation.
 * This is to prepare GSL for what's to come.
 * Must be called before interpolation!
 */
void interp2_init_interpolation(magfield_t *B) {
	/* Prepare the accelerators for both r and z */
	ra = gsl_interp_accel_alloc();
	za = gsl_interp_accel_alloc();

	/* Create interpolation objects */
	Br = gsl_spline2d_alloc(gsl_interp2d_bicubic, B->nz, B->nr);
	Bphi = gsl_spline2d_alloc(gsl_interp2d_bicubic, B->nz, B->nr);
	Bz = gsl_spline2d_alloc(gsl_interp2d_bicubic, B->nz, B->nr);

	/* Prepare the interpolation objects for our situation */
	gsl_spline2d_init(Br, B->z, B->r, B->Br[0], B->nz, B->nr);
	gsl_spline2d_init(Bphi, B->z, B->r, B->Bphi[0], B->nz, B->nr);
	gsl_spline2d_init(Bz, B->z, B->r, B->Bz[0], B->nz, B->nr);

	rmin = B->r[0];
	rmax = B->r[B->nr-1];
	zmin = B->z[0];
	zmax = B->z[B->nz-1];

	B_interp = vnew(3);
	jacobian = malloc(3*sizeof(double*));
	jacobian[0] = malloc(3*sizeof(double));
	jacobian[1] = malloc(3*sizeof(double));
	jacobian[2] = malloc(3*sizeof(double));
}

/**
 * main interpolation function
 *
 * xyz: The point (in cartesian coordinates) in which the field
 * strength should be evaluated.
 */
vector* interp2_interpolate(double x, double y, double z) {
	/* Transform vector coordinates from cartesian to cylindrical */
	double  r = sqrt(x*x + y*y);

	/* Make sure the point is within our defined area */
	if (r < rmin) r = rmin+EPS;
	else if (r > rmax) r = rmax-EPS;

	if (z < zmin) z = zmin+EPS;
	else if (z > zmax) z = zmax-EPS;

	/* Interpolate */
	double Br_interp, Bphi_interp, Bz_interp;
	Br_interp = gsl_spline2d_eval(Br, z, r, za, ra);
	Bphi_interp = gsl_spline2d_eval(Bphi, z, r, za, ra);
	Bz_interp = gsl_spline2d_eval(Bz, z, r, za, ra);

	if (gsl_isnan(Br_interp)) Br_interp = 0;
	if (gsl_isnan(Bphi_interp)) Bphi_interp = 0;
	if (gsl_isnan(Bz_interp)) Bz_interp = 0;

	/* Transform field to cartesian coordinates */
	double sinp = y/r, cosp = x/r;
	double Bx_interp = Br_interp * cosp - Bphi_interp * sinp;
	double By_interp = Br_interp * sinp + Bphi_interp * cosp;

	/* Store interpolation values in vector */
	B_interp->val[0] = Bx_interp;
	B_interp->val[1] = By_interp;
	B_interp->val[2] = Bz_interp;

	return B_interp;
}

/**
 * Calculate the Jacobian of the B-field
 * at the given point
 *
 * xyz: The point (in cartesian coordinates) in which the field
 * strength should be evaluated.
 */
double **interp2_jacobian(double x, double y, double z) {
	/* Transform vector coordinates from cartesian to cylindrical */
	double  r = sqrt(x*x + y*y);

	/* Make sure the point is within our definied are */
	if (r < rmin) r = rmin+EPS;
	else if (r > rmax) r = rmax-EPS;

	if (z < zmin) z = zmin+EPS;
	else if (z > zmax) z = zmax-EPS;

	double cylJ[3][2] = {{0.,0.}, {0.,0.}, {0.,0.}};
	/* dBr/dr */
	cylJ[0][0] = gsl_spline2d_eval_deriv_y(Br, z, r, za, ra);
	/* dBr/dz */
	cylJ[0][1] = gsl_spline2d_eval_deriv_x(Br, z, r, za, ra);
	/* dB0/dr */
	cylJ[1][0] = gsl_spline2d_eval_deriv_y(Bphi, z, r, za, ra);
	/* dB0/dz */
	cylJ[1][1] = gsl_spline2d_eval_deriv_x(Bphi, z, r, za, ra);
	/* dBz/dr */
	cylJ[2][0] = gsl_spline2d_eval_deriv_y(Bz, z, r, za, ra);
	/* dBz/dz */
	cylJ[2][1] = gsl_spline2d_eval_deriv_x(Bz, z, r, za, ra);

	double
		sin0 = y/r, cos0 = x/r, sin20 = sin0*sin0, cos20 = cos0*cos0, sc0 = sin0*cos0,
		dsin0_dx = -sc0/r,  dcos0_dx = sin20/r,
		dsin0_dy = cos20/r, dcos0_dy = -sc0/r,
		/* Partial derivatives of cylindrical components */
		dBr_dr = cylJ[0][0], dBr_dz = cylJ[0][1],
		dB0_dr = cylJ[1][0], dB0_dz = cylJ[1][1],
		dBz_dr = cylJ[2][0], dBz_dz = cylJ[2][1],
		/* Make B components more readable */
		_Br = gsl_spline2d_eval(Br, z, r, za, ra),
		_B0 = gsl_spline2d_eval(Bphi, z, r, za, ra);

	/* Bx */
	jacobian[0][0] = cos20*dBr_dr + _Br*dcos0_dx - sc0  *dB0_dr - _B0*dsin0_dx;
	jacobian[0][1] = sc0  *dBr_dr + _Br*dcos0_dy - sin20*dB0_dr - _B0*dsin0_dy;
	jacobian[0][2] = cos0 *dBr_dz - sin0 *dB0_dz;
	/* By */
	jacobian[1][0] = sc0  *dBr_dr + _Br*dsin0_dx + cos20*dB0_dr + _B0*dcos0_dx;
	jacobian[1][1] = sin20*dBr_dr + _Br*dsin0_dy + sc0  *dB0_dr + _B0*dcos0_dy;
	jacobian[1][2] = sin0 *dBr_dz + cos0 *dB0_dz;
	/* Bz */
	jacobian[2][0] = cos0*dBz_dr;
	jacobian[2][1] = sin0*dBz_dr;
	jacobian[2][2] = dBz_dz;

	return jacobian;
}
