/* Distribution function interpolation */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_version.h>
#include <math.h>
#include "distfunc.h"

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

gsl_interp_accel **df_interp_accelerators_t,
				 **df_interp_accelerators_p;
gsl_spline2d **df_interp_splines;
double df_interp_dr,
	   df_interp_rmin, df_interp_cosmin, df_interp_pmin,
	   df_interp_rmax, df_interp_cosmax, df_interp_pmax;
distfunc *df_interp_f;

#pragma omp threadprivate(df_interp_accelerators_t,df_interp_accelerators_p,df_interp_splines)

void df_interp_init(distfunc *f) {
	gsl_set_error_handler(df_interp_error_handler);
	df_interp_dr = 0;
	df_interp_f = f;

	df_interp_rmin = f->rmin, df_interp_rmax = f->rmax;
	df_interp_cosmin = f->ximin, df_interp_cosmax = f->ximax;
	df_interp_pmin = f->pmin, df_interp_pmax = f->pmax;

	if (f->nr > 1) df_interp_dr = (f->rmax-f->rmin)/(f->nr-1);
}
void df_interp_init_run(void) {
	df_interp_accelerators_t = (gsl_interp_accel**)malloc(sizeof(gsl_interp_accel*)*df_interp_f->nr);
	df_interp_accelerators_p = (gsl_interp_accel**)malloc(sizeof(gsl_interp_accel*)*df_interp_f->nr);
	df_interp_splines = (gsl_spline2d**)malloc(sizeof(gsl_spline*)*df_interp_f->nr);

	size_t r;
	for (r = 0; r < df_interp_f->nr; r++) {
		df_interp_accelerators_t[r] = gsl_interp_accel_alloc();
		df_interp_accelerators_p[r] = gsl_interp_accel_alloc();
		df_interp_splines[r] = gsl_spline2d_alloc(gsl_interp2d_bicubic, df_interp_f->np, df_interp_f->nxi);

		gsl_spline2d_init(df_interp_splines[r], df_interp_f->p, df_interp_f->xi, df_interp_f->value[r], df_interp_f->np, df_interp_f->nxi);
	}
}

double df_interp_eval(double r, double xi, double p) {
	/* Check input */
	if (r < df_interp_rmin || r > df_interp_rmax) {
		fprintf(stderr, "FATAL: Particle located outside distribution function! r = %e, rmin = %e, rmax = %e\n", r, df_interp_rmin, df_interp_rmax);
		exit(EXIT_FAILURE);
	}
	if (xi < df_interp_cosmin || xi > df_interp_cosmax) {
		fprintf(stderr, "FATAL: Particle pitch angle is not within distribution function! cos(theta) = %e\n", xi);
		exit(EXIT_FAILURE);
	}
	if (p < df_interp_pmin || p > df_interp_pmax) {
		fprintf(stderr, "FATAL: Particle momentum is not within distribution function! p = %e, pmin = %e, pmax = %e\n", p, df_interp_pmin, df_interp_pmax);
		exit(EXIT_FAILURE);
	}

	/* Find closest matches for r and xi */
	int ir = (int)(round((r - df_interp_f->rmin)/df_interp_dr));

	if (ir >= (signed int)df_interp_f->nr || ir < 0) {
		fprintf(stderr, "WARNING: Attempting to evaluate distribution function outside it's bounds! (r=%e)", r);
		ir = 0;
	}

	double v = gsl_spline2d_eval(df_interp_splines[ir], p, xi, df_interp_accelerators_t[ir], df_interp_accelerators_p[ir]);
	//printf("p = %e, theta = %e, f = %e\n", p, acos(xi), v);
	//printf("p = %e, pitch = %e,  f = %e\n", p, acos(xi), v);
	return fabs(v);
}

void df_interp_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
	switch (gsl_errno) {
		case GSL_EDOM:
			fprintf(stderr, "FATAL: Particle defined in distribution function is not within particle map!\n");
			break;
		default:
			fprintf(stderr, "FATAL: %s:%d: %s\n", file, line, reason);
			break;
	}

	exit(EXIT_FAILURE);
}

