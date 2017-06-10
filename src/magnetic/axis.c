/* Code for finding the magnetic axis
 * of a magnetic field
 *
 * This function is used to determine where the magnetic
 * axis is, in order to drop particles aligned with it
 * in the z direction.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "magnetic_axis.h"
#include "magnetic_field.h"
#include "util.h"

double magnetic_axis_r, magnetic_axis_z;
int magnetic_axis_set=0;

/**
 * Compute the magnetic axis of the current
 * magnetic field. Initial guess in z is always
 * 0.
 *
 * rguess: Initial guess in r.
 */
double *magnetic_axis_find(double rguess) {
	double x1=0., x2=0., f1=0., f2=0.,
		   a=0., b=.1, gradr, gradz,
		   rold = 0., zold = 0.,
		   rnew = 0., znew = 0.;
	double tau = (sqrt(5)-1)/2;	/* Golden ratio */
	double tol = 1e-5;	/* Tolerance (should be the same for both step length and descent) */
	vector *v;

	rnew = rguess;	/* Start at rguess, we always choose z = 0 */

	diff_data *dd = magnetic_field_diff_notor(rnew, 0., znew);

	do {
		rold = rnew, zold = znew;

		gradr = dd->gradB->val[0];
		gradz = dd->gradB->val[2];

		/* Prepare for finding step length */
		a = 0., b = 1.;
		x1 = a + (1-tau)*(b-a);
		x2 = a + tau*(b-a);

		v = magnetic_field_get(rold - x1*gradr, 0., zold - x1*gradz);
		f1 = hypot(v->val[0], v->val[2]);

		v = magnetic_field_get(rold - x2*gradr, 0., zold - x2*gradz);
		f2 = hypot(v->val[0], v->val[2]);

		/* Find most appropriate step length */
		while ((b-a) > tol) {
			if (f1 > f2) {
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + tau * (b-a);
				rnew = rold - x2*gradr;
				znew = zold - x2*gradz;
				v = magnetic_field_get(rnew, 0., znew);
				f2 = hypot(v->val[0], hypot(v->val[1], v->val[2]));
			} else {
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = a + (1-tau)*(b-a);
				rnew = rold - x1*gradr;
				znew = zold - x1*gradz;
				v = magnetic_field_get(rnew, 0., znew);
				f1 = hypot(v->val[0], v->val[2]);
			}
		}

		//printf("(r,z) = %e, %e;   step = %e\n", rnew, znew, x1);

		dd = magnetic_field_diff_notor(rnew, 0., znew);
	} while (fabs(rold-rnew) > tol || fabs(zold-znew) > tol);

	double *retval = malloc(sizeof(double)*2);
	magnetic_axis_set_loc(rnew, znew);
	retval[0] = rnew;
	retval[1] = znew;

	return retval;
}
void magnetic_axis_set_loc(double r, double z) {
	magnetic_axis_r = r;
	magnetic_axis_z = z;
	magnetic_axis_set = 1;

	printf("Magnetic axis found at (%e, %e)\n", r, z);
}

void magnetic_axis_test(void) {
	magnetic_init();
	magnetic_handler *circular = magnetic_handler_select("circular");
	struct general_settings *mset_circular = malloc(sizeof(struct general_settings));
	mset_circular->name = "circular";
	mset_circular->n = 3;
	mset_circular->setting = malloc(mset_circular->n*sizeof(char*));
	mset_circular->value = malloc(mset_circular->n*sizeof(char*));
	mset_circular->setting[0] = setname("B0");
	mset_circular->value[0]   = setname("6");
	mset_circular->setting[1] = setname("major_radius");
	mset_circular->value[1]   = setname("6.2");
	mset_circular->setting[2] = setname("safety_factor");
	mset_circular->value[2]   = setname("1");

	circular->init(mset_circular);
	circular->init_run();
	circular->init_particle();

	double *rz = magnetic_axis_find(6.);
	printf("Magnetic axes:\n");
	printf(" :: Circular = (%e, %e)\n", rz[0], rz[1]);

	/********/
	/* ITER */
	/********/
	magnetic_handler *numeric = magnetic_handler_select("numeric");
	struct general_settings *mset_numeric = malloc(sizeof(struct general_settings));
	mset_numeric->name = "numeric";
	mset_numeric->n = 1;
	mset_numeric->setting = malloc(mset_numeric->n*sizeof(char*));
	mset_numeric->value = malloc(mset_numeric->n*sizeof(char*));
	mset_numeric->setting[0] = setname("file");
	mset_numeric->value[0]   = setname("../resources/ITER/iter2d.bkg");

	numeric->init(mset_numeric);
	numeric->init_run();
	numeric->init_particle();

	rz = magnetic_axis_find(6.);
	printf(" :: ITER     = (%e, %e)\n", rz[0], rz[1]);

	/**********/
	/* DIII-D */
	/**********/
	numeric = magnetic_handler_select("numeric");
	mset_numeric = malloc(sizeof(struct general_settings));
	mset_numeric->name = "numeric";
	mset_numeric->n = 1;
	mset_numeric->setting = malloc(mset_numeric->n*sizeof(char*));
	mset_numeric->value = malloc(mset_numeric->n*sizeof(char*));
	mset_numeric->setting[0] = setname("file");
	mset_numeric->value[0]   = setname("../resources/DIII-D/DIII-D.bkg");

	numeric->init(mset_numeric);
	numeric->init_run();
	numeric->init_particle();

	rz = magnetic_axis_find(1.3);
	printf(" :: DIII-D   = (%e, %e)\n", rz[0], rz[1]);
}
