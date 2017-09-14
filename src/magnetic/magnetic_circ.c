/* Circular magnetic field handler */
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "domain.h"
#include "magnetic_axis.h"
#include "magnetic_circ.h"

struct magnetic_circ_params *mcpar=NULL;
vector *magnetic_circ_retval, *magnetic_circ_gradB, *magnetic_circ_curlB;
double magnetic_circ_oldx=0., magnetic_circ_oldy=0., magnetic_circ_oldz=0.;
diff_data *magnetic_circ_diff_dd=NULL;

#pragma omp threadprivate(magnetic_circ_retval,magnetic_circ_gradB,magnetic_circ_curlB,magnetic_circ_oldx,magnetic_circ_oldy,magnetic_circ_oldz,magnetic_circ_diff_dd)

void magnetic_circ_init(struct general_settings *set) {
	if (set == NULL) {
		fprintf(stderr, "ERROR: The circular magnetic field handler requires settings to be given!\n");
		exit(-1);
	}

	/* Loop over all settings */
	int i;
	mcpar = malloc(sizeof(struct magnetic_circ_params));
	mcpar->B0 = 1.0;
	mcpar->Rm = 1.0;
	mcpar->r  = 0.5;
	mcpar->q  = 1.0;

	for (i = 0; i < set->n; i++) {
		if (!strcmp(set->setting[i], "B0")) {
			mcpar->B0 = atof(set->value[i]);
		} else if (!strcmp(set->setting[i], "major_radius")) {
			mcpar->Rm = atof(set->value[i]);
			if (mcpar->Rm <= 0.0) {
				fprintf(stderr, "ERROR: Major radius must be greater than zero.\n");
				exit(EXIT_FAILURE);
			}
		} else if (!strcmp(set->setting[i], "minor_radius")) {
			mcpar->r = atof(set->value[i]);
			if (mcpar->r <= 0.0) {
				fprintf(stderr, "ERROR: Minor radius must be greater than zero.\n");
				exit(EXIT_FAILURE);
			}
		} else if (!strcmp(set->setting[i], "safety_factor")) {
			mcpar->q = atof(set->value[i]);
			if (mcpar->q <= 0.0) {
				fprintf(stderr, "ERROR: Safety factor must be greater than zero.\n");
				exit(EXIT_FAILURE);
			}
		} else {
			fprintf(stderr, "Invalid magnetic field handler 'circular' setting: %s!", set->setting[i]);
			exit(-1);
		}
	}

	magnetic_axis_set_loc(mcpar->Rm, 0.);

	/* Generate domain */
	int n = 100;
	double *r = malloc(sizeof(double)*n),
		   *z = malloc(sizeof(double)*n);
	for (i = 0; i < n; i++) {
		r[i] = mcpar->r * cos(2.0*PI * i / (double)n);
		z[i] = mcpar->r * sin(2.0*PI * i / (double)n);
	}

	domain_set(r, z, n);
}

void magnetic_circ_init_run(void) {
	magnetic_circ_retval = vnew(3);
	magnetic_circ_gradB = vnew(3);
	magnetic_circ_curlB = vnew(3);
	magnetic_circ_diff_dd = malloc(sizeof(diff_data));
}
void magnetic_circ_init_particle(void) { }
vector *magnetic_circ_eval(double x, double y, double z) {
	/* If we computed the B field in this point
	 * the last time, magnetic_circ_retval
	 * already contains the appropriate value */
	if (x == magnetic_circ_oldx &&
		y == magnetic_circ_oldy &&
		z == magnetic_circ_oldz)
		return magnetic_circ_retval;

	double x2y2 = mcpar->Rm - hypot(x, y);
	double r = hypot(x2y2, z);
	/*double theta = atan2(z, x2y2);
	double phi = atan2(y, x);*/

	double sintheta = z/r, costheta = x2y2/r;
	double sinphi = y/(mcpar->Rm-r*costheta), cosphi = x/(mcpar->Rm-r*costheta);
	/*double sinphi = sin(phi), cosphi = cos(phi);
	double sintheta = sin(theta), costheta = cos(theta);*/

	double pf = mcpar->B0 / (1 - r*costheta/mcpar->Rm);
	double t1 = r/(mcpar->q*mcpar->Rm);
	magnetic_circ_retval->val[0] = pf * (t1 * (sintheta*cosphi) + sinphi);
	magnetic_circ_retval->val[1] = pf * (t1 * (sintheta*sinphi) - cosphi);
	magnetic_circ_retval->val[2] = pf * (t1 * costheta);

	return magnetic_circ_retval;
}
diff_data *magnetic_circ_diff(double x, double y, double z) {
	double x2y2 = mcpar->Rm - hypot(x, y);

	/* Switch to toroidal coordinates */
	double r = hypot(z, x2y2);

	double sintheta = z/r, costheta = x2y2/r;
	double sinphi = y/(mcpar->Rm-r*costheta), cosphi = x/(mcpar->Rm-r*costheta);

	double sub = mcpar->Rm/(mcpar->Rm-r*costheta);
	double sub2 = sub*sub;
	double r_qRm = r/(mcpar->q*mcpar->Rm), r_qRm2 = r_qRm*r_qRm;
	double sqr = sqrt(r_qRm2 + 1);

	/* Derivatives of field strength */
	double dB_dr = mcpar->B0*sub*(costheta*sub/mcpar->Rm*sqr + r_qRm2/(r*sqr));
	double dB_dt =-mcpar->B0*sintheta*sub2*sqr/mcpar->Rm;

	magnetic_circ_gradB->val[0] = -dB_dr*costheta*cosphi + dB_dt*sintheta*cosphi;
	magnetic_circ_gradB->val[1] = -dB_dr*costheta*sinphi + dB_dt*sintheta*sinphi;
	magnetic_circ_gradB->val[2] = dB_dr*sintheta + dB_dt*costheta;

	/* Derivatives of field components */
	double curlB_p = mcpar->B0/(mcpar->q*mcpar->Rm*mcpar->Rm)*(2*mcpar->Rm-r*costheta)*sub2;

	magnetic_circ_curlB->val[0] =-curlB_p*sinphi;
	magnetic_circ_curlB->val[1] = curlB_p*cosphi;
	//magnetic_circ_curlB->val[2] = 2*mcpar->B0*sub2/(mcpar->Rm-r*costheta);
	magnetic_circ_curlB->val[2] = 0.;

	magnetic_circ_diff_dd->B = magnetic_circ_eval(x,y,z);
	magnetic_circ_diff_dd->Babs = mcpar->B0*sub*sqr;
	magnetic_circ_diff_dd->gradB = magnetic_circ_gradB;
	magnetic_circ_diff_dd->curlB = magnetic_circ_curlB;

	return magnetic_circ_diff_dd;
}
diff_data *magnetic_circ_diff_notor(double x, double y, double z) {
	double x2y2 = mcpar->Rm - hypot(x, y);

	/* Switch to toroidal coordinates */
	double r2 = z*z + x2y2*x2y2;
	double r = sqrt(r2);

	double sintheta = z/r, costheta = x2y2/r;
	double sinphi = y/(mcpar->Rm-r*costheta), cosphi = x/(mcpar->Rm-r*costheta);

	double sub = 1/(mcpar->Rm-r*costheta);
	double sub2 = sub*sub;
	double r_qRm = r/(mcpar->q*mcpar->Rm), r_qRm2 = r_qRm*r_qRm;
	double sqr = sqrt(r_qRm2 + 1);

	/* Derivatives of field strength */
	double dB_dr = mcpar->B0*mcpar->Rm/mcpar->q*sub2;
	double dB_dt = mcpar->B0/mcpar->q*sub2*r*sintheta;

	magnetic_circ_gradB->val[0] =-dB_dr*costheta*cosphi + dB_dt*sintheta*cosphi;
	magnetic_circ_gradB->val[1] =-dB_dr*costheta*sinphi + dB_dt*sintheta*sinphi;
	magnetic_circ_gradB->val[2] = dB_dr*sintheta + dB_dt*costheta;

	/* Derivatives of field components */
	double curlB_p = mcpar->B0/(mcpar->q*mcpar->Rm*mcpar->Rm)*(2*mcpar->Rm-r*costheta)*sub2;

	magnetic_circ_curlB->val[0] =-curlB_p*sinphi;
	magnetic_circ_curlB->val[1] = curlB_p*cosphi;
	magnetic_circ_curlB->val[2] = 0.;

	magnetic_circ_diff_dd->B = magnetic_circ_eval(x,y,z);
	magnetic_circ_diff_dd->Babs = mcpar->B0*sub*sqr;
	magnetic_circ_diff_dd->gradB = magnetic_circ_gradB;
	magnetic_circ_diff_dd->curlB = magnetic_circ_curlB;

	return magnetic_circ_diff_dd;
}
