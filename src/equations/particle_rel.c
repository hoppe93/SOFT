#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "equation_particle_rel.h"
#include "global.h"
#include "interp2.h"
#include "IO_data.h"
#include "magnetic_field.h"
#include "quantities.h"
#include "readfile.h"
#include "rkf45.h"
#include "tools.h"
#include "vector.h"

particle *particle_rel_particle=NULL;
/* Pre-allocated vectors */
vector *particle_rel_result;

#pragma omp threadprivate(particle_rel_particle,particle_rel_result)

void equation_particle_rel_init(settings *set) {
	/* Setup result vector */
	particle_rel_result = vnew(6);
}

void equation_particle_rel_init_run(particle *part, ode_solution *solver_object) {
	particle_rel_particle = part;

	double c=299792458;		// speed of light
	double vx = part->v0[0],
		   vy = part->v0[1],
		   vz = part->v0[2];
	double v2 = vx*vx + vy*vy + vz*vz;
	double gamma = 1/sqrt(1-v2/(c*c));

	if (part->gc_position) {
		vector *B = magnetic_field_get(part->r0[0], part->r0[1], part->r0[2]);
		double rho[3], bx, by, bz, Babs,
			   a[3], c[3], aabs, sn, cs,
			   rhoabs, pperp[3], pabs2,
			   pperpabs, pparabs;

		Babs = vnorm3(B);
		bx = B->val[0] / Babs;
		by = B->val[1] / Babs;
		bz = B->val[2] / Babs;

		aabs = hypot(bx, by);
		a[0] = by / aabs;
		a[1] =-bx / aabs;
		a[2] = 0;

		c[0] = a[1]*bz;	/* a[2] = 0 */
		c[1] =-a[0]*bz;	/* - " - */
		c[2] = a[0]*by - a[1]*bx;

		sn = sin(part->zeta0), cs = cos(part->zeta0);
		pabs2 = part->v0[0]*part->v0[0] + part->v0[1]*part->v0[1] + part->v0[2]*part->v0[2];
		pparabs = part->v0[0]*bx + part->v0[1]*by + part->v0[2]*bz;
		pperpabs = sqrt(pabs2 - pparabs*pparabs);

		rhoabs = fabs(pperpabs / (part->charge * Babs));
		rho[0] = rhoabs * (a[0]*cs + c[0]*sn);
		rho[1] = rhoabs * (a[1]*cs + c[1]*sn);
		rho[2] = rhoabs * (a[2]*cs + c[2]*sn);

		pperp[0] = pperpabs * (c[0]*cs - a[0]*sn);
		pperp[1] = pperpabs * (c[1]*cs - a[1]*sn);
		pperp[2] = pperpabs * (c[2]*cs - a[2]*sn);

		solver_object->Z->val[0] = part->r0[0] + rho[0];
		solver_object->Z->val[1] = part->r0[1] + rho[1];
		solver_object->Z->val[2] = part->r0[2] + rho[2];
		solver_object->Z->val[3] = pparabs*bx + pperp[0];
		solver_object->Z->val[4] = pparabs*by + pperp[1];
		solver_object->Z->val[5] = pparabs*bz + pperp[2];
	} else {
		solver_object->Z->val[0] = part->r0[0];
		solver_object->Z->val[1] = part->r0[1];
		solver_object->Z->val[2] = part->r0[2];
		solver_object->Z->val[3] = part->v0[0]*part->mass*gamma;
		solver_object->Z->val[4] = part->v0[1]*part->mass*gamma;
		solver_object->Z->val[5] = part->v0[2]*part->mass*gamma;
	}
}

/**
 * Initialize output data structure
 */
solution_data *equation_particle_rel_init_output(solution_data *output) {
	output->nvars = 6;
	output->labels = malloc(sizeof(char*)*output->nvars);
	output->labels[0] = "x";
	output->labels[1] = "y";
	output->labels[2] = "z";
	output->labels[3] = "px";
	output->labels[4] = "py";
	output->labels[5] = "pz";

	return output;
}

/**
 * Equation for the charged particle motion
 * Lorentz force, only magnetic field.
 * 
 * T: time, not used here but needed for ode.c
 * Z: Pointer to vector containing particle position
 * coordinates in the first three values, and particle
 * velocity in the next three.
 *
 * RETURNS: vector of values of function f 
 *
 * Used as the first argument to ode_solve in ode.c if no GCM is selected
 */
vector * equation_particle_rel_eq(double T, vector* Z) {
	particle *part = particle_rel_particle;
	double m=part->mass; // particle mass
	double e=part->charge; // particle charge
	double c=299792458;		// speed of light

	double x=Z->val[0],
		   y=Z->val[1],
		   z=Z->val[2];
	/* Save  particle velocity*/
	double p1=Z->val[3],
		   p2=Z->val[4],
		   p3=Z->val[5];

	/* Calculate particle scalar velocity squared */
	double psqr=p1*p1 + p2*p2 + p3*p3;

	/* Get magnetic field in point of particle */
	vector *B = magnetic_field_get(x,y,z);

	/* Extract x, y and z values of magnetic field */
	double B1=B->val[0],
	B2=B->val[1],
	B3=B->val[2];

	/* Calculate gamma-dot factor */
	double gamma= sqrt(1+psqr/(c*c*m*m));

	/* Calculate each function (f) value */
	double f1=p1/(gamma*m),
	f2=p2/(gamma*m),
	f3=p3/(gamma*m),
	f4=(e/(gamma*m))*(p2*B3-p3*B2),
	f5=(e/(gamma*m))*(p3*B1-p1*B3),
	f6=(e/(gamma*m))*(p1*B2-p2*B1);

	double Bmag = sqrt(B1*B1 + B2*B2 + B3*B3);
	double bhatx=B1/Bmag, bhaty=B2/Bmag, bhatz=B3/Bmag;

	double ppar = p1*bhatx + p2*bhaty + p3*bhatz,
		   pperp = sqrt(psqr - ppar*ppar),
		   pmag = sqrt(psqr);

	tool_record_data(
		T, ppar*ppar, pperp*pperp, Bmag,
		x, y, z,
		p1/pmag, p2/pmag, p3/pmag,
		ppar/(gamma*m), pperp/(gamma*m)
	);

	/* Save in vector and return */
	vector *result = particle_rel_result;
	result->val[0]=f1;
	result->val[1]=f2;
	result->val[2]=f3;
	result->val[3]=f4;
	result->val[4]=f5;
	result->val[5]=f6;

	return result;
}

