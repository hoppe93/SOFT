#include <math.h>
#include <omp.h>
#include <stdlib.h>

#include "equation_particle.h"
#include "global.h"
#include "interp2.h"
#include "IO_data.h"
#include "magnetic_field.h"
#include "quantities.h"
#include "readfile.h"
#include "rkf45.h"
#include "tools.h"
#include "vector.h"

particle *particle_particle=NULL;
/* Pre-allocated vectors */
vector *particle_result=NULL;

#pragma omp threadprivate(particle_particle,particle_result)

void equation_particle_init(settings *set) {
	/* Pre-allocate vectors */
	particle_result = vnew(6);
}

void equation_particle_init_run(particle *part, ode_solution *solver_object) {
	particle_particle = part;

	/* Setup solution vector */
	solver_object->Z->val[0] = part->r0[0];
	solver_object->Z->val[1] = part->r0[1];
	solver_object->Z->val[2] = part->r0[2];
	solver_object->Z->val[3] = part->v0[0];
	solver_object->Z->val[4] = part->v0[1];
	solver_object->Z->val[5] = part->v0[2];
}

/**
 * Initialize output data structure
 */
solution_data *equation_particle_init_output(solution_data *output) {
	output->nvars = 6;
	output->labels = malloc(sizeof(char*)*output->nvars);
	output->labels[0] = "x";
	output->labels[1] = "y";
	output->labels[2] = "z";
	output->labels[3] = "vx";
	output->labels[4] = "vy";
	output->labels[5] = "vz";

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
vector * equation_particle_eq(double T, vector* Z){
	particle *part = particle_particle;
	double m=part->mass; // particle mass
	double e=part->charge; // particle charge

	double x=Z->val[0],
		   y=Z->val[1],
		   z=Z->val[2];
	/* Save  particle velocity*/
	double v1=Z->val[3],
		   v2=Z->val[4],
		   v3=Z->val[5];

	/* Get magnetic field in point of particle */
	vector *B = magnetic_field_get(x,y,z);

	/* Extract x, y and z values of magnetic field */
	double B1=B->val[0],
		   B2=B->val[1],
		   B3=B->val[2];
	/* Calculate each function (f) value */
	double f1=v1,
		   f2=v2,
		   f3=v3,
		   f4=(e/m)*(v2*B3-v3*B2),
		   f5=(e/m)*(v3*B1-v1*B3),
		   f6=(e/m)*(v1*B2-v2*B1);

	double Bmag = hypot(B1, hypot(B2, B3));
	double vpar = v1*B1/Bmag + v2*B2/Bmag + v3*B3/Bmag;
	double v = hypot(v1, hypot(v2, v3));
	double vperp = sqrt(v*v - vpar*vpar);
	double ppar2 = m*m*vpar*vpar;
	double pperp2 = m*m*vperp*vperp;
	tool_record_data(
		T,
		ppar2, pperp2,
		Bmag,
		x, y, z,
		v1/v, v2/v, v3/v,
		vpar, vperp
	);

	/* Save in vector and return */
	vector *result = particle_result;
	result->val[0]=f1;
	result->val[1]=f2;
	result->val[2]=f3;
	result->val[3]=f4;
	result->val[4]=f5;
	result->val[5]=f6;

	return result;
}

