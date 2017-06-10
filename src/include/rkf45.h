#ifndef _ODE_H
#define _ODE_H

#define REDO_STEP 1
#define OK_STEP 0

#include "vector.h"

typedef struct {
	vector *Z;     /* Solution points, z = z(t) */
	vector *result;/* Result of solver step */
	double step;   /* Optimal step */
	double actualstep;/* Actual stepsize taken */
	int flag;	   /* indicates if REDO_STEP (1) or OK_STEP (0) */
} ode_solution;

void ode_init(unsigned int);
/**
 * Solve an Initial Value Problem (IVP ODE)
 *
 * RETURNS a solution to the equation as a pointer to an ode_solution consisting
 * of Z, step, flag.
 * Now Z contains the new solution points, new step size and new flag value.
 *
 */
ode_solution* ode_solve(vector *(equation)(double, vector*),ode_solution*, double);
/**
 * Calculates new solution points
 *
 */
vector * ode_step(vector *(*equation)(double, vector*),ode_solution*, double ,unsigned int);

void ode_test(void);

extern double ode_tolerance, ode_maxtimestep;

#endif/*_ODE_H*/
