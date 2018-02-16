/* ODE Solver */
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "ctsv.h"
#include "equations.h"
#include "rkf45.h"
#include "vector.h"

/* 1e-7 seems pretty much the best.
 */
#define EPS0 1e-12                /* error tolerance */
#define SAFETY_FACTOR 0.9	/* Safety factor beta */
#define NUMBER_OF_TESTPOINTS 5000 /* for ode_test */
#define ORDER1 4
#define ORDER2 5

/* used if explicit time dependence */
//double c[]={1.0/5,3.0/10,3.0/5.0,1,7.0/8};
double c[]={0,1.0/4.0,3.0/8.0,12.0/13.0,1,1.0/2.0};

/* Stores Fehlberg coefficients A */
double A[6][6] = {
	{0,0,0,0,0,0},
	{1.0/4.0,0,0,0,0,0},
	{3.0/32.0,9/32.0,0,0,0,0},
	{1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0, 0, 0},
	{439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0, 0},
	{-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0}
};
/* Stores Fehlberg coefficients. Contains
   b in first row and bhat in second */
double B[2][6] = {
	{25.0/216.0, 0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0},
	{16.0/135.0, 0, 6656.0/12825.0,28561.0/56430.0,-9.0/50.0,2.0/55.0}
};

/* Pre-allocated vectors */
vector *ode_K=NULL, *ode_sum=NULL, *ode_sum1=NULL, *ode_sum2=NULL;
double ode_tolerance=EPS0, ode_maxtimestep=0;

#pragma omp threadprivate(ode_K,ode_sum,ode_sum1,ode_sum2)

/**
 * Initialize the ODE solver
 */
void ode_init(unsigned int variables) {
	int i;
	ode_K = malloc(sizeof(vector)*(ORDER2+1));
	for (i = 0; i <= ORDER2; i++) {
		ode_K[i].val = malloc(sizeof(double)*variables);
		ode_K[i].n = variables;
	}

	ode_sum = vnew(variables);
	ode_sum1= vnew(variables);
	ode_sum2= vnew(variables);
}

/**
 * Solve an Initial Value Problem (IVP ODE)
 *
 * equation: Pointer to function giving the equation,
 * that is f in "z' = f(t,z)" where t is time
 * (scalar) and z is the unknown (vector)
 *
 * solver_object: pointer to a defined type ode_solution consisting 
 * of Z, step, flag.
 * Z contains the solution values from the previous iteration.
 * step is the optimal step size. flag indicates REDO_STEP (1) or OK_STEP (0).
 *
 * T: time
 *
 * RETURNS a solution to the equation as a pointer to an ode_solution consisting
 * of Z, step, flag.
 * Now Z contains the new solution points, new step size and new flag value.
 *
 * Called from main
 */
ode_solution* ode_solve(vector *(*equation)(double, vector*), ode_solution *solver_object, double T) {
	/* Choose order of iteration */
	unsigned int order1=ORDER1;
	unsigned int order2=ORDER2;
	double eps0=ode_tolerance, eps, epst; /* Tolerance parameter */
	double beta=SAFETY_FACTOR; /* Safety parameter */
	/* Stores flag indicating whether the the iteration should be redone */
	int flag = solver_object->flag;
	/* current solution points */
	vector* Z = solver_object->Z;
	/* current step size */
	double h=solver_object->step;
	solver_object->actualstep = h;

	/* To store optimal steplenght*/
	double hopt;
	/* Loop variable */
	unsigned int i, j;

	/* Calculate next point */
	vector* K = ode_step(equation, solver_object, T, order2);

	/* Calculate sum to be used in next point for Z_next and Zhat */
	/* Initialize sum */
	vector *sum1=ode_sum1, *sum2=ode_sum2;
	for (i = 0; i < sum1->n; i++) {
		sum1->val[i] = 0;
		sum2->val[i] = 0;
	}

	/* Sum k1 to k5 (4th order method has "5 stages" (e.g. 5 k's) */
	for (i=0; i<=order1; i++){
		for (j = 0; j < sum1->n; j++)
			sum1->val[j] += h*B[0][i]*K[i].val[j];
	}

	/* Calculate next point (add Z to sum1 vector) */
	vector *Z_next = vaddf(sum1, Z);

	/* Sum k1 to k6 (5th order method has "6 stages" (e.g. 6 k's) */
	for (i=0; i <= order2; i++){
		for (j = 0; j < sum1->n; j++)
			sum2->val[j] += h*B[1][i]*K[i].val[j];
	}

	/* Add sum2 and Z (into sum2 vector) */
	vector *Zhat = vaddf(sum2, Z);

	/* calculate relative error epsilon for each component,
	 * save the largest error value */
	/* Calculate the maximum error made */
    double Zhateps = Zhat->val[0];
	eps = fabs((Zhat->val[0] - Z_next->val[0]) /  Z_next->val[0]);
	for (i = 1; i < Z_next->n; i++) {
		epst = fabs((Zhat->val[i] - Z_next->val[i]) / Z_next->val[i]);
		if (epst > eps) {
            if (eps / epst < 1e-6 && fabs(Zhateps / Zhat->val[i]) > 1e6) continue;
            else {
                eps = epst;
                Zhateps = Zhat->val[i];
            }
		}
	}
	hopt = h; /* Optimal step size */

	/* Choose optimal step */
	if (eps >= eps0) {
		hopt=beta*h*pow(eps0/eps,0.20);
		flag=REDO_STEP;
	} else if (eps == 0) {
		hopt = h;
		flag=OK_STEP;
	} else {
		hopt=h*pow(eps0/eps,0.25);
		flag=OK_STEP;
	}

	if (ode_maxtimestep > 0 && hopt > ode_maxtimestep)
		hopt = ode_maxtimestep;
	else if (hopt <= 0.0) {
		hopt = ode_maxtimestep;
		flag = REDO_STEP;
	}

	/* Save and return calculated values */
	solver_object->step = hopt;
	solver_object->flag = flag;
	vequate(solver_object->result, Z_next);

	/* Return the solver object */
	return solver_object;
}
/**
 * Calculates new solution points
 *
 * equation: Pointer to function giving the equation,
 * that is f in "z' = f(t,z)" where t is time
 * (scalar) and z is the unknown (vector)
 *
 * solver_object: pointer to a defined type ode_solution consisting 
 * of Z, step, flag.
 * Z contains the solution values from the previous iteration.
 * step is the optimal step size. flag indicates REDO_STEP (1) or OK_STEP (0).
 *
 * T: time
 * order: order of RK-algorithm
 *
 * RETURNS array with K values
 *
 * Called from ode_solve.
 */
vector *ode_step(vector *(*equation)(double, vector*), ode_solution *solver_object, double T, unsigned int order){
	vector* Z = solver_object->Z;
	double h = solver_object->step;
	/* loop variables*/
	unsigned int i, j, k;
	vector *K = ode_K, *sum = ode_sum;

	/* Cacluate first K */
	vector *eq = (*equation)(T,Z);
	vequate(ode_K, eq);

	/* Calculate each K up to order. Start from K2 (i=1).
	* Note that a 4th order method contains 5 k's and
	* a 5th order method contains 6 k's */
	for (i=1; i <= order; i++) {
		/* Initialize sum */
		for (j = 0; j < Z->n; j++)
			sum->val[j] = 0;

		/* Compute sum to use in argument */
		for (j=0; j < i; j++) {
			for (k = 0; k < K->n; k++) {
				sum->val[k] += h*A[i][j] * K[j].val[k];
			}
		}

		/* Calculate K */
		vaddf(sum, Z);
		vequate(K+i, (*equation)(T+c[i]*h, sum));
	}

	return K;
}

/**
   Test function for this module 
*/
/*
void ode_test(void) {
	/ * Iteration variable * /
	unsigned int i=0;

	/ * Predator prey model * /
	/ * Initiate vector to store calculated points * /
	vector* coordinates;
	unsigned int points = NUMBER_OF_TESTPOINTS;
	coordinates=malloc(sizeof(vector)*(points+1));

	/ * Set initial point * /
	coordinates->val = malloc(sizeof(double)*2);
	coordinates->n = 2;

	/ * Initial condition:nbr of predators and prey * /
	coordinates->val[0] = 34.91;
	coordinates->val[1] = 3.857;

	/ * To store time * /
	double *t = malloc(sizeof(double)*(points+1));
	t[0] = 0;
	/ * To store dummy ''Energy'' * /
	double *E = malloc(sizeof(double)*(points+1));
	for (i = 0; i < points; i++) {
		E[i] = 0;
	}

	/ * Choose starting steplenght * /
	double h=0.3;

	/ * Save everything in type 'ode_solution' * /
	ode_solution *param;
	param = malloc(sizeof(ode_solution));
	param->step=h;

	/ * First value of flag is 'OK'* /
	param->flag=0;

	/ * Iterate and save!* /
	/ * Choose size of intervall * /

	for (i = 0; i < points; i++) {
		t[i+1]=t[i]+param->step;

		/ * Iterate once * /
		param->Z = coordinates+i;
		//    ode_solve(equation_predator_prey, param, t[i]);

		/ * Check if iteration needs to be re-done * /
		/ * flag=0 -> ok -> continue * /
		if (param->flag!=0) i=i-1; // Redo step with new calculated h in param
	}

	solution_data output;
	output.T=t;
	output.v=coordinates;
	output.labels=malloc(sizeof(char *)*2);
	output.labels[0]="x";
	output.labels[1]="y";
	output.points=points;
	output.nvars=2;

	ctsv_write("Output.csv",',',&output, NULL);
}*/

