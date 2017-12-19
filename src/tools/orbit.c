/* Orbit following tool */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arguments.h"
#include "constants.h"
#include "ctsv.h"
#include "IO_data.h"
#include "global.h"
#include "magnetic_field.h"
#include "orbit.h"
#include "rkf45.h"
#include "tools.h"

double max_time=0, tstart=0, orbit_mass=0, *simtime, **solution_values=NULL;
unsigned int time_index=0, allocated=0, points=0, orbit_variables=0,
	solution_values_length=0, *solution_values_sizes=NULL,
	QUANTITY_KINEN, QUANTITY_XI, QUANTITY_B, QUANTITY_BX, QUANTITY_BY,
	QUANTITY_BZ, QUANTITY_GAMMA, QUANTITY_MU;
vector *solution=NULL;
char *orbit_output_file=NULL;
#define DEFAULT_POINTS 1000

void orbit_init(struct general_settings *set, struct general_settings *outset, int nouts) {
	/* Parse settings */
	int i;
	for (i = 0; i < set->n; i++) {
		if (!strcmp(set->setting[i], "output"))
			orbit_output_file = set->value[i];
		else {
			fprintf(stderr, "ERROR: Unknown setting for tool orbit: %s!\n", set->setting[i]);
			exit(-1);
		}
	}
}
void orbit_init_run(unsigned int variables) {
	orbit_variables = variables;
}
ode_solution *orbit_init_particle(particle *p) {
	/* Allocate various structures sitting in memory */
	orbit_allocate();

	tstart = p->t0;
	simtime[0] = p->t0;
	orbit_mass = p->mass;
	max_time = p->tend;

	quantities_init();

	/* Define quantity */
	QUANTITY_KINEN=quantities_define("Energy");
	QUANTITY_XI=quantities_define("xi");
	QUANTITY_B=quantities_define("B");
	QUANTITY_BX=quantities_define("BX");
	QUANTITY_BY=quantities_define("BY");
	QUANTITY_BZ=quantities_define("BZ");
	QUANTITY_GAMMA=quantities_define("gamma");
	QUANTITY_MU=quantities_define("mu");

	/* Report quantities */
	//double kinenergy = orbit_mass * sqrt(p->v0[0]*p->v0[0] + p->v0[1]*p->v0[1] + p->v0[2]*p->v0[2]) / 2;
	double gamma = 1/sqrt(1-(p->v0[0]*p->v0[0] + p->v0[1]*p->v0[1] + p->v0[2]*p->v0[2])/(LIGHTSPEED*LIGHTSPEED));
	double kinenergy = orbit_mass * LIGHTSPEED*LIGHTSPEED * (gamma - 1);
	quantities_report(QUANTITY_KINEN, kinenergy, 0);

	/* Find xi */
	double vtot = sqrt(p->vpar*p->vpar + p->vperp*p->vperp);
	quantities_report(QUANTITY_XI, p->vpar/vtot, 0);

	/* Save magnetic field strength */
	vector *mf = magnetic_field_get(p->r0[0], p->r0[1], p->r0[2]);
	double B = sqrt(mf->val[0]*mf->val[0] + mf->val[1]*mf->val[1] + mf->val[2]*mf->val[2]);
	quantities_report(QUANTITY_B, B, 0);

	/* Save gamma factor */
	quantities_report(QUANTITY_GAMMA, gamma, 0);

	/* Save mu */
	double pperp = gamma*orbit_mass * p->vperp;
	quantities_report(QUANTITY_MU, pperp*pperp / (2*orbit_mass*B), 0);

	/* Allocate solution vector list */
	solution[0].val = solution_values[solution_values_length-1];
	solution[0].n   = orbit_variables;
	solution[1].val = solution_values[solution_values_length-1]+orbit_variables;
	solution[1].n   = orbit_variables;

	/* Create the ode_solution object and point it to the solution */
	ode_solution *solver_object;
	solver_object = malloc(sizeof(ode_solution));
	solver_object->Z = solution;
	solver_object->result = solution+1;
	solver_object->step = 1e-8;
	points++;

	time_index++;

	return solver_object;
}
void orbit_deinit_run(void) {}
/* There is no stop condition for the orbit tool. Always return FALSE */
int orbit_stop_condition(void) { return 0; }
/* The orbit tool has no 'in-step' function. This is just a dummy */
void orbit_step(ode_solution *solver_object, step_data *sd) {
	/* Check if more memory needs to be allocated */
	if (time_index+1 >= allocated)
		orbit_allocate();

	/* Initialize new result vector */
	int minus = 0;
	if (solution_values_length > 1)
		minus = solution_values_sizes[solution_values_length-2];

	int next = orbit_variables*(time_index+1-minus);
	solution[time_index+1].val = solution_values[solution_values_length-1]+next;
	solution[time_index+1].n   = orbit_variables;

	/* Advance the step */
	solver_object->Z = solution+time_index;		/* Result becomes new input */
	solver_object->result = solution+time_index+1;
	points++;

	/* Save energy */
	//double kinen = orbit_mass*(sd->vpar*sd->vpar + sd->vperp*sd->vperp)/2;
	double vtot = (sd->vpar*sd->vpar + sd->vperp*sd->vperp);
	double gamma = 1/sqrt(1-vtot/(LIGHTSPEED*LIGHTSPEED));
	double kinen = orbit_mass * LIGHTSPEED*LIGHTSPEED * (gamma - 1);
	double pperp = gamma * orbit_mass * sd->vperp;
	double mu = pperp*pperp / (2*orbit_mass*sd->B);
	vector *Bvec = magnetic_field_get(sd->x, sd->y, sd->z);

	quantities_report(QUANTITY_KINEN, kinen, time_index);
	quantities_report(QUANTITY_XI, sd->vpar/sqrt(vtot), time_index);
	quantities_report(QUANTITY_B, sd->B, time_index);
	quantities_report(QUANTITY_BX, Bvec->val[0], time_index);
	quantities_report(QUANTITY_BY, Bvec->val[1], time_index);
	quantities_report(QUANTITY_BZ, Bvec->val[2], time_index);
	quantities_report(QUANTITY_GAMMA, gamma, time_index);
	quantities_report(QUANTITY_MU, mu, time_index);

	/* Store time */
	simtime[time_index++] = sd->time;
}
void orbit_output(equation *eq) {
	solution_data *output;
	output = malloc(sizeof(solution_data));
	output->T = simtime;
	output->v = solution;
	output->points = points;
	output->quantities = quantities_get();
	output->no_quantities = quantities_get_no();

	output = eq->output(output);
	ctsv_write(orbit_output_file, ',', output, NULL);
}

void orbit_allocate(void) {
	unsigned int np = DEFAULT_POINTS;
	if (time_index!=0) {
		if (max_time > tstart)
			np = (unsigned int)ceil(allocated * max_time/(simtime[time_index-1]-tstart));
		else	/* We have no way to predict how many points are needed. Increase linearly. */
			np = allocated + DEFAULT_POINTS;
	}

	solution = realloc(solution, np*(sizeof(vector)));
	solution_values = realloc(solution_values, sizeof(double*)*(solution_values_length+1));
	solution_values_sizes = realloc(solution_values_sizes, sizeof(unsigned int)*(solution_values_length+1));

	solution_values[solution_values_length] = malloc((np-allocated)*(sizeof(double)*orbit_variables));
	solution_values_sizes[solution_values_length] = np;
	solution_values_length++;

	simtime = realloc(simtime, np*sizeof(double));
	quantities_expand(np);

	allocated = np;
}
