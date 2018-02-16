/* Tool manager */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "orbit.h"
#include "sycamera.h"
#include "tools.h"
#include "util.h"

const int NUMBER_OF_TOOLS=2;
tool *all_tools=NULL, *selected_tool=NULL;

step_data *TOOL_STEPDATA=NULL;
#pragma omp threadprivate(TOOL_STEPDATA)

void tool_init_handler(void) {
	all_tools = malloc(sizeof(tool)*NUMBER_OF_TOOLS);
	
	all_tools[0].name = setname("orbit");
	all_tools[0].tolerance = 1e-7;
	all_tools[0].require_jacobian = 0;
	all_tools[0].init = orbit_init;
	all_tools[0].init_run = orbit_init_run;
	all_tools[0].init_particle = orbit_init_particle;
	all_tools[0].deinit_run = orbit_deinit_run;
	all_tools[0].step = orbit_step;
	all_tools[0].stop_condition = orbit_stop_condition;
	all_tools[0].output = orbit_output;

	all_tools[1].name = setname("sycamera");
	all_tools[1].tolerance = 1e-12;
	all_tools[1].require_jacobian = 1;
	all_tools[1].init = sycamera_init;
	all_tools[1].init_run = sycamera_init_run;
	all_tools[1].init_particle = sycamera_init_particle;
	all_tools[1].deinit_run = sycamera_deinit_run;
	all_tools[1].step = sycamera_step;
	all_tools[1].stop_condition = sycamera_stop_condition;
	all_tools[1].output = sycamera_output;

	//TOOL_STEPDATA = malloc(sizeof(step_data));
}

void tool_prepare_run(void) {
	TOOL_STEPDATA = malloc(sizeof(step_data));
}

tool *tool_select(char *name) {
	int i;
	for (i = 0; i < NUMBER_OF_TOOLS; i++) {
		if (!strcmp(all_tools[i].name, name)) {
			printf("Selecting tool '%s'...\n", all_tools[i].name);
			selected_tool = all_tools+i;
			ode_tolerance = all_tools[i].tolerance;
			return selected_tool;
		}
	}

	return NULL;
}

void tool_record_data(
	/* Current time */
	double T,
	/* Parallel and perpendicular momentum squared */
	double ppar2, double pperp2,
	/* Magnetic field strength */
	double B,
	/* Particle/guiding-center position */
	double x, double y, double z,
	/* Particle/guiding-center velocity */
	double vhatx, double vhaty, double vhatz,
	/* Particle/guiding-center parallel speed */
	double vpar, double vperp
) {
	step_data *t = TOOL_STEPDATA;
	t->time = T;
	t->ppar2 = ppar2;
	t->pperp2= pperp2;
	t->B = B;
	t->x = x;
	t->y = y;
	t->z = z;
	t->vx = vhatx;
	t->vy = vhaty;
	t->vz = vhatz;
	t->vpar = vpar;
	t->vperp = vperp;
}
void tool_update_position(double x, double y, double z) {
	step_data *t = TOOL_STEPDATA;
	t->x = x;
	t->y = y;
	t->z = z;
}
void tool_step(ode_solution *solver_object, double Jdtdrho) {
	TOOL_STEPDATA->Jdtdrho = Jdtdrho;
    TOOL_STEPDATA->dt = solver_object->actualstep;
	selected_tool->step(solver_object, TOOL_STEPDATA);
}

