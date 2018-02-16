#ifndef _TOOLS_H
#define _TOOLS_H

#include "arguments.h"
#include "equations.h"
#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"

typedef struct {
	double time;
	double ppar2;
	double pperp2;
	double B;
	double x, y, z;
	double vx, vy, vz;		/* Guiding-center velocity vector */
	double vpar, vperp;		/* Particle velocity components */
	double Jdtdrho;			/* Jacobian * timestep * radial step */
    double dt;              /* Time-step */
} step_data;

typedef struct {
	char *name;						/* Tool name */
	double tolerance;				/* Recommended ODE Solver accuracy */
	int require_jacobian;			/* 1 = Requires the Jacobian connecting (R, Z) and (t, rho)
									   to be computed ni each timestep Two orbits must then be solved */
	void (*init)(struct general_settings*,struct general_settings*,int);	/* Initialization function */
	void (*init_run)(unsigned int);
	ode_solution *(*init_particle)(particle *p);
	void (*step)(ode_solution*, step_data*);
	void (*deinit_run)(void);
	int (*stop_condition)(void);	/* Returns true if simulation should be terminated */
	void (*output)(equation*);			/* Generate output function */
} tool;

void tool_init_handler(void);
void tool_prepare_run(void);
tool *tool_select(char*);
/*void tool_in_step(double,
				  double, double, double,
				  double, double, double,
				  double, double, double);*/
void tool_update_position(double, double, double);
void tool_record_data(double, double, double, double,
					  double, double, double,
					  double, double, double,
					  double, double);
void tool_step(ode_solution*, double);

#endif/*_TOOLS_H*/
