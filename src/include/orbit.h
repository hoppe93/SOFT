#ifndef _ORBITS_H
#define _ORBITS_H

#include "arguments.h"
#include "equations.h"
#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"
#include "tools.h"

void orbit_init(settings*, struct general_settings*, struct general_settings*, int);
void orbit_init_run(unsigned int);
ode_solution *orbit_init_particle(particle*);
void orbit_deinit_run(void);
void orbit_step(ode_solution*, step_data*);
int orbit_stop_condition(void);
void orbit_output(equation*);
void orbit_allocate(void);

#endif/*_ORBITS_H*/
