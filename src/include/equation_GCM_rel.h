#ifndef _EQUATION_GCM_REL_H
#define _EQUATION_GCM_REL_H

#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"

void equation_GCM_rel_init(settings*);
void equation_GCM_rel_init_run(particle*, ode_solution*);
solution_data *equation_GCM_rel_init_output(solution_data*);
vector *equation_GCM_rel_eq(double, vector*);

#endif/*_EQUATION_GCM_REL_H*/
