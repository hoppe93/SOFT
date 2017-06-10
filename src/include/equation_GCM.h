#ifndef _EQUATION_GCM_H
#define _EQUATION_GCM_H

#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"

void equation_GCM_init(settings*);
void equation_GCM_init_run(particle*, ode_solution*);
solution_data *equation_GCM_init_output(solution_data*);
vector *equation_GCM_eq(double, vector*);

#endif/*_EQUATION_GCM_H*/
