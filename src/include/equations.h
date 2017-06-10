#ifndef _EQUATIONS_H
#define _EQUATIONS_H

#include "equation_predprey.h"
#include "equation_particle.h"
#include "equation_GCM.h"
#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"
#include "vector.h"

typedef struct {
	char *name;
	unsigned int variables;
	vector *(*eq)(double, vector*);
	void (*init)(settings*);
	void (*init_run)(particle*,ode_solution*);
	solution_data *(*output)(solution_data*);
} equation;

/* Initialize equation handler */
void equation_handler_init();

equation *select_equation(char*);
//ode_solution *init_equation(equation*, vector*, initial_data*);

#endif/*_EQUATIONS_H*/
