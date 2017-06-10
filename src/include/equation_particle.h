#ifndef _EQUATION_PARTICLE_H
#define _EQUATION_PARTICLE_H

#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"
#include "vector.h"

void equation_particle_init(settings*);
/* Save input particle in local variable in equation */
void equation_particle_init_run(particle*, ode_solution*);

/* Initialize output structure */
solution_data *equation_particle_init_output(solution_data*);

/* Equation for the charged particle motion. Lorentz force, magnetic field only. */
vector * equation_particle_eq(double, vector* );

#endif/*_EQUATION_PARTICLE_H*/
