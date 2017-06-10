/* Equation handler
 *
 * When you have implemented a new equation it
 * should be registered in this file to allow the
 * solver to use it. Simply add four lines of code
 * to the 'equation_handler_init()' function
 * (it should be fairly obvious how from the
 * already implemented equations). Also, don't
 * forget to change the value of the global
 * variable 'NUMBER_OF_EQUATIONS' to the
 * total number of equations available. It's
 * declared at the top of the document.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "equations.h"
#include "equation_GCM.h"
#include "equation_GCM_rel.h"
#include "equation_particle.h"
#include "equation_particle_rel.h"
#include "equation_predprey.h"
#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"
#include "util.h"
#include "vector.h"

/** DON'T FORGET TO CHANGE THIS!!! */
const int NUMBER_OF_EQUATIONS=4;
equation *all_equations;

/**
 * CREATE AN ENTRY FOR YOUR EQUATION IN
 * THE FUNCTION BELOW
 **/
void equation_handler_init(void) {
	all_equations = malloc(sizeof(equation)*NUMBER_OF_EQUATIONS);

	all_equations[0].name 		= setname("guiding-center");
	all_equations[0].variables	= 5;
	all_equations[0].eq 		= equation_GCM_eq;
	all_equations[0].init 		= equation_GCM_init;
	all_equations[0].init_run	= equation_GCM_init_run;
	all_equations[0].output		= equation_GCM_init_output;

	all_equations[1].name		= setname("guiding-center-relativistic");
	all_equations[1].variables	= 5;
	all_equations[1].eq			= equation_GCM_rel_eq;
	all_equations[1].init		= equation_GCM_rel_init;
	all_equations[1].init_run	= equation_GCM_rel_init_run;
	all_equations[1].output		= equation_GCM_rel_init_output;

	all_equations[2].name 		= setname("particle");
	all_equations[2].variables	= 6;
	all_equations[2].eq 		= equation_particle_eq;
	all_equations[2].init		= equation_particle_init;
	all_equations[2].init_run	= equation_particle_init_run;
	all_equations[2].output		= equation_particle_init_output;

	all_equations[3].name		= setname("particle-relativistic");
	all_equations[3].variables	=6;
	all_equations[3].eq			= equation_particle_rel_eq;
	all_equations[3].init		= equation_particle_rel_init;
	all_equations[3].init_run	= equation_particle_rel_init_run;
	all_equations[3].output		= equation_particle_rel_init_output;
}

/**
 * NO NEED TO CARE ABOUT WHAT'S HERE
 **/
equation *select_equation(char *name) {
	int i;
	for (i = 0; i < NUMBER_OF_EQUATIONS; i++) {
		if (!strcmp(name, all_equations[i].name)) {
			printf("Selecting equation '%s'...\n", all_equations[i].name);
			return all_equations+i;
		}
	}

	return NULL;
}
