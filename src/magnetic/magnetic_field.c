/* Magnetic field handler */

#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "magnetic_circ.h"
#include "magnetic_field.h"
#include "magnetic_num.h"
#include "util.h"

const int NUMBER_OF_HANDLERS=2;
magnetic_handler *all_handlers=NULL, *selected_handler=NULL;

void magnetic_init(void) {
	all_handlers = malloc(sizeof(magnetic_handler)*NUMBER_OF_HANDLERS);

	all_handlers[0].name = setname("numeric");
	all_handlers[0].init = magnetic_num_init;
	all_handlers[0].init_run = magnetic_num_init_run;
	all_handlers[0].init_particle = magnetic_num_init_particle;
	all_handlers[0].eval = magnetic_num_get;
	all_handlers[0].diff = magnetic_num_diff;
	all_handlers[0].diff_notor = magnetic_num_diff_notor;

	all_handlers[1].name = setname("circular");
	all_handlers[1].init = magnetic_circ_init;
	all_handlers[1].init_run = magnetic_circ_init_run;
	all_handlers[1].init_particle = magnetic_num_init_particle;
	all_handlers[1].eval = magnetic_circ_eval;
	all_handlers[1].diff = magnetic_circ_diff;
	all_handlers[1].diff_notor = magnetic_circ_diff_notor;
}

magnetic_handler *magnetic_handler_select(char *name) {
	int i;
	for (i = 0; i < NUMBER_OF_HANDLERS; i++) {
		if (!strcmp(name, all_handlers[i].name)) {
			printf("Selecting magnetic field handler '%s'...\n", name);
			selected_handler = all_handlers+i;
			return selected_handler;
		}
	}

	return NULL;
}

vector *magnetic_field_get(double x, double y, double z) {
	return selected_handler->eval(x,y,z);
}
diff_data *magnetic_field_diff(double x, double y, double z) {
	return selected_handler->diff(x, y, z);
}
diff_data *magnetic_field_diff_notor(double x, double y, double z) {
	return selected_handler->diff_notor(x, y, z);
}
magnetic_handler *magnetic_handler_get(void) {
	return selected_handler;
}
