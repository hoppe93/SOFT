/* Quantity check */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "global.h"
#include "quantities.h"

/* Contains interesting quantities (such as energy, xi, etc.) */
quantity *quantities;
int quantities_no;

#pragma omp threadprivate(quantities,quantities_no)

void quantities_init(void) {
	quantities = NULL;
	quantities_no = 0;
}

quantity *quantities_get(void) {
	return quantities;
}
int quantities_get_no(void) {
	return quantities_no;
}

/**
 * Define a new quantity that will
 * be reported during simulation.
 *
 * name: Name of the quantity (as
 *  displayed in output file)
 * RETURNS the index which must be
 * passed as a parameter when
 * reporting the quantity.
 */
int quantities_define(char *name) {
	int index = quantities_no;

	quantities_no++;
	quantities = realloc(quantities, sizeof(quantity)*quantities_no);
	quantities[index].values = malloc(sizeof(double)*NUMBER_OF_SIMULATION_POINTS);
	quantities[index].name = malloc(strlen(name)+1);
	strcpy(quantities[index].name, name);

	return index;
}

/**
 * Append an interesting quantity to
 * the `quantities' array.
 * 
 * quantity: Index of the quantity in 'quantities' array
 * value: The value of the quantity
 */
void quantities_report(int quant, double value, unsigned int time_index) {
	quantities[quant].values[time_index] = value;
}
void quantities_expand(int newsize) {
	int i;
	for (i = 0; i < quantities_no; i++) {
		quantities[i].values = realloc(quantities[i].values, sizeof(double)*newsize);
	}
}

/**
 * Clear all quantities
 */
void quantities_clear(void) {
	free(quantities);
}
