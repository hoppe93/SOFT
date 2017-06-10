#ifndef _IO_DATA_H
#define _IO_DATA_H

#include "vector.h"
#include "quantities.h"

/**
 * Structure defining the type solution_data, containing results from a 
 * run of the integrator.
 */
typedef struct {
	double* T; 				// Time points
	vector *v;				// List of coordinates to output 
	char **labels;			//List of labels 
	unsigned int points;	// Number of points 
	unsigned int nvars;		// Number of coordinate-variables (or columns) 

	quantity *quantities;	// Other interesting quantities
	unsigned int no_quantities;// Number of quantities
} solution_data;

#endif/*_IO_DATA_H*/
