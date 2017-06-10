#ifndef _MAGNETIC_FIELD_H
#define _MAGNETIC_FIELD_H

#include "vector.h"
#include "settings.h"

typedef struct {
	vector *B;		/* Magnetic field in point */
	vector *gradB;	/* Gradient of absolute value of magnetic field */
	vector *curlB;	/* Curl of unit vector b^ (!!! not of full magnetic field !!!) */
	double Babs;	/* Absolute value of B, ||B|| */
} diff_data;

typedef struct {
	char *name;				/* Name of magnetic field handler */
	void (*init)(struct general_settings*);
	void (*init_run)(void);
	void (*init_particle)(void);
	vector* (*eval)(double, double, double);		/* Evaluate the field in the given point (cartesian coordinates) */
	diff_data* (*diff)(double, double, double);		/* Differentiate field in the given point */
	diff_data* (*diff_notor)(double, double, double);/* Differentiate field in the given point, but ignore field toroidal component */
} magnetic_handler;

void magnetic_init(void);
magnetic_handler *magnetic_handler_select(char*);
vector *magnetic_field_get(double, double, double);
diff_data *magnetic_field_diff(double, double, double);
diff_data *magnetic_field_diff_notor(double, double, double);
magnetic_handler *magnetic_handler_get(void);		/* Get current magnetic field handler */

#endif/*_MAGNETIC_FIELD_H*/
