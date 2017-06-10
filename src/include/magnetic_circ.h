#ifndef _MAGNETIC_CIRC_H
#define _MAGNETIC_CIRC_H
/* Circular magnetic field handler */

#include "magnetic_field.h"
#include "settings.h"

/* Circular magnetic field parameters */
struct magnetic_circ_params {
	double B0;		/* Magnetic field strength */
	double Rm;		/* Major radius */
	double r;		/* Minor radius */
	double q;		/* Safety factor (constant) */
};

void magnetic_circ_init(struct general_settings*);
void magnetic_circ_init_run(void);
void magnetic_circ_init_particle(void);
vector *magnetic_circ_eval(double, double, double);
diff_data *magnetic_circ_diff(double, double, double);
diff_data *magnetic_circ_diff_notor(double, double, double);

#endif/*_MAGNETIC_CIRC_H*/
