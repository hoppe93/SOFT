#ifndef _MAGNETIC_NUM_H
#define _MAGNETIC_NUM_H

#include <stdio.h>
#include "magnetic_field.h"
#include "settings.h"
#include "vector.h"

/* Returns the magnetic field strength in a point of cartesian coordinates */
vector* magnetic_num_get(double, double, double);

void magnetic_num_init(struct general_settings*);
void magnetic_num_init_run(void);
void magnetic_num_init_particle(void);
diff_data *magnetic_num_diff(double,double,double);
diff_data *magnetic_num_diff_notor(double,double,double);

#endif/*_MAGNETIC_NUM_H*/
