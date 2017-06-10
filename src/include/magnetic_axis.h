#ifndef _MAGNETIC_AXIS_H
#define _MAGNETIC_AXIS_H
/* Routine for finding magnetic axis */

double *magnetic_axis_find(double);
void magnetic_axis_set_loc(double, double);
void magnetic_axis_test(void);

extern double magnetic_axis_r, magnetic_axis_z;
extern int magnetic_axis_set;

#endif/*_MAGNETIC_AXIS_H*/
