#ifndef _INTERP2_H
#define _INTERP2_H

#include "magnetic_num.h"
#include "magfield.h"
#include "vector.h"

vector* interp2_interpolate(double, double, double);
/* Must be used before interpolation! */
void interp2_init_interpolation(magfield_t*);

double **interp2_jacobian(double, double, double);
double diff_BrDr(double);

#endif/*_INTERP2_H*/
