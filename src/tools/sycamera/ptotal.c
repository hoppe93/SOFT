/* This function determines which spectrum weighting
 * function that should be used.
 */

#include <omp.h>
#include <stdio.h>
#include "sycamera.h"
#include "tools.h"

void sycamera_spectrum_init(double lambda0, double lambda1, int res) {
    /* We initialize all methods */
    //sycamera_pas2_init(lambda0, lambda1, res);
    sycamera_pcyl_init(lambda0, lambda1, res);
}
void sycamera_spectrum_init_run(void) {
    sycamera_pcyl_init_run();
}
double sycamera_spectrum_weight(step_data *sd, double mass) {
/*
    if (sycamera_pas2_valid(sd->ppar2, sd->pperp2, sd->B, mass))
        return sycamera_pas2_int();
    else*/ return sycamera_pcyl_int(sd->ppar2, sd->pperp2, sd->B, mass);
}
double *sycamera_spectrum_get_wavelengths(void) {
	return sycamera_pcyl_get_wavelengths();
}
double *sycamera_spectrum_get(void) {
/*
    if (sycamera_pas2_valid(sd->ppar2, sd->pperp2, sd->B, mass))
        return sycamera_pas2_int();
    else*/ return sycamera_pcyl_get_spectrum();
}
int sycamera_spectrum_length(void) {
	return sycamera_pcyl_get_spectrum_length();
}
