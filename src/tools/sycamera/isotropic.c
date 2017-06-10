/* Uniform emission from spheres */

#include <omp.h>
#include <stdlib.h>
#include "sycamera.h"

void isotropic_init(
	enum sycamera_radiation_type radt, enum sycamera_polarization_type polt,
	double *lambdas, int spectrum_resolution, int integral_resolution
) {}

void isotropic_init_run(void) {}
void isotropic_init_particle(particle *p) {}
void isotropic_init_step(step_data *sd) {}
double isotropic_intensity(
	step_data *sd, vector *rcp, vector *vhat,
	vector *empty1, vector *empty2, vector *empty3
) {
	double rcpl2 = rcp->val[0]*rcp->val[0] + rcp->val[1]*rcp->val[1] + rcp->val[2]*rcp->val[2];
	return 1.0/rcpl2;
}

double *isotropic_get_wavelengths(void) {
	return NULL;
}
double *isotropic_get_spectrum(void) {
	return NULL;
}
