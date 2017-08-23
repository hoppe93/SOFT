/* Uniform emission from spheres */

#include <omp.h>
#include <stdlib.h>
#include "sycamera.h"

void sphere_init(
	enum sycamera_radiation_type radt, enum sycamera_polarization_type polt,
	double *lambdas, int spectrum_resolution, int integral_resolution
) {}

void sphere_init_run(void) {}
void sphere_init_particle(particle *p) {}
void sphere_init_step(step_data *sd) {}
double sphere_intensity(
	step_data *sd, vector *rcp, vector *vhat,
	vector *empty1, vector *empty2, vector *empty3
) {
	return 1.0;
}

double *sphere_get_wavelengths(void) {
	return NULL;
}
double *sphere_get_spectrum(void) {
	return NULL;
}
int sphere_get_spectrum_length(void) {
	return 0;
}
