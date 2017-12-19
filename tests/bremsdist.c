/* Testing bremsstrahlung spectrum calculation */

#include <math.h>
#include <stdio.h>
#include <sycamera.h>
#include "test.h"

int test_bsdist(void) {
	double
		E0 = 9.0,
		theta0 = 0.05,
		pitch = 0.13,
		Zeff = 1.0,
		k1 = 5.0,
		k2 = 10.0,

		/* Computed quantities */
		mu = theta0+pitch,
		E02 = E0*E0,
		p0 = sqrt(E02 - 1),
		beta = p0 / E0,
		sinmu = sin(mu), cosmu = cos(mu),
		sinp = sin(pitch), cosp = cos(pitch);
	
	sycamera_zeff = Zeff;
	sycamera_bremsdist_init(k1, k2, 100);
	sycamera_bremsdist_init_run();
	sycamera_bremsdist_init_particle(ELECTRONMASS, p0);

	double val = sycamera_bremsdist_int(E0, E02, beta, sinmu, cosmu, sinp, cosp);
	printf("%e\n", val);

	return 0;
}
