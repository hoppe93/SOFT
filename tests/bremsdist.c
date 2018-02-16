/* Testing bremsstrahlung spectrum calculation */

#include <math.h>
#include <stdio.h>
#include <sycamera.h>
#include "test.h"

void sycamera_bremsdist_evalZ(double, double);

int test_bsdist(void) {
	double
		E0 = 37.658246484419998,
		theta0 = 0.0,
		pitch = 0.052006199537590,
		Zeff = 1.0,
		k1 = 0.5,
		k2 = 0.5,

		/* Computed quantities */
		mu = theta0+pitch,
		E02 = E0*E0,
		p0 = sqrt(E02 - 1),
		beta = p0 / E0,
		sinmu = sin(mu), cosmu = cos(mu),
		sinp = sin(pitch), cosp = cos(pitch);
    mu = 1.584621531520000;
	
	sycamera_zeff = Zeff;
	sycamera_bremsdist_init(k1, k2, 100);
	sycamera_bremsdist_init_run();
	sycamera_bremsdist_init_particle(ELECTRONMASS, p0*ELECTRONMASS*LIGHTSPEED);

	//double val = sycamera_bremsdist_int(E0, E02, beta, sinmu, cosmu, sinp, cosp);
	double val = sycamera_bremsdist_int(E0, E02, beta, sin(mu), cos(mu), sinp, cosp);
	printf("%e\n", val);

    sycamera_bremsdist_evalZ(0.051252121601571, 0.026925466310344);

	return 0;
}
