/* Equation handler
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "magnetic_field.h"
#include "vector.h"

/**
 * Allocate a string and copy the given string
 * to it
 */
char *setname(const char *n) {
	char *s = malloc(strlen(n)+1);
	strcpy(s, n);
	return s;
}
/* ASCII to double pointer (n-component) */
double *atodpn(const char *s, int n, double *dp) {
	char buffer[1024];
	int bi=0, di=0;

	if (dp == NULL)
		dp = malloc(sizeof(double)*n);

	/* Allow for easier syntax below */
	s--;

	/* Turn comma-separated list of numbers into
	* actual double list */
	while (*(++s) && di < n) {
		if (*s == ',') {
			buffer[bi] = 0;
			dp[di++] = atof(buffer);
			bi = 0;
		} else buffer[bi++] = *s;
	}

	/* Handle last number */
	if (di < n) {
		buffer[bi] = 0;
		dp[di] = atof(buffer);
	}

	return dp;
}

/**
 * Convert velocity given as v_|| and v_\perp
 * to a 3D velocity component.
 */
double *vpp2v(double vpar, double vperp, double *xyz, double *v) {
	/* Get bhat in xyz */
	double x=xyz[0], y=xyz[1], z=xyz[2];
	vector *bhat = magnetic_field_get(x,y,z);
	double Babs = sqrt(vdot(bhat,bhat));

	bhat->val[0] /= Babs;
	bhat->val[1] /= Babs;
	bhat->val[2] /= Babs;

	double bhatox, bhatoy, bhatoz;
	if (bhat->val[2] > bhat->val[1] && bhat->val[2] > bhat->val[0]) {
		bhatoz = -(bhat->val[0] + bhat->val[1])/bhat->val[2];
		bhatoy = 1;
		bhatox = 1;
	} else if (bhat->val[1] > bhat->val[0]) {
		bhatox = 1;
		bhatoy = -(bhat->val[2]+bhat->val[0])/bhat->val[1];
		bhatoz = 1;
	} else {
		bhatox = 1;
		bhatoy = 1;
		bhatoz = -(bhat->val[0]+bhat->val[1])/bhat->val[2];
	}

	/* Normalize orthogonal vector */
	double bhato_abs = sqrt(bhatox*bhatox + bhatoy*bhatoy + bhatoz*bhatoz);
	bhatox = bhatox/bhato_abs;
	bhatoy = bhatoy/bhato_abs;
	bhatoz = bhatoz/bhato_abs;

	/* Calculate velocity */
	if (v == NULL) v = malloc(sizeof(double)*3);

	v[0] = vpar*bhat->val[0] + vperp*bhatox;
	v[1] = vpar*bhat->val[1] + vperp*bhatoy;
	v[2] = vpar*bhat->val[2] + vperp*bhatoz;

	return v;
}

/**
 * Convert velocity given as vx, vy, vz
 * into a v_||/v_\perp component pair
 *
 * Returns a double array with 2 elements.
 * The first is vpar, the second is vperp.
 */
double *v32vpp(double vx, double vy, double vz, double *xyz) {
	double x=xyz[0], y=xyz[1], z=xyz[2];
	vector *bhat = magnetic_field_get(x,y,z);
	double Babs = sqrt(vdot(bhat,bhat));

	bhat->val[0] /= Babs;
	bhat->val[1] /= Babs;
	bhat->val[2] /= Babs;

	double v2 = vx*vx + vy*vy + vz*vz;
	double *v = malloc(sizeof(double)*2);
	v[0] = vx*bhat->val[0] + vy*bhat->val[1] + vz*bhat->val[2];
	v[1] = sqrt(v2 - v[0]*v[0]);

	return v;
}

void print_timediff(int mpi_rank, int i, int pos) {
	static struct timespec debug_clock;
	static int initialized=0;

	struct timespec newclock;

	if (!initialized) {
		clock_gettime(CLOCK_REALTIME, &debug_clock);
		initialized = 1;
	} else {
		if (!clock_gettime(CLOCK_REALTIME, &newclock)) {
			long ds = newclock.tv_sec - debug_clock.tv_sec;
			long dn = newclock.tv_nsec - debug_clock.tv_nsec;
			if (dn < 0) { ds -= 1; dn += 1000000000; }

			printf("[%d]: %d(%d) -- %ld seconds, %ld ns\n", mpi_rank, i, pos, ds, dn);
			debug_clock = newclock;
		} else
			printf("[%d]: %d(%d) -- Unable to read CLOCK_REALTIME\n", mpi_rank, i, pos);
	}
}

