/* Domain manager */

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "constants.h"
#include "domain.h"
#include "magnetic_axis.h"
#include "readfile.h"

extern int errno;
double *domain_r, *domain_z;
int domain_npoints;
int domain_no_outer_wall=0;
double domain_major_radius=0.;

/**
 * Get the domain object
 */
domain *domain_get(void) {
	domain *d = malloc(sizeof(domain));
	d->name = malloc(sizeof("Unknown")+1);
	strcpy(d->name, "Unknown");
	d->n = domain_npoints;
	d->r = domain_r;
	d->z = domain_z;

	return d;
}

/**
 * Get maximum bounds of domain.
 *
 * RETURNS bound_r and bound_z
 */
void domain_get_bounds(double *rmax, double *rmin, double *zmax, double *zmin) {
	if (rmax != NULL) *rmax = domain_r[0];
	if (rmin != NULL) *rmin = domain_r[0];
	if (zmax != NULL) *zmax = domain_z[0];
	if (zmin != NULL) *zmin = domain_z[0];

	int i;
	for (i = 1; i < domain_npoints; i++) {
		if (rmax != NULL && domain_r[i] > *rmax) *rmax = domain_r[i];
		if (rmin != NULL && domain_r[i] < *rmin) *rmin = domain_r[i];
		if (zmax != NULL && domain_z[i] > *zmax) *zmax = domain_z[i];
		if (zmin != NULL && domain_z[i] < *zmin) *zmin = domain_z[i];
	}
}

/**
 * Get the major radius of the device
 */
double domain_get_major_radius(void) {
	if (domain_major_radius > 0)
		return domain_major_radius;

	int i;
	for (i = 0; i < domain_npoints; i++) {
		if (domain_major_radius < domain_r[i])
			domain_major_radius = domain_r[i];
	}

	return domain_major_radius;
}

/**
 * Read the numeric data from a domain file
 *
 * f: The C file object to read through
 * d: Domain object to load data into
 * n: Number of points to read
 *
 * Called from domain_load
 */
void domain_read_data(FILE *f, domain *d, unsigned int n) {
	unsigned int i;
	char *p;
	d->r = malloc(sizeof(double)*n);
	d->z = malloc(sizeof(double)*n);

	for (i = 0; i < n; i++) {
		/* Read R value */
		p = readfile_word(f);
		sscanf(p, "%lf", d->r+i);
		/* Read Z value */
		p = readfile_word(f);
		sscanf(p, "%lf", d->z+i);
		/* Skip flag */
		readfile_word(f);
	}
}

/**
 * Sets the 'no_outer_wall' parameter
 * which ignores the outer wall, i.e.
 * all wall tiles at R > R_maj, where
 * R_maj is the major radius of the
 * device.
 **/
void domain_remove_outer_wall(void) {
	domain_no_outer_wall = 1;
}

/**
 * Load a file containing data points giving
 * the domain of the problem.
 *
 * filename: Name of file to load domain from
 *
 * RETURNS the domain object loaded
 */
domain *domain_load(char *filename) {
	FILE *f;
	domain *d=NULL;

	f = fopen(filename, "r");
	/* File error */
	if (!f) {
		perror("ERROR");
		fprintf(stderr, "Unable to open file: %s\n", filename);
		exit(EXIT_FAILURE);
	}

	d = malloc(sizeof(domain));
	/* Memory error */
	if (d == NULL) {
		fprintf(stderr, "ERROR: Memory error!\n");
		exit(EXIT_FAILURE);
	}

	/* Copy name */
	d->name = malloc(strlen(filename)+1);
	strcpy(d->name, filename);

	/* Skip first 2 lines */
	readfile_skip_lines(2, f);
	/* Get number of points */
	char *np = readfile_word(f);
	d->n = strtol(np, (char**)NULL, 10);

	/* Skip 2 lines */
	readfile_skip_lines(2, f);

	/* Read the actual data */
	domain_read_data(f, d, d->n);

	/* DONE */
	//domain_d = d;
	return d;
}

void domain_set(double *r, double *z, int n) {
	domain_r = r;
	domain_z = z;
	domain_npoints = n;
}

/**
 * Check if a points lies within or outside of the domain.
 *
 * d: Domain for the problem
 * r: change in Radial coordinate, two points
 * z: change in Z-coordinate, two points
 *
 * RETURNS DOMAIN_WITHIN (0) or
 * DOMAIN_OUTSIDE (1) depending on whether the given point
 * is inside or outside the given domain
 */
int domain_check(double *r, double *z) {
	/* A point inside the contour */
	double x00=r[0], y00=z[0];
	double x01=r[1]-x00, y01=z[1]-y00;
	double x10, y10, x11, y11, s, t, det;

	int i;
	for (i=0;i < domain_npoints-1; i++) {
		/* Ignore wall segments which are obviously not intersected */
		/*if ((d->z[i] < z[0] && d->z[i+1] < z[0] && d->z[i] < z[1] && d->z[i+1] < z[1]) ||
			(d->z[i] > z[0] && d->z[i+1] > z[0] && d->z[i] > z[1] && d->z[i+1] > z[1]))
			return DOMAIN_WITHIN;*/

		x10=domain_r[i];
		x11=domain_r[i+1]-x10;
		y10=domain_z[i];
		y11=domain_z[i+1]-y10;

		/* Check if matrix is zero */
		if (x00-x10==0 && y00-y10==0) return DOMAIN_OUTSIDE;

		/* Calculates the determinant */
		det=x11*y01-x01*y11;
		/* Check if determinant is zero */
		if (det==0) return DOMAIN_WITHIN;

		/* Calculates s and t */
		s=(1/det)*((x00-x10)*y01-(y00-y10)*x01);
		t=(1/det)*(-(-(x00-x10)*y11+(y00-y10)*x11));
		/* If s and t are between 0 and 1 => intersection */
		if (s>=0 && s<=1 && t>=0 && t<=1) return DOMAIN_OUTSIDE;
	}

	return DOMAIN_WITHIN;
}

/**
 * Check whether the given line (given by its start- and end-points)
 * passes through the outside of the device.
 *
 * d: Domain object
 * xp: X coordinate of particle
 * yp: Y coordinate of particle
 * zp: Z coordinate of particle
 * xc: X coordinate of camera
 * yc: Y coordinate of camera
 * zc: Z coordinate of camera
 */
int domain_check3d(
	double xp, double yp, double zp,
	double xc, double yc, double zc
) {
	int retval = DOMAIN_WITHIN;
	int i;
	for (i=0; i < domain_npoints-1; i++) {
		double dr=domain_r[i], drp=domain_r[i+1];
		double dz=domain_z[i], dzp=domain_z[i+1];

		/* Ignore all outer wall tiles (so the camera can
		 * be placed outside the device) */
		if (domain_no_outer_wall &&
			dr > magnetic_axis_r && drp > magnetic_axis_r) {
			//return DOMAIN_WITHIN;
			continue;
		}

		if ((dz < zp && dzp < zp && dz < zc && dzp < zc) ||
			(dz > zp && dzp > zp && dz > zc && dzp > zc))
			continue;

		if (domain_check3d_one(
				drp, dzp, dr, dz,
				xp, yp, zp,
				xc, yc, zc) == DOMAIN_OUTSIDE) {
			retval = DOMAIN_OUTSIDE;
		}
	}

	return retval;
}
int domain_check3d_one(
	double xnp, double znp, double xn, double zn,
	double xp, double yp, double zp,
	double xc, double yc, double zc
) {
	double t1=0., t2=0.,
		   s1=0., s2=0.;
	if (znp == zn) {
		if (zc == zp) return DOMAIN_WITHIN;		/* Equation becomes singular => line of sight does not cross wall section */
		t1 = t2 = (zn-zp)/(zc-zp);
		double ic = 1/(xn-xnp);
		double term1 = xp + (xc-xp)*t1, term2 = yp + (yc-yp)*t1;
		double square = term1*term1+term2*term2;
		double sqr = sqrt(square);
		s1 = ic*(xn + sqr);
		s2 = ic*(xn - sqr);

		if ((s1 >= 0. && s1 <= 1.) || (s2 >= 0. && s2 <= 1.))
			return DOMAIN_OUTSIDE;
		else return DOMAIN_WITHIN;
	} else {
		double izn=1/(znp-zn);
		double kz=(zc-zp)*izn;
		double mz=(zp-zn)*izn;
		double a=xc-xp, b=yc-yp, c=(xnp-xn)*kz;
		double m=xn+mz*(xnp-xn);

		double abc = a*a + b*b - c*c;
		if (abc == 0.) return DOMAIN_OUTSIDE;
		double iabc = 1/abc;

		double term = (c*m - a*xp - b*yp)*iabc;
		double q = (m*m-xp*xp-yp*yp)*iabc;
		double square = term*term + q;
		if (square < 0.) return DOMAIN_WITHIN;	/* no real solutions */
		double sqr = sqrt(square);
		t1 = term + sqr, t2 = term - sqr;
		s1 = kz*t1, 	 s2 = kz*t2;
		if ((t1 >= 0. && t1 <= 1.) || (t2 >= 0. && t2 <= 1.))
			return DOMAIN_OUTSIDE;
	}

	return DOMAIN_WITHIN;
}

/*
void domain_test3d(void) {
	domain_load("d3d.wall_2d");

	#define POINTS 7
	const int points = POINTS;
	double testpoints[POINTS][4] = {
		/ * x      y       z        seen * /
		{0.8802,  1.4920, 0.,      (double)(DOMAIN_OUTSIDE)},	// 1
		{1.1230, -1.3260, 0.,      (double)(DOMAIN_WITHIN)},	// 2
		{1.5450,  0.4549, 0.,      (double)(DOMAIN_WITHIN)},	// 3
		{0.9397,  1.4020, 0.,      (double)(DOMAIN_OUTSIDE)},	// 4
		{0.8686,  1.4540, 0.,      (double)(DOMAIN_OUTSIDE)},	// 5
		{1.4780,  0.8701, 0.,      (double)(DOMAIN_OUTSIDE)},	// 6
		{8.965684e-01,1.343647e+00,1.273679e-01,(double)(DOMAIN_OUTSIDE)}	// 7
	};

	int i, fail=0;
	for (i = 0; i < points; i++) {
		if (domain_check3d(
			testpoints[i][0], testpoints[i][1], testpoints[i][2],
			0., -2.20, 0.) != (int)testpoints[i][3]) {
			printf("Point %d marked as %s though it should be %s!\n",(i+1),(testpoints[i][3]==(double)DOMAIN_OUTSIDE?"detected":"invisible"),(testpoints[i][3]==(double)DOMAIN_OUTSIDE?"invisible":"detected"));
			fail = 1;
		}
	}

	if (!fail) printf("All points we're detected correctly!\n");
}
*/
/**
 * Function for testing the module
 */
/*
void domain_test(void) {
	/ * Check wether reading file is OK * /
	domain *d = domain_load("iter.wall_2d");

	srand(time(NULL));
	int i = rand() % (d->n);

	printf("Number of points: %d\n", d->n);
	printf("First point:      i=0 , r=%f, z=%f\n", d->r[0], d->z[0]);
	printf("Random point:     i=%2d, r=%f, z=%f\n", i, d->r[i], d->z[i]);
	printf("Last point:       i=%2d, r=%f, z=%f\n", d->n-1, d->r[d->n-1], d->z[d->n-1]);

	/ * Test points * /
	double r[2];
	double z[2];

	/ * To convert from int to certain message * /
	char *location[]={"Inside","Outside"};
	/ * Indicating the result * /
	int is;
	/ * Indicating what the result should be* /
	int should;

	/ * TEST BEGINS * /
	printf(" ********* TEST BEGINS ********\n");

	/ * TESTPOINT 1 * /
	r[0]=4.79839;r[1]=6.92137;
	z[0]=1.78125; z[1]=0.4375;
	should=0;

	is=domain_check(d, r, z);

	printf("Should be %s, is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * TESTPOINT 2 * /
	r[0]=6.03226;r[1]=6.47681;
	z[0]=-4.8125;z[1]=-0.78125;

	should=1;
	is=domain_check(d, r, z);

	printf("Should be %s, is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * TESTPOINT 3 * /
	r[0]=6.07762;r[1]=6.64012;
	z[0]=4.125;z[1]=3.65625;
	should=0;
	is=domain_check(d, r, z);

	printf("Should be %s is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * TESTPOINT 4 * /
	r[0]=7.91028;r[1]=7.21169;
	z[0]=-1.46875;z[1]=-1.78125;
	should=1;
	is=domain_check(d, r, z);

	printf("Should be %s is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * TESTPOINT 5 * /
	r[0]=4.16331;r[1]=4.3629;
	z[1]=-1.96875;z[1]=-3.75;
	should=1;
	is=domain_check(d, r, z);

	printf("Should be %s is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * TESTPOINT 6 * /
	r[0]=5.44254; r[1]=5.59677;
	z[0]=-4.03125;z[1]=-4.71875;
	should=1;
	is=domain_check(d, r, z);

	printf("Should be %s is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * TESTPOINT 7 * /
	r[0]=4.39012;r[1]=4.57157;
	z[0]=-2.3125;z[1]=3.3125;
	should=0;
	is=domain_check(d, r, z);

	printf("Should be %s is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * TESTPOINT 8 * /
	r[0]=7.88717;r[1]=7.99016;
	z[0]=-1.34802;z[1]=-1.16703;
	should=1;
	is=domain_check(d, r, z);

	printf("Should be %s is %s\n",location[should],location[is]);
	if (should!=is) {printf("INCORRECT!!!\n"); return;}
	printf("CORRECT!!!\n\n");

	/ * END OF TEST * /
	printf(" ********* END OF TEST ********\n");
	printf("TEST DONE, WELL DONE!\n");
}
*/
