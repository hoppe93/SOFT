#ifndef _DOMAIN_H
#define _DOMAIN_H

/* Return values of `domain_check' */
#define DOMAIN_WITHIN 0
#define DOMAIN_OUTSIDE 1

/**
 * The `domain' data type contains the
 * coordinates for the points that make
 * up the walls of the reactor.
 */
typedef struct {
	char *name;
	double *r, *z;		/* Coordinates */
	int n;				/* Number of points */
} domain;

domain *domain_get(void);
void domain_get_bounds(double*, double*, double*, double*);
double domain_get_major_radius(void);
/* Load the domain coordinates from file */
domain *domain_load(char*); // *
void domain_set(double*, double*, int);
/* Check if the given points lies within or
 * outside the domain
 */
int domain_check(double* , double* );
int domain_check3d(double, double, double, double, double, double);
int domain_check3d_one(double, double, double, double, double, double, double, double, double, double);
void domain_remove_outer_wall(void);

/* Functions for testing the module */
void domain_test3d(void);
void domain_test(void);

#endif/*_DOMAIN_H*/
