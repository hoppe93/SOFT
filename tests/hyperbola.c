/* Testing hyperbola calculation */

#include <math.h>
#include <stdio.h>
#include <sycamera.h>
#include "test.h"

void cone_delta_hyperbola_intersections(
	double, double, double, double, double,
	double*, double*, double*, double*,
	int*, int*, double, double
);
void cone_delta_find_intersections(
	double*, double*, double*, int*,
	double, double, double, double, double,
	double, double, void (*intersect)(
		double, double, double, double, double,
		double*, double*, double*, double*,
		int*, int*, double, double
	)
);

void test_hyp_single_solution(int branch, double a, double b, double ola, double x0, double y0, double rdet, double xi) {
	double Ax,Bx,Ay,By,d,c,t1,t2,x1=0,x2=0,y1=0,y2=0;
	int ok1=0, ok2=0;

	switch (branch) {
		case 1:
			Ax = a*cos(xi), Bx = b*sin(xi);
			Ay = a*sin(xi), By = b*cos(xi);
			c = y0;
			d = x0-rdet/2.0;
			break;
		case 2:
			Ax = a*cos(xi), Bx = b*sin(xi);
			Ay = a*sin(xi), By = b*cos(xi);
			c = y0;
			d = x0+rdet/2.0;
			break;
		case 3:
			Ax =-a*sin(xi), Bx = b*cos(xi);
			Ay =-a*cos(xi), By = b*sin(xi);
			c = x0;
			d = y0-rdet/2.0;
			break;
		case 4:
			Ax =-a*sin(xi), Bx = b*cos(xi);
			Ay =-a*cos(xi), By = b*sin(xi);
			c = x0;
			d = y0+rdet/2.0;
			break;
		default:
			fprintf(stderr, "ERROR: Invalid branch chosen for single hyperbola test: %d.\n", branch);
			return;
	}

	cone_delta_hyperbola_intersections(Ax, Bx, Ay, By, d, &t1, &t2, &y1, &y2, &ok1, &ok2, c, rdet/2.0);

	printf("ok1 = %d,  ok2 = %d\n", ok1, ok2);
	if (ok1) {
		printf("t1 = %e,  x1 = %e,  y1 = %e\n", t1, x1, y1);
	}
	if (ok2) {
		printf("t2 = %e,  x2 = %e,  y2 = %e\n", t2, x2, y2);
	}
}
void test_hyp_multi(double a, double b, double ola, double x0, double y0, double xi, double rdet) {
	double t[8], xval[8], yval[8];
	int tindex=0, i;
	
	cone_delta_find_intersections(t, xval, yval, &tindex, x0, y0, a, b, rdet/2.0, cos(xi), sin(xi), &cone_delta_hyperbola_intersections);

	printf("Solutions found = %d\n", tindex);
	printf("Solution parameters =");
	for (i=0;i<tindex;i++) printf("  %e", t[i]);
	printf("\nSolution values (x) =");
	for (i=0;i<tindex;i++) printf("  %e", xval[i]);
	printf("\nSolution values (y) =");
	for (i=0;i<tindex;i++) printf("  %e", yval[i]);
	printf("\n");
}
int test_hyp(void) {
	double rdet = 0.2,		/* Detector aperture */
		   thetap = 0.5,	/* Pitch angle */
		   phi = 1.1,		/* Inclination between velocity and camera normal */
		   X = .12,			/* Distance from focus point to particle */
		   xi = 0;			/* Rotation of curve in camera plane */
	
	double cms = cos(thetap)*cos(thetap) - sin(phi)*sin(phi),
		a = fabs((cos(thetap)*sin(thetap)*X*cos(phi))/cms),
		b = fabs(sin(thetap)*X*cos(phi))/sqrt(fabs(cms)),
		ola = -sin(thetap)*sin(thetap)*X*sin(phi) / cms,
		x0 = -ola*cos(xi),
		y0 = ola*sin(xi);
	
	if (cms >= 0) {
		fprintf(stderr, "ERROR: The specified parameters give 'cms' = %e >= 0.\n", cms);
		return 2;
	}

	test_hyp_multi(a, b, ola, x0, y0, xi, rdet);
	//test_hyp_single_solution(3, a, b, ola, x0, y0, rdet, xi);

	return 0;
}
