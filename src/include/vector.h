#ifndef _VECTOR_H
#define _VECTOR_H

/* Define the vector type */
typedef struct {
	double *val;		/* The actual vector */
	unsigned int n;		/* Number of dimensions */
} vector;

/* The following macro(s) are useful for debugging
 * memory leaks, which are not memory leaks in the
 * way valgrind defines it.
 */
//#define vnew(n) vnew_dbg(n, __FILE__, __LINE__)
#define vnew(n) vnew_int(n)
vector *vnew_int(unsigned int);
vector *vnew_dbg(unsigned int, char const*, int);
vector *vinit(unsigned int, ...);
void vfree(vector*);
vector *vadd(vector*, vector*, vector*);
vector *vaddf(vector*, vector*);
vector *vmuls(double, vector*, vector*);
vector *vmulsf(double, vector*);
void vequate(vector*, vector*);
double vdot(vector*, vector*);
double vdot3(vector*, vector*);
double vnorm(vector*);
double vnorm3(vector*);


#endif/*_VECTOR_H*/
