/* Vector operations */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>	/* For obtaining backtrace */
#include "vector.h"

/*
 * Create a new n-dimensional vector
 *
 * n: Number of dimensions
 * RETURNS: A new n-dimensional vector
 */
vector *vnew_int(unsigned int n) {
	/* Declare vector */
	vector *v;
	v = malloc(sizeof(vector));

	/* Set size of vector */
	v->n = n;

	/* Allocate memory for matrix */
	v->val = malloc(sizeof(double)*n);
	if (v->val == NULL) printf("MEMERR\n");

	return v;
}
vector *vnew_dbg(unsigned int n, char const *caller_name, int line) {
	fprintf(stderr, "vnew called from %s:%d\n", caller_name, line);
	return vnew_int(n);
}
/*
 * Initialize vector with values from a given set of double's
 *
 * n: Number of vector dimensions
 * VARIABLE ARGUMENTS: Vector elements. All must be double's!!
 * RETURNS: A new matrix containing the same elements
 * as v
 */
vector *vinit(unsigned int n, ...) {
	/* Create a variable argument list */
	va_list valist;
	/* Allocate new vector */
	vector *v = vnew(n);

	/* Initialize VA list */
	va_start(valist, n);

	/* Copy arguments to vector in order */
	unsigned int i;
	for (i = 0; i < n; i++) {
		double val = va_arg(valist, double);
		v->val[i] = val;
	}

	/* Clean memory */
	va_end(valist);

	return v;
}
/*
 * Free the memory occupied by a vector
 *
 * vec: Vector to free
 */
void vfree(vector *vec) {
	free(vec->val);
	free(vec);
}
/*
 * Add two vectors
 *
 * a: First vector to add
 * b: Second vector to add
 * RETURNS: A new vector with
 * elements equal to a+b
 */
vector *vadd(vector *a, vector *b, vector *result) {
	/* Make sure a and b are of the same size */
	if (a->n != b->n) {
		printf("ERROR: Given vectors are of different size!\n");
		exit(EXIT_FAILURE);
	}

	/* Check sum matrix */
	if (result->n != a->n) {
		printf("ERROR: The result vector is of a different size than the term vectors!\n");
		exit(EXIT_FAILURE);
	}

	/* Sum over all elements */
	unsigned int i;
	for (i = 0; i < a->n; i++) {
		result->val[i] = a->val[i] + b->val[i];
	}

	return result;
}

/*
 * Add two vectors, but DO NOT allocate
 * space for a new vector. Store result in `a'.
 *
 * a: First vector to add
 * b: Second vector to add
 * RETURNS: A new vector with
 * elements equal to a+b
 */
vector *vaddf(vector *a, vector *b) {
	/* Make sure a and b are of the same size */
	if (a->n != b->n) {
		printf("ERROR: Given vectors are of different size!\n");
		exit(EXIT_FAILURE);
	}

	/* Sum over all elements */
	unsigned int i;
	for (i = 0; i < a->n; i++) {
		a->val[i] += b->val[i];
	}

	return a;
}

/*
 * Multiply each element of a vector
 * by a scalar.
 *
 * scalar: Scalar to multiply with
 * a: Vector to multiply with
 */
vector *vmuls(double scalar, vector *a, vector *result) {
	if (result->n != a->n) {
		printf("ERROR: Result vector and a are of different size!\n");
		exit(EXIT_FAILURE);
	}

	unsigned int i;
	for (i = 0; i < a->n; i++) {
		result->val[i] = scalar * a->val[i];
	}

	return result;
}

/*
 * Multiply each element of a vector
 * by a scalar, but DO NOT allocate
 * new memory for resulting vector.
 *
 * scalar: Scalar to multiply with
 * a: Vector to multiply with
 */
vector *vmulsf(double scalar, vector *a) {
	unsigned int i;
	for (i = 0; i < a->n; i++) {
		a->val[i] = scalar * a->val[i];
	}

	return a;
}

/*
 * Caclulates scalar product between vectors.
 * u: Vector
 * v: Vector
 * Must be the same size
 */

/*
 * Caclulates scalar product between vectors.
 * u: Vector
 * v: Vector
 * Must be the same size
 */
double vdot(vector *v, vector *u) {
	/* Check if same size */
	double scalar=0;

	unsigned int i;
	for (i = 0; i < v->n; i++) {
		scalar = scalar + v->val[i]*u->val[i];
	}

	return scalar;
}
double vdot3(vector *v, vector *u) {
	return (v->val[0]*u->val[0] + v->val[1]*u->val[1] + v->val[2]*u->val[2]);
}

/**
 * Computes the 2-norm of the given vector
 */
double vnorm(vector *v) {
	return sqrt(vdot(v,v));
}
double vnorm3(vector *v) {
	return sqrt(vdot3(v,v));
}

/**
 * Copy the values of vector b to
 * the values of vector a
 */
void vequate(vector *a, vector *b) {
	if (a->n != b->n) {
		printf("ERROR: Vectors a and b must be same length! (%d != %d)\n", a->n, b->n);
		exit(EXIT_FAILURE);
	}

	unsigned int i;
	for (i = 0; i < a->n; i++) {
		a->val[i] = b->val[i];
	}
}

/**
 * Clone a vector
 *
 * v: Vector to clone
 */
vector *vclone(vector *v) {
	vector *p = vnew(v->n);
	unsigned int i;
	for (i = 0; i < v->n; i++) {
		p->val[i] = v->val[i];
	}

	return p;
}
