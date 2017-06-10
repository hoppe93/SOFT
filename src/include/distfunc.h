#ifndef _DISTFUNC_H
#define _DISTFUNC_H

#include <stdio.h>
#include <stdlib.h>

typedef struct {
	double **value;		/* Values of the distribution function */
	double *r;			/* List of radial coordinates */
	double rmin, rmax;	/* Smallest/largest value of radial coordinate */
	double *p;			/* List of momentum */
	double pmin, pmax;	/* Smallest/largest value of momentum */
	double *costheta;	/* List of pitch angles */
	double costmin, costmax;/* Minimum/maximum cosine of pitch angle */
	int nr, np, ncostheta;/* Number of radial values, number of momenum values, number of pitch angle values */
} distfunc;

#define DISTFUNC_BUFFER_SIZE 1024
void distfunc_init_run(void);
void distfunc_load(const char*);
double distfunc_eval(double, double, double);

/* Readfile helpers */
/* Macro for checking that the most recently read value was not empty */
#define BUFCHECK(b) if (*b == 0) {fprintf(stderr, "ERROR: Badly formatted line in map!\n"); exit(EXIT_FAILURE);}

char df_readfile_stopped_at_newline(void);
int df_readfile_eof(void);
void df_readfile_reset(void);
void df_readfile_loaddbl(double*);
void df_readfile_loadint(int*);
/* Read one word from file (returns empty string if empty line) */
char *df_readfile_value(void);
/* Skip the given number of lines */
void df_readfile_skip_lines(int);
void df_readfile_unload(void);
void df_readfile_load(const char*);

/* Interpolation routines */
void df_interp_init(distfunc*);
void df_interp_init_run(void);
double df_interp_eval(double, double, double);
void df_interp_error_handler(const char*, const char*, int, int);

#endif/*_DISTFUNC_H*/
