#ifndef _SETTINGS_H
#define _SETTINGS_H

#include "particles.h"

/*struct run_args {
	particle *particles;/ * List of particles to run * /
	int n;				/ * Number of particles * /
};*/
struct run_args {
	struct particlespec *spec;
};

struct general_settings {
	char *name;			/* Name of sycout which these settings apply to */
	char **setting;		/* List of setting names */
	char **value;		/* List of setting values */
	int n;				/* Number of settings */
};

typedef struct {
	double tolerance;	/* RKF45 solver tolerance */
	double maxtimestep;	/* Largest allowed timestep */
	char *magfield;		/* Name of magnetic field handler */
	char *domain;		/* Domain file */
	char *equation;		/* Name of equation to solve */
	char *tool;			/* Name of tool to use */
	struct general_settings *magnetic;/* Magnetic field handler settings */
	int nmagnetic;		/* Number of mf handler settings (for different handlers) */
	struct general_settings *tools;/* Tool settings */
	int ntools;			/* Number of tool settings (for different tools) */
	struct general_settings *sycouts;/* Sycamera outputs */
	int nsycouts;		/* Number of sycamera outputs */
	particle *particles;/* Particle settings */
	int nparts;			/* Number of particles */
	struct particlespec *particlespec;/* Recipe for generating particles */
	int threads;		/* Number of threads to use */
	int warnoncollision;/* Warn when particles collide with the wall */
	double interptimestep;/* 0 = Do not linearly interpolate timesteps,
							 n = Interpolate with n points between timesteps (int),
							 x = Interpolate so that the timestep is at most x (double < 1) */
	int nodrifts;		/* If 1, drop drift terms from the GC EOM's */
} settings;

enum settings_token_type {
	TKN_INVALID,
	TKN_NAME,			/* Setting name */
	TKN_VALUE,			/* Setting value */
	TKN_MAGNETIC,		/* Magnetic keyword */
	TKN_TOOL,			/* Tool keyword */
	TKN_SYCOUT,			/* Sycout keyword */
	TKN_PARTICLE,		/* Particle keyword */
	TKN_PARTICLES,		/* Particles keyword (distinct from 'particle', which is now deprecated) */
	TKN_BLOCK_START,	/* Block start character ({) */
	TKN_BLOCK_END,		/* Block end character (}) */
	TKN_EOF				/* End-of-file marker */
};

typedef struct {
	enum settings_token_type type;
	char *val;
} settings_token;

extern char DEBUG_OUTPUT;

settings *load_settings(const char*);

#endif/*_SETTINGS_H*/
