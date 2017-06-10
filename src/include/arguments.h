#ifndef _ARGUMENTS_H
#define _ARGUMENTS_H

typedef struct {
	double tstart, tend,	/* Start and end times */
		   vpar, vperp;		/* Initial parallel/perpendicular velocities */
	double *r0, *v0;		/* Initial position and velocity */
	char *magfield_file,	/* Name of magnetic field file */
		 *domain_file,		/* Name of domain file */
		 *output_file;		/* Name of output file */
	double particle_mass;	/* Mass of particle */
	double particle_charge;	/* Charge of particle */
	char print_settings,	/* Wether or not to print settings at top
							   of output file */
		 vpar_set,			/* 1 if vpar has been set, 0 otherwise */
		 vperp_set;			/* 1 if vperp has been set, 0 otherwise */
	char *equation;
	double *detector_direction;/* Detector direction */
	double *detector_position;/* Detector X position */
	double detector_radius;	/* Detector radius */
	int detector_pixels;	/* Number of pixels of detector (per dimension, true number is pixels^2) */
	double detector_vision_angle;/* Detector viewing angle */
	char *tool;				/* Name of tool to use */
} arguments;

#define PROBLEM_GC 1		/* Solve the guiding center problem */
#define PROBLEM_NO_GC 0 	/* Solve the particle motion problem */

arguments *arguments_default(void);
arguments *parse_args(int, char*[]);

#endif/*_ARGUMENTS_H*/
