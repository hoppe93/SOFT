#ifndef _MAGFIELD_H
#define _MAGFIELD_H

#include "sfile.h"

typedef struct {
	char *name;					/* Name of magnetic field */
	char *desc;					/* Optional description of field */
	double axis_r, axis_z;		/* Location of magnetic axis */
	int nr, nz;					/* Number of grid points in r and z */
	double *r, *z;				/* r and z points (1-D arrays) */
	double **Br, **Bphi, **Bz;	/* Magnetic field data (2-D arrays of size nr x nz ) */
	int nwall;					/* Number of wall points */
	double *wall_r, *wall_z;	/* Wall r and z coordinates */
	int nsep;					/* Number of separatrix points */
	double *sep_r, *sep_z;		/* Separatrix r and z coordinates */
} magfield_t;

enum magfield_wall_type {
	MAGFIELD_WALL,
	MAGFIELD_SEPARATRIX,
	MAGFIELD_ANY
};

magfield_t *magfield_load(const char*, enum sfile_type);

#endif/*_MAGFIELD_H*/
