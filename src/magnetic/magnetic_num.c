/* Magnetic field reader */

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "domain.h"
#include "interp2.h"
#include "magnetic_axis.h"
#include "magnetic_field.h"
#include "magfield.h"
#include "magnetic_num.h"
#include "vector.h"
#include "particles.h"
#include "readfile.h"
#include "sfile.h"
#include "util.h"

magfield_t *magnetic_num_field;
enum magfield_wall_type magnetic_num_wall=MAGFIELD_ANY;

vector *magnetic_num_gradB, *magnetic_num_curlB;
diff_data *magnetic_num_diff_dd;
#pragma omp threadprivate(magnetic_num_gradB,magnetic_num_curlB,magnetic_num_diff_dd)

void magnetic_num_init(struct general_settings *set) {
	double *axis = NULL;
	char *magfield_name = NULL;
	magnetic_num_field = NULL;
	enum sfile_type fformat = FILETYPE_UNKNOWN;

	if (set == NULL) {
		fprintf(stderr, "ERROR: The numeric magnetic field handler requires settings to be given!\n");
		exit(-1);
	}

	/* Loop over all settings */
	int i;
	for (i = 0; i < set->n; i++) {
		if (!strcmp(set->setting[i], "file")) {
			magfield_name = set->value[i];
		} else if (!strcmp(set->setting[i], "wall")) {
			if (!strcmp(set->value[i], "separatrix"))
				magnetic_num_wall = MAGFIELD_SEPARATRIX;
			else if (!strcmp(set->value[i], "wall"))
				magnetic_num_wall = MAGFIELD_WALL;
			else if (!strcmp(set->value[i], "any"))
				magnetic_num_wall = MAGFIELD_ANY;
			else {
				fprintf(stderr, "ERROR: Unrecognized wall type specified for magnetic field: '%s'\n", set->setting[i]);
				exit(-1);
			}
		} else if (!strcmp(set->setting[i], "axis")) {
			axis = atodpn(set->value[i], 2, NULL);
			if (axis[0] <= 0.0) {
				fprintf(stderr, "ERROR: The magnetic axis radial location must be greater than zero.\n");
				exit(EXIT_FAILURE);
			}
		} else if (!strcmp(set->setting[i], "format")) {
			if (!strcmp(set->value[i], "auto"))
				fformat = FILETYPE_UNKNOWN;
			else if (!strcmp(set->value[i], "hdf5"))
				fformat = FILETYPE_HDF5;
			else if (!strcmp(set->value[i], "mat"))
				fformat = FILETYPE_MATLAB;
			else {
				fprintf(stderr, "ERROR: Unrecognized magnetic field format specified: '%s'.\n", set->value[i]);
				exit(-1);
			}
		} else {
			fprintf(stderr, "Invalid magnetic field handler 'numeric' setting: %s!", set->setting[i]);
			exit(-1);
		}
	}

	/* Was name of magnetic field given? */
	if (magfield_name == NULL) {
		fprintf(stderr, "No magnetic field has been specified! Unable to continue.\n");
		exit(-1);
	}

	/* Load magnetic field */
	magnetic_num_field = magfield_load(magfield_name, fformat);

	/* Successfully loaded? */
	if (magnetic_num_field == NULL) {
		fprintf(stderr, "ERROR: Unable to load magnetic field data: '%s'.\n", magfield_name);
		exit(-1);
	}

	/* Set magnetic axis manually if requested */
	if (axis != NULL) {
		magnetic_axis_set_loc(axis[0], axis[1]);
		free(axis);
	} else {
		/* Set magnetic axis */
		magnetic_axis_set_loc(magnetic_num_field->axis_r, magnetic_num_field->axis_z);
	}

	/* Select domain */
	if (magnetic_num_wall == MAGFIELD_WALL) {
		if (magnetic_num_field->wall_r == NULL) {
			fprintf(stderr, "The magnetic field provided does not contain the device wall.\n");
			exit(-1);
		}

		domain_set(magnetic_num_field->wall_r, magnetic_num_field->wall_z, magnetic_num_field->nwall);
	} else if (magnetic_num_wall == MAGFIELD_SEPARATRIX) {
		if (magnetic_num_field->sep_r == NULL) {
			fprintf(stderr, "The magnetic field provided does not contain the separatrix.\n");
			exit(-1);
		}

		domain_set(magnetic_num_field->sep_r, magnetic_num_field->sep_z, magnetic_num_field->nwall);
	} else {
		if (magnetic_num_field->wall_r != NULL)
			domain_set(magnetic_num_field->wall_r, magnetic_num_field->wall_z, magnetic_num_field->nwall);
		else if (magnetic_num_field->sep_r != NULL)
			domain_set(magnetic_num_field->sep_r, magnetic_num_field->sep_z, magnetic_num_field->nwall);
		else {
			fprintf(stderr, "The magnetic field provided contains neither the wall nor separatrix.\n");
			exit(-1);
		}
	}
}
void magnetic_num_init_run(void) {
	interp2_init_interpolation(magnetic_num_field);

	magnetic_num_gradB = vnew(3);
	magnetic_num_curlB = vnew(3);
	magnetic_num_diff_dd = malloc(sizeof(diff_data));

/*
	if (!magnetic_axis_set)
		magnetic_axis_find((magnetic_num_field->rmin + magnetic_num_field->rmax)/2);
*/
}
void magnetic_num_init_particle(void) {}

/**
 * Calculates the magnetic field strength in a given point (x,y,z).
 *
 * B: The magnetic field
 * xyz: The point (in cartesian coordinates) in which the field
 *   strength should be evaluated.
 *
 * RETURNS the field strength at the given point in
 * cartesian coordinates
 */
vector* magnetic_num_get(double x, double y, double z) {
	vector *B_interp = interp2_interpolate(x,y,z);
    return B_interp;
}

diff_data *magnetic_num_diff(double x, double y, double z) {
	/* Calculate Jacobian */
	double **J = interp2_jacobian(x,y,z);
	double
		dBx_dx=J[0][0],  dBx_dy=J[0][1],  dBx_dz=J[0][2],
		dBy_dx=J[1][0],  dBy_dy=J[1][1],  dBy_dz=J[1][2],
		dBz_dx=J[2][0],  dBz_dy=J[2][1],  dBz_dz=J[2][2];

	vector *B = magnetic_num_get(x,y,z);
	/* Quicker access to the components of B */
	double Bx = B->val[0], By = B->val[1], Bz = B->val[2];
	/* Absolute value of B, ||B|| */
	double Babs = sqrt(Bx*Bx + By*By + Bz*Bz);

	/* Calculate components of gradient (grad ||B||) */
	double
		gradBx = (Bx*dBx_dx + By*dBy_dx + Bz*dBz_dx) / Babs,
		gradBy = (Bx*dBx_dy + By*dBy_dy + Bz*dBz_dy) / Babs,
		gradBz = (Bx*dBx_dz + By*dBy_dz + Bz*dBz_dz) / Babs;

	/* Create gradient vector */
	magnetic_num_gradB->val[0] = gradBx;
	magnetic_num_gradB->val[1] = gradBy;
	magnetic_num_gradB->val[2] = gradBz;

	/* Calculate curl bhat */
	double
		curlBx = (Babs*dBz_dy - Bz*gradBy - Babs*dBy_dz + By*gradBz) / (Babs*Babs),
		curlBy = (Babs*dBx_dz - Bx*gradBz - Babs*dBz_dx + Bz*gradBx) / (Babs*Babs),
		curlBz = (Babs*dBy_dx - By*gradBx - Babs*dBx_dy + Bx*gradBy) / (Babs*Babs);

	magnetic_num_curlB->val[0] = curlBx;
	magnetic_num_curlB->val[1] = curlBy;
	magnetic_num_curlB->val[2] = curlBz;

	magnetic_num_diff_dd->B = B;
	magnetic_num_diff_dd->gradB = magnetic_num_gradB;
	magnetic_num_diff_dd->curlB = magnetic_num_curlB;
	magnetic_num_diff_dd->Babs = Babs;

	return magnetic_num_diff_dd;
}
diff_data *magnetic_num_diff_notor(double x, double y, double z) {
	/* Calculate Jacobian (we set By = 0) */
	double **J = interp2_jacobian(x,y,z);
	double
		dBx_dx=J[0][0],  dBx_dy=J[0][1],  dBx_dz=J[0][2],
		dBy_dx=0.,       dBy_dy=0.,       dBy_dz=0.,
		dBz_dx=J[2][0],  dBz_dy=J[2][1],  dBz_dz=J[2][2];

	vector *B = magnetic_num_get(x,y,z);
	/* Quicker access to the components of B (set By = 0) */
	double Bx = B->val[0], By = 0., Bz = B->val[2];
	/* Absolute value of B, ||B|| */
	double Babs = sqrt(Bx*Bx + By*By + Bz*Bz);

	/* Calculate components of gradient (grad ||B||) */
	double
		gradBx = (Bx*dBx_dx + By*dBy_dx + Bz*dBz_dx) / Babs,
		gradBy = (Bx*dBx_dy + By*dBy_dy + Bz*dBz_dy) / Babs,
		gradBz = (Bx*dBx_dz + By*dBy_dz + Bz*dBz_dz) / Babs;

	/* Create gradient vector */
	magnetic_num_gradB->val[0] = gradBx;
	magnetic_num_gradB->val[1] = gradBy;
	magnetic_num_gradB->val[2] = gradBz;

	/* Calculate curl bhat */
	double
		curlBx = (Babs*dBz_dy - Bz*gradBy - Babs*dBy_dz + By*gradBz) / (Babs*Babs),
		curlBy = (Babs*dBx_dz - Bx*gradBz - Babs*dBz_dx + Bz*gradBx) / (Babs*Babs),
		curlBz = (Babs*dBy_dx - By*gradBx - Babs*dBx_dy + Bx*gradBy) / (Babs*Babs);

	magnetic_num_curlB->val[0] = curlBx;
	magnetic_num_curlB->val[1] = curlBy;
	magnetic_num_curlB->val[2] = curlBz;

	magnetic_num_diff_dd->B = B;
	magnetic_num_diff_dd->gradB = magnetic_num_gradB;
	magnetic_num_diff_dd->curlB = magnetic_num_curlB;
	magnetic_num_diff_dd->Babs = Babs;

	return magnetic_num_diff_dd;
}

