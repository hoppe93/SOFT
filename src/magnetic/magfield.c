/* Magnetic field handler
 *
 * This module loads magnetic equilibrium
 * data and domain information in one. The input
 * should be in HDF5 format.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "magfield.h"
#include "sfile.h"

/**
 * Transpose a matrix
 */
double **magfield_transpose(double **a, int rows, int cols) {
	int i, j;
	double **r = malloc(sizeof(double*)*cols),
			*p = malloc(sizeof(double)*rows*cols);

	for (i = 0; i < cols; i++) {
		r[i] = p + i*rows;
		for (j = 0; j < rows; j++) {
			r[i][j] = a[j][i];
		}
	}

	free(a[0]);
	free(a);

	return r;
}
/**
 * Load magnetic field data from the HDF5
 * file name "filename". The file is expected
 * to contain certain datasets.
 */
magfield_t *magfield_load(const char *filename, enum sfile_type ftype) {
	double **temp;
	sfilesize_t dims[2];
	magfield_t *data;
	data = malloc(sizeof(magfield_t));
	sFILE *s;

	/* If filetype has not been specified, try
	 * to find out the file type */
	if (ftype == FILETYPE_UNKNOWN) {
		ftype = sfile_get_filetype(filename);
		if (ftype == FILETYPE_UNKNOWN) {
			fprintf(stderr, "ERROR: Unable to determine file format of magnetic field data: '%s'.\n", filename);
			exit(-1);
		}
	}

	s = sfile_init(ftype);
	s->open(s, filename, SFILE_MODE_READ);

	data->name = s->get_string(s, "name");
	data->desc = s->get_string(s, "desc");
	if (data->name == NULL) {fprintf(stderr, "The magnetic equilibrium data has no name!\n"); exit(-1);}
	if (data->desc == NULL) {fprintf(stderr, "The magnetic equilibrium data has no description!\n"); exit(-1);}

	temp = s->get_doubles(s, "maxis", dims);
	if (temp == NULL) {fprintf(stderr, "The magnetic equilibrium data contains no information about the magnetic axis (maxis)!\n"); exit(-1);}
	data->axis_r = temp[0][0];
	data->axis_z = temp[0][1];
	free(temp);

	temp = s->get_doubles(s, "r", dims);
	if (temp == NULL) {fprintf(stderr, "The magnetic equilibrium data contains no radial grid data!\n"); exit(-1);}
	data->r = *temp; free(temp);
	data->nr = dims[1];
	temp = s->get_doubles(s, "z", dims);
	if (temp == NULL) {fprintf(stderr, "The magnetic equilibrium data contains no vertical grid data!\n"); exit(-1);}
	data->z = *temp; free(temp);
	data->nz = dims[1];

	/* Load and, if necessary, transpose the magnetic
	 * field components. Different file formats tend
	 * to save the magnetic field in different ways,
	 * and so sometimes we may have to transpose them.
	 */
	data->Br = s->get_doubles(s, "Br", dims);
	if (dims[0] != (unsigned int)data->nr) data->Br = magfield_transpose(data->Br, dims[0], dims[1]);
	data->Bphi=s->get_doubles(s, "Bphi", dims);
	if (dims[0] != (unsigned int)data->nr) data->Bphi = magfield_transpose(data->Bphi, dims[0], dims[1]);
	data->Bz = s->get_doubles(s, "Bz", dims);
	if (dims[0] != (unsigned int)data->nr) data->Bz = magfield_transpose(data->Bz, dims[0], dims[1]);

	if (data->Br == NULL) {fprintf(stderr, "The magnetic equilibrium data has no radial component!\n"); exit(-1);}
	if (data->Bphi==NULL) {fprintf(stderr, "The magnetic equilibrium data has no toroidal component!\n"); exit(-1);}
	if (data->Bz == NULL) {fprintf(stderr, "The magnetic equilibrium data has no vertical component!\n"); exit(-1);}

	temp = s->get_doubles(s, "wall", dims);
	if (temp != NULL) {
		data->nwall = dims[1];
		data->wall_r = temp[0];
		data->wall_z = temp[1];
		free(temp);
	} else {
		data->nwall = 0;
		data->wall_r = NULL;
		data->wall_z = NULL;
	}

	temp = s->get_doubles(s, "separatrix", dims);
	if (temp != NULL) {
		data->nsep = dims[1];
		data->sep_r = temp[0];
		data->sep_z = temp[1];
	} else {
		data->nsep = 0;
		data->sep_r = NULL;
		data->sep_z = NULL;
	}

	s->close(s);

	return data;
}
