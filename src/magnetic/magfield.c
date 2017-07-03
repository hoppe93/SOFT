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

//char *(*magfield_get_string)(hid_t, const char*);
//double **(*magfield_get_doubles)(hid_t, const char*, hsize_t);

/**
 * Reads a string from the dataset with name "name".
 *
 * fileid: HDF5 file descriptor of file to read
 * name: Name of dataset to load string from
 */
/*
char *magfield_get_string(hid_t fileid, const char *name) {
	char *s;
	size_t length;
	hid_t filetype, memtype, dset;

	dset = H5Dopen(fileid, name, H5P_DEFAULT);
	filetype = H5Dget_type(dset);
	length = H5Tget_size(filetype);
	memtype = H5Tcopy(H5T_C_S1);
	H5Tset_size(memtype, length);

	s = malloc(sizeof(char)*(length+1));
	H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, s);
	s[length] = 0;

	H5Tclose(filetype);
	H5Tclose(memtype);
	H5Dclose(dset);

	return s;
}
*/
/**
 * Reads an array of C doubles from the dataset
 * with name "name". The dimensions of the array are
 * returned in "dims".
 *
 * fileid: HDF5 file descriptor to read from
 * name: Name of dataset to read
 * dims: Array with two elements. Will contain dimensions of returned value on return.
 *
 * Returns a 1-D array (logically 2-D) which contains
 * the data of the dataset.
 */
/*
double **magfield_get_doubles(hid_t fileid, const char *name, hsize_t *dims) {
	double *data, **pointers;
	hid_t dset, space;
	int ndims;
	unsigned int i;

	if (H5Lexists(fileid, name, H5P_DEFAULT) <= 0)
		return NULL;
	
	dset = H5Dopen(fileid, name, H5P_DEFAULT);

	space = H5Dget_space(dset);
	ndims = H5Sget_simple_extent_dims(space, dims, NULL);

	if (dims == NULL) {
		data = malloc(sizeof(double)*ndims);
		pointers = malloc(sizeof(double*));
		pointers[0] = data;
	} else if (ndims == 1) {
		data = malloc(sizeof(double)*dims[0]);
		pointers = malloc(sizeof(double*));
		pointers[0] = data;
	} else {
		data = malloc(sizeof(double)*dims[0]*dims[1]);
		pointers = malloc(sizeof(double*)*dims[0]);
		for (i = 0; i < dims[0]; i++) {
			pointers[i] = data+(i*dims[1]);
		}
	}
	
	H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Sclose(space);
	H5Dclose(dset);

	return pointers;
}
*/

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
	//(*sfile_open)(filename, SFILE_MODE_READ);
	(*(s->open))(s, filename, SFILE_MODE_READ);

	data->name = (*(s->get_string))(s, "name");
	data->desc = (*(s->get_string))(s, "desc");
	if (data->name == NULL) {fprintf(stderr, "The magnetic equilibrium data has no name!\n"); exit(-1);}
	if (data->desc == NULL) {fprintf(stderr, "The magnetic equilibrium data has no description!\n"); exit(-1);}

	temp = (*(s->get_doubles))(s, "maxis", dims);
	if (temp == NULL) {fprintf(stderr, "The magnetic equilibrium data contains no information about the magnetic axis (maxis)!\n"); exit(-1);}
	data->axis_r = temp[0][0];
	data->axis_z = temp[0][1];
	free(temp);

	temp = (*(s->get_doubles))(s, "r", dims);
	if (temp == NULL) {fprintf(stderr, "The magnetic equilibrium data contains no radial grid data!\n"); exit(-1);}
	data->r = *temp; free(temp);
	data->nr = dims[1];
	temp = (*(s->get_doubles))(s, "z", dims);
	if (temp == NULL) {fprintf(stderr, "The magnetic equilibrium data contains no vertical grid data!\n"); exit(-1);}
	data->z = *temp; free(temp);
	data->nz = dims[1];

	data->Br = (*(s->get_doubles))(s, "Br", dims);
	data->Bphi=(*(s->get_doubles))(s, "Bphi", dims);
	data->Bz = (*(s->get_doubles))(s, "Bz", dims);
	if (data->Br == NULL) {fprintf(stderr, "The magnetic equilibrium data has no radial component!\n"); exit(-1);}
	if (data->Bphi==NULL) {fprintf(stderr, "The magnetic equilibrium data has no toroidal component!\n"); exit(-1);}
	if (data->Bz == NULL) {fprintf(stderr, "The magnetic equilibrium data has no vertical component!\n"); exit(-1);}

	temp = (*(s->get_doubles))(s, "wall", dims);
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

	temp = (*(s->get_doubles))(s, "separatrix", dims);
	if (temp != NULL) {
		data->nsep = dims[1];
		data->sep_r = temp[0];
		data->sep_z = temp[1];
	} else {
		data->nsep = 0;
		data->sep_r = NULL;
		data->sep_z = NULL;
	}

	(*(s->close))(s);

	return data;
}
