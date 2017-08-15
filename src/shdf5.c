/**
 * A SOFT HDF5 interface for simplified I/O of
 * HDF5 files.
 */

#include <hdf5.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "sfile.h"
#include "shdf5.h"

hid_t fileid;
#pragma omp threadprivate(fileid)

/**
 * Open HDF5 file.
 *
 * filename: Name of file to open
 */
int shdf5_open(sFILE *s, const char *filename, enum sfile_mode mode) {
	s->identifier = malloc(sizeof(hid_t));
	s->mode = mode;

	switch (mode) {
		case SFILE_MODE_READ:
			*((hid_t*)(s->identifier)) = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
			break;
		case SFILE_MODE_UPDATE:
			*((hid_t*)(s->identifier)) = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
			break;
		case SFILE_MODE_WRITE:
			*((hid_t*)(s->identifier)) = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			break;
		default:
			fprintf(stderr, "Unrecognized option for opening HDF5 file: %d.\n", mode);
			free(s->identifier);
			s->identifier = NULL;
			return 0;
	}

	if (*((hid_t*)(s->identifier)) < 0) {
		fprintf(stderr, "Unable to open HDF5 file: %s\n", filename);
		return 0;
	}

	return 1;
}

/**
 * Close the currently open HDF5 file.
 */
void shdf5_close(sFILE *s) {
	H5Fclose(*((hid_t*)(s->identifier)));
}

/*******************************
 ************ INPUT ************
 *******************************/
/**
 * Reads a string from the dataset with name "name".
 *
 * name: Name of dataset to load string from
 */
char *shdf5_get_string(sFILE *s, const char *name) {
	char *str;
	size_t length;
	hid_t filetype, memtype, dset;

	dset = H5Dopen(*((hid_t*)(s->identifier)), name, H5P_DEFAULT);
	filetype = H5Dget_type(dset);
	length = H5Tget_size(filetype);
	memtype = H5Tcopy(H5T_C_S1);
	H5Tset_size(memtype, length);

	str = malloc(sizeof(char)*(length+1));
	H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, str);
	str[length] = 0;

	H5Tclose(filetype);
	H5Tclose(memtype);
	H5Dclose(dset);

	return str;
}
/**
 * Reads an array of C doubles from the dataset
 * with name "name". The dimensions of the array are
 * returned in "dims".
 *
 * name: Name of dataset to read
 * dims: Array with two elements. Will contain dimensions of returned value on return.
 *
 * Returns a 1-D array (logically 2-D) which contains
 * the data of the dataset. If the named dataset does not exist,
 * NULL is returned.
 */
double **shdf5_get_doubles(sFILE *s, const char *name, sfilesize_t *dims) {
	double *data, **pointers;
	hid_t dset, space;
	sfilesize_t ndims;
	size_t i;
	hid_t fileid = *((hid_t*)(s->identifier));

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

/*******************************
 *********** OUTPUT ************
 *******************************/
/**
 * Write a dataset consisting of a string
 *
 * fileid: hid_t id of file to write to
 * name: Name of dataset to create
 * str: String to write
 * length: Length of string to write, or <= 0
 *   to automatically determine the length.
 */
void shdf5_write_string(sFILE *s, const char *name, const char *str, size_t length) {
	hid_t dsetid, filetype, memtype, spaceid;

	hsize_t dim[1] = {1};
	filetype = H5Tcopy(H5T_C_S1);
	memtype = H5Tcopy(H5T_C_S1);

	if (length <= 0)
		length = strlen(str);

	H5Tset_size(filetype, length); 
	H5Tset_size(memtype, length);
	
	spaceid = H5Screate_simple(1, dim, NULL);
	dsetid = H5Dcreate(*((hid_t*)s->identifier), name, filetype, spaceid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dsetid, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, str);

	H5Dclose(dsetid);
	H5Sclose(spaceid);
	H5Tclose(filetype);
	H5Tclose(memtype);
}

/**
 * Write a 2-D array to the target file.
 *
 * fileid: hid_t id of the HDF5 file
 * name: Name of the dataset to create
 * arr: 2-D array data to write
 * rows: Number of rows of array
 * cols: Number of columns of data
 */
void shdf5_write_array(sFILE *s, const char *name, double **arr, size_t rows, size_t cols) {
	hid_t spaceid, dsetid;

	hsize_t dims[2] = {rows, cols};
	spaceid = H5Screate_simple(2, dims, NULL);

	dsetid = H5Dcreate2(*((hid_t*)(s->identifier)), name, H5T_IEEE_F64LE, spaceid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);

	H5Dclose(dsetid);
	H5Sclose(spaceid);
}
/**
 * Write an image to the HDF5 file. This
 * function is just a wrapper for
 * 'shdf5_write_array()'.
 *
 * fileid: hid_t id of file to write to
 * name: Name of dataset to create
 * image: Image data to write
 * n: Number of pixels to write (image is assumed square)
 */
void shdf5_write_image(sFILE *s, const char *name, double **image, size_t n) {
	shdf5_write_array(s, name, image, n, n);
}
/**
 * Write a simple list (1-D array) of data
 * to the HDF5 file. This function is just
 * a wrapper for 'shdf5_write_array()'.
 *
 * fileid: hid_t id of file to write to
 * name: Name of dataset to create
 * list: Data to write
 * n: Number of elements in list
 */
void shdf5_write_list(sFILE *s, const char *name, double *list, size_t n) {
	shdf5_write_array(s, name, &list, 1, n);
}

/**
 * Add a scalar attribute to a HDF5 dataset.
 *
 * fileid: hid_t id of file to write to
 * dsetname: Name of dataset to apply this attribute to
 * name: Name of the attribute to create
 * val: Value to give the attribute
 */
void shdf5_write_attribute_scalar(sFILE *s, const char *dsetname, const char *name, double val) {
	hid_t dsetid, attspaceid, attid;

	dsetid = H5Dopen2(*((hid_t*)(s->identifier)), dsetname, H5P_DEFAULT);
	attspaceid = H5Screate(H5S_SCALAR);

	attid = H5Acreate2(dsetid, name, H5T_IEEE_F64LE, attspaceid, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attid, H5T_NATIVE_DOUBLE, &val);

	H5Aclose(attid);
	H5Sclose(attspaceid);
	H5Dclose(dsetid);
}

/**
 * Add a string attribute to a HDF5 dataset.
 *
 * fileid: hid_t id of file to write to
 * dsetname: Name of dataset to apply this attribute to
 * name: Name of the attribute to create
 * str: String to write to the attribute
 * length: Length of string to write, or <= 0
 *   to automatically determine length
 */
void shdf5_write_attribute_string(sFILE *s, const char *dsetname, const char *name, const char *str, size_t length) {
	hid_t dsetid, filetype, memtype, spaceid, attid;

	hsize_t dim[1] = {1};
	filetype = H5Tcopy(H5T_C_S1);
	memtype = H5Tcopy(H5T_C_S1);

	if (length <= 0)
		length = strlen(str);

	H5Tset_size(filetype, length); 
	H5Tset_size(memtype, length);
	
	dsetid = H5Dopen2(*((hid_t*)(s->identifier)), dsetname, H5P_DEFAULT);
	spaceid = H5Screate_simple(1, dim, NULL);
	attid = H5Acreate(dsetid, name, filetype, spaceid, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attid, memtype, str);

	H5Aclose(attid);
	H5Dclose(dsetid);
	H5Sclose(spaceid);
	H5Tclose(memtype);
	H5Tclose(filetype);
}

