/**
 * A SOFT MAT interface for simplified I/O of
 * MATLAB .MAT files.
 */

#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "mat.h"
#include "sfile.h"
#include "smat.h"

MATFile *smat_mfp;

#pragma omp threadprivate(smat_mfp)

/**
 * Closes the currently open MATLAB file.
 */
void smat_close(sFILE *s) {
	matClose((MATFile*)s->identifier);
}

/**
 * Opens a MATLAB '.mat' file for read/write.
 *
 * filename: Name of file to open
 * ot: Purpose for opening file (read or write)
 */
int smat_open(sFILE *s, const char *filename, enum sfile_mode mode) {
	s->mode = mode;
	s->identifier = NULL;

	switch (mode) {
		case SFILE_MODE_READ:
			s->identifier = matOpen(filename, "r");
			break;
		case SFILE_MODE_UPDATE:
			s->identifier = matOpen(filename, "u");
			break;
		case SFILE_MODE_WRITE:
			s->identifier = matOpen(filename, "w");
			break;
		default:
			fprintf(stderr, "Unrecognized option for opening MATLAB file: %d.\n", mode);
			return 0;
	}

	if (s->identifier == NULL) {
		fprintf(stderr, "Unable to create or open '%s'.\n", filename);
		return 0;
	}

	return 1;
}

/******************************
 ************ INPUT ***********
 ******************************/
/**
 * Loads a variable from the MAT file as
 * a string.
 *
 * name: Name of variable to load
 *
 * Returns the variable as a C string. If
 * the string does not exist, or if the
 * variable cannot be interpreted as a
 * string, NULL is returned.
 */
char *smat_get_string(sFILE *s, const char *name) {
	MATFile *mfp = (MATFile*)s->identifier;
	mxArray *arr = matGetVariable(mfp, name);

	if (arr == NULL || mxIsEmpty(arr)) return NULL;
	if (mxGetClassID(arr) != mxCHAR_CLASS)
		return NULL;
	
	int length = mxGetN(arr)+1;
	char *str = malloc(length);
	mxGetString(arr, str, length);
	
	mxDestroyArray(arr);

	return str;
}

/**
 * Reads an array of C doubles from the MAT file.
 * The dimensions of the array are returned in
 * "dims".
 *
 * name: Name of variable to read
 * dims: Array with two elements. Contains the
 *  number of rows and columns in the returned
 *  array on return.
 *
 * RETURNS a 1-D array (logically 2-D) which contains
 * the data of the variable. If the named variable
 * does not exist, or is not a 2-D matrix, NULL
 * is returned.
 */
double **smat_get_doubles(sFILE *s, const char *name, sfilesize_t *dims) {
	int i, j, rows, cols;
	MATFile *mfp = (MATFile*)s->identifier;
	mxArray *arr = matGetVariable(mfp, name);

	if (arr == NULL || mxIsEmpty(arr)) return NULL;
	if (!mxIsDouble(arr))
		return NULL;
	
	rows = mxGetM(arr);
	cols = mxGetN(arr);

	if (dims != NULL) {
		dims[0] = rows;
		dims[1] = cols;
	}

	double **data = malloc(sizeof(double*)*rows);
	data[0] = malloc(sizeof(double)*rows*cols);

	for (i = 1; i < rows; i++) {
		data[i] = data[i-1] + cols;
	}

	double *tdata = mxGetPr(arr);
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			data[i][j] = tdata[j*rows + i];
		}
	}

	mxDestroyArray(arr);
	return data;
}

/******************************
 *********** OUTPUT ***********
 ******************************/
void smat_write_string(sFILE *s, const char *name, const char *str, int length) {
	int status;
	mxArray *ms;
	MATFile *mfp = (MATFile*)s->identifier;

	ms = mxCreateString(str);
	if (ms == NULL) {
		fprintf(stderr, "ERROR: Unable to create MATLAB string '%s'.\n", name);
		return;
	}
	status = matPutVariable(mfp, name, ms);
	mxDestroyArray(ms);

	if (status != 0)
		fprintf(stderr, "ERROR: Unable to write string '%s' to MATLAB file.\n", name);
}
void smat_write_array(sFILE *s, const char *name, double **arr, int rows, int cols) {
	int status, i, j;
	double *t;
	mxArray *ma;
	MATFile *mfp = (MATFile*)s->identifier;

	ma = mxCreateDoubleMatrix(rows, cols, mxREAL);
	if (ma == NULL) {
		fprintf(stderr, "ERROR: Unable to allocate MATLAB array for '%s'.\n", name);
		return;
	}

	t = mxGetPr(ma);
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			t[i*cols + j] = arr[i][j];
		}
	}

	status = matPutVariable(mfp, name, ma);
	mxDestroyArray(ma);
	if (status != 0) {
		fprintf(stderr, "ERROR: Unable to write variable '%s' to MATLAB file.\n", name);
	}
}
void smat_write_image(sFILE *s, const char *name, double **image, int n) {
	smat_write_array(s, name, image, n, n);
}
void smat_write_list(sFILE *s, const char *name, double *list, int n) {
	smat_write_array(s, name, &list, 1, n);
}

char *_smat_get_attribute_name(const char *dsetname, const char *name) {
	int l1 = strlen(dsetname), l2 = strlen(name);
	char *nname = malloc(sizeof(char)*(l1+l2+2));

	strcpy(nname, dsetname);
	nname[l1] = '_';
	strcpy(nname+l1+1, name);
	
	return nname;
}

/**
 * Since MATLAB files don't have attributes for
 * variables, we instead give attributes names
 * according to 'dsetname_name'.
 *
 * dsetname: Dataset name
 * name: Name of attribute
 * q: Value of scalar to write
 */
void smat_write_attribute_scalar(sFILE *s, const char *dsetname, const char *name, double q) {
	char *nname = _smat_get_attribute_name(dsetname, name);
	smat_write_list(s, nname, &q, 1);
	free(nname);
}
/**
 * Since MATLAB files don't have attributes for
 * variables, we instead give attributes names
 * according to 'dsetname_name'.
 *
 * dsetname: Dataset name
 * name: Name of attribute
 * str: Value of attribute
 * len: Length of attribute string
 */
void smat_write_attribute_string(sFILE *s, const char *dsetname, const char *name, const char *str, int len) {
	char *nname = _smat_get_attribute_name(dsetname, name);
	smat_write_string(s, nname, str, len);
	free(nname);
}

