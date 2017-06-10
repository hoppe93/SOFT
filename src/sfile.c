/* SOFT advanced file management */

#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#ifdef USE_HDF5
#	include "shdf5.h"
#endif
#ifdef USE_MATLAB
#	include "smat.h"
#endif
#include "ssdt.h"
#include "sfile.h"

/*
void (*_sfile_close)(void);
double** (*_sfile_get_doubles)(const char*, sfilesize_t*);
char* (*_sfile_get_string)(const char*);
void (*_sfile_open)(const char*, enum sfile_mode);
void (*_sfile_write_array)(const char*, double**, int, int);
void (*_sfile_write_attribute_scalar)(const char*, const char*, double);
void (*_sfile_write_attribute_string)(const char*, const char*, const char*, int);
void (*_sfile_write_image)(const char*, double**, int);
void (*_sfile_write_list)(const char*, double*, int);
void (*_sfile_write_string)(const char*, const char*, int);

#pragma omp threadprivate(_sfile_close,_sfile_get_doubles,_sfile_get_string,_sfile_open,_sfile_write_array, \
	_sfile_write_attribute_scalar,_sfile_write_attribute_string,_sfile_write_image,_sfile_write_list, \
	_sfile_write_string)

void sfile_close(void){(*_sfile_close)();}
double **sfile_get_doubles(const char *n, sfilesize_t *s){return (*_sfile_get_doubles)(n, s);}
char *sfile_get_string(const char *n){return (*_sfile_get_string)(n);}
void sfile_open(const char *n, enum sfile_mode m){(*_sfile_open)(n, m);}
void sfile_write_array(const char *n, double **d, int i, int j){(*_sfile_write_array)(n, d, i, j);}
void sfile_write_attribute_scalar(const char *n, const char *m, double q){(*_sfile_write_attribute_scalar)(n, m, q);}
void sfile_write_attribute_string(const char *n, const char *m, const char *q, int i){(*_sfile_write_attribute_string)(n, m, q, i);}
void sfile_write_image(const char *n, double **d, int i){(*_sfile_write_image)(n, d, i);}
void sfile_write_list(const char *n, double *d, int i){(*_sfile_write_list)(n, d, i);}
void sfile_write_string(const char *n, const char *m, int i){(*_sfile_write_string)(n, m, i);}
*/

/**
 * Identify the filetype from the filename
 *
 * filename: Name of file
 */
enum sfile_type sfile_get_filetype(const char *filename) {
	/* Find dot (go backwards to allow relative paths) */
	int i = strlen(filename);
	for (;i>=0&&*(filename+i)!='.';i--);
	i++;

	return sfile_name2filetype(filename+i);
}
/**
 * Convert a filetype name to 'enum sfile_type'.
 *
 * filetype: Name of filetype
 */
enum sfile_type sfile_name2filetype(const char *filetype) {
	if (!strcmp(filetype, "h5")) return FILETYPE_HDF5;
	if (!strcmp(filetype, "hdf5")) return FILETYPE_HDF5;
	if (!strcmp(filetype, "mat")) return FILETYPE_MATLAB;
	if (!strcmp(filetype, "sdt")) return FILETYPE_SDT;
	if (!strcmp(filetype, "out")) return FILETYPE_SDT;
	
	return FILETYPE_UNKNOWN;
}

void sfile_deinit(sFILE *s) {
	free(s);
}

sFILE *sfile_init(enum sfile_type ftype) {
	sFILE *s = malloc(sizeof(sFILE));
	s->identifier = NULL;

	switch (ftype) {
		case FILETYPE_HDF5:
#ifdef USE_HDF5
			s->close = shdf5_close;
			s->get_doubles = shdf5_get_doubles;
			s->get_string = shdf5_get_string;
			s->open = shdf5_open;
			s->write_array = shdf5_write_array;
			s->write_attribute_scalar = shdf5_write_attribute_scalar;
			s->write_attribute_string = shdf5_write_attribute_string;
			s->write_image = shdf5_write_image;
			s->write_list = shdf5_write_list;
			s->write_string = shdf5_write_string;
#else
			fprintf(stderr, "ERROR: This version of SOFT was not compiled with HDF5 support.\n");
			exit(-1);
#endif
			break;
		case FILETYPE_MATLAB:
#ifdef USE_MATLAB
			s->close = smat_close;
			s->get_doubles = smat_get_doubles;
			s->get_string = smat_get_string;
			s->open = smat_open;
			s->write_array = smat_write_array;
			s->write_attribute_scalar = smat_write_attribute_scalar;
			s->write_attribute_string = smat_write_attribute_string;
			s->write_image = smat_write_image;
			s->write_list = smat_write_list;
			s->write_string = smat_write_string;
#else
			fprintf(stderr, "ERROR: This version of SOFT was not compiled with MATLAB support.\n");
			exit(-1);
#endif
			break;
		case FILETYPE_SDT:
			s->close = ssdt_close;
			s->get_doubles = ssdt_get_doubles;
			s->get_string = ssdt_get_string;
			s->open = ssdt_open;
			s->write_array = ssdt_write_array;
			s->write_attribute_scalar = ssdt_write_attribute_scalar;
			s->write_attribute_string = ssdt_write_attribute_string;
			s->write_image = ssdt_write_image;
			s->write_list = ssdt_write_list;
			s->write_string = ssdt_write_string;
			break;
		default:
			fprintf(stderr, "ERROR: Trying to open file of unrecognized format.\n");
			exit(-1);
	}

	return s;
}
