#ifndef _SHDF5_H
#define _SHDF5_H

#include <hdf5.h>
#include "sfile.h"

void shdf5_close(sFILE*);
int shdf5_open(sFILE*, const char*, enum sfile_mode);
char *shdf5_get_string(sFILE*, const char*);
double **shdf5_get_doubles(sFILE*, const char*, sfilesize_t*);
void shdf5_write_string(sFILE*, const char*, const char*, int);
void shdf5_write_array(sFILE*, const char*, double**, int, int);
void shdf5_write_image(sFILE*, const char*, double**, int);
void shdf5_write_list(sFILE*, const char*, double*, int);
void shdf5_write_attribute_scalar(sFILE*, const char*, const char*, double);
void shdf5_write_attribute_string(sFILE*, const char*, const char*, const char*, int);

#endif/*_SHDF5_H*/
