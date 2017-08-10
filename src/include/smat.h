#ifndef _SMAT_H
#define _SMAT_H

#include "sfile.h"

void smat_close(sFILE*);
char *smat_get_string(sFILE*, const char*);
double **smat_get_doubles(sFILE*, const char*, sfilesize_t*);
void smat_open(sFILE*, const char*, enum sfile_mode);
void smat_write_array(sFILE*, const char*, double**, size_t, size_t);
void smat_write_attribute_scalar(sFILE*, const char*, const char*, double);
void smat_write_attribute_string(sFILE*, const char*, const char*, const char*, size_t);
void smat_write_image(sFILE*, const char*, double**, size_t);
void smat_write_list(sFILE*, const char*, double*, size_t);
void smat_write_string(sFILE*, const char*, const char*, size_t);

#endif/*_SMAT_H*/
