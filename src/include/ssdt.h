#ifndef _SSDT_H
#define _SSDT_H

#include <stdio.h>
#include "sfile.h"

typedef struct {
	char *name;			/* Name of object */
	size_t m, n;			/* Object size */
	char *location;		/* Pointer to start of data */
} ssdt_key;
struct ssdt_keylist {
	size_t nkeys;
	ssdt_key *keys;
	char *filedata;		/* Contents of SDT file */
};

int ssdt_open(sFILE *, const char*, enum sfile_mode);
void ssdt_close(sFILE*);
struct ssdt_keylist *_ssdt_load(FILE*);
ssdt_key *ssdt_locate(sFILE*, const char*);
char *ssdt_get_string(sFILE*, const char*);
double **ssdt_get_doubles(sFILE*, const char*, sfilesize_t*);
void ssdt_write_array(sFILE*, const char*, double**, size_t, size_t);
void ssdt_write_attribute_scalar(sFILE*, const char*, const char*, double);
void ssdt_write_attribute_string(sFILE*, const char*, const char*, const char*, size_t);
void ssdt_write_image(sFILE*, const char*, double**, size_t);
void ssdt_write_list(sFILE*, const char*, double*, size_t);
void ssdt_write_string(sFILE*, const char*, const char*, size_t);

#endif/*_SSDT_H*/
