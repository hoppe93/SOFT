#ifndef _SSDT_H
#define _SSDT_H

#include <stdio.h>
#include "sfile.h"

typedef struct {
	char *name;			/* Name of object */
	int m, n;			/* Object size */
	char *location;		/* Pointer to start of data */
} ssdt_key;
struct ssdt_keylist {
	int nkeys;
	ssdt_key *keys;
};

int ssdt_open(sFILE *, const char*, enum sfile_mode);
void ssdt_close(sFILE*);
struct ssdt_keylist *_ssdt_load(FILE*);
ssdt_key *ssdt_locate(sFILE*, const char*);
char *ssdt_get_string(sFILE*, const char*);
double **ssdt_get_doubles(sFILE*, const char*, sfilesize_t*);
void ssdt_write_array(sFILE*, const char*, double**, int, int);
void ssdt_write_attribute_scalar(sFILE*, const char*, const char*, double);
void ssdt_write_attribute_string(sFILE*, const char*, const char*, const char*, int);
void ssdt_write_image(sFILE*, const char*, double**, int);
void ssdt_write_list(sFILE*, const char*, double*, int);
void ssdt_write_string(sFILE*, const char*, const char*, int);

#endif/*_SSDT_H*/
