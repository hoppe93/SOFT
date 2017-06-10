#ifndef _CTSV_H
#define _CTSV_H

#include "IO_data.h"
#include "settings.h"

/* Write data to CSV or TSV file */
void ctsv_write(char*, char, solution_data*, particle*);

/* Test function for this module */
void ctsv_test(void);

#endif/*_CTSV_H*/
