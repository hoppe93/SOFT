#ifndef _READFILE_H
#define _READFILE_H

#include <stdio.h>

/* Read one word from file (returns empty string if empty line) */
char *readfile_word(FILE*);
/* Skip the given number of lines */
void readfile_skip_lines(int, FILE*);

#endif/*_READFILE_H*/
