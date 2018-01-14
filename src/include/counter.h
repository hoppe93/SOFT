#ifndef _COUNTER_H
#define _COUNTER_H

#include <stdlib.h>

struct counter {
	char *name;
	size_t n;
};

void counter_init(void);
int counter_create(const char*);
void counter_inc(int);
void counter_merge_all(void);
void counter_print(void);

#endif/*_COUNTER_H*/
