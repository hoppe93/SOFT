#ifndef _QUANTITIES_H
#define _QUANTITIES_H

typedef struct {
	double *values;
	char *name;
} quantity;

quantity *quantities_get(void);
void quantities_init(void);
int quantities_get_no(void);
int quantities_define(char*);
void quantities_report(int, double, unsigned int);
void quantities_expand(int);

#endif/*_QUANTITIES_H*/
