#ifndef _DIAG_H
#define _DIAG_H

/* SOFT internal diagnostics */

void diag_init(const char*);
void diag_deinit(void);
int diag_register(const char*);
double diag_get(int);
void diag_update(int, double);
void diag_increase(int);
void diag_trigger(void);

#endif/*_DIAG_H*/
