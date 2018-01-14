/* A SOFT counter utility
 *
 * This utility can be used to gather statistics
 * about how many times various lines are passed
 * in the code.
 *
 *   1. Initialize the utility with a call to
 *      'counter_init()' (once per thread).
 *   2. Register a counter with 'counter_create()'.
 *      This is done per thread. A counter index
 *      is returned which must be provided when
 *      working with the counter.
 *   3. To increase the counter, call 'counter_inc'
 *      with the index of the counter to increase.
 *   4. Before the thread exits, call 'counter_merge'
 *      to merge all counters with global counters of
 *      the same name.
 *
 * At the end of the program, 'counter_print' can
 * be used to dump all counters to screen.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "counter.h"

int counter_ncounters, counter_nglobcounters=0;
struct counter *counter_counters,
			   *counter_global_counters=NULL;

#pragma omp threadprivate(counter_counters,counter_ncounters)

void counter_init(void) {
	counter_counters = NULL;
	counter_ncounters = 0;
}

/**
 * Create a counter with name 'name'.
 *
 * name: Name of counter.
 *
 * RETURNS the index of counter, which must
 * be passed when interacting with a particular
 * counter.
 */
int counter_create(const char *name) {
	/* Add the counter locally */
	counter_counters = realloc(counter_counters, sizeof(struct counter)*(counter_ncounters+1));
	counter_counters[counter_ncounters].name = malloc(strlen(name)+1);
	strcpy(counter_counters[counter_ncounters].name, name);
	counter_counters[counter_ncounters].n = 0;
	counter_ncounters++;

	#pragma omp critical
	{
		/* Add global counter if it doesn't already exists */
		int i, exists=0;
		for (i = 0; i < counter_nglobcounters; i++) {
			if (!strcmp(counter_global_counters[i].name, name)) {
				exists = 1;
				break;
			}
		}

		if (!exists) {
			counter_global_counters = realloc(counter_global_counters, sizeof(struct counter)*(counter_nglobcounters+1));
			counter_global_counters[counter_nglobcounters].name = malloc(strlen(name)+1);
			strcpy(counter_global_counters[counter_nglobcounters].name, name);
			counter_global_counters[counter_nglobcounters].n = 0;
			counter_nglobcounters++;
		}
	}

	return counter_ncounters-1;
}

/**
 * Increment counter
 */
void counter_inc(int counter_index) {
	counter_counters[counter_index].n++;
}

/**
 * Merge all local counters to main thread
 */
void counter_merge_all(void) {
	#pragma omp critical
	{
		int i;
		for (i = 0; i < counter_nglobcounters; i++)
			counter_global_counters[i].n += counter_counters[i].n;
	}
}

/**
 * Print all counters
 */
void counter_print(void) {
	int i;
	for (i = 0; i < counter_nglobcounters; i++) {
		printf("%20s = %zu\n", counter_global_counters[i].name, counter_global_counters[i].n);
	}
}

