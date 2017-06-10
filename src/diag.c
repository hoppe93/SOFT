/* A SOFT diagnostic utility
 *
 * This utility can be used to monitor various
 * parameters during a SOFT run. Consider for example
 * a case where you only want output when certain
 * parameters take on certain values. In this case
 * outputting data in every step is very inefficient
 * and gives huge amounts of clearly uninteresting
 * data. To avoid this the "diag" utility can be
 * used. The diag utility is used in four "steps":
 *
 *   1. Initialize the utility when the program starts (once per thread)
 *   2. Register a parameter by its name. The utility returns
 *      a parameter index which is used to update the parameters
 *      value later on.
 *   3. Update the parameters you would like to monitor using
 *      the "diag_update()" routine, passing the parameter index
 *      and the parameter value to the function.
 *   4. When the condition for having "interesting" data is
 *      met (must be checked separately), call the function
 *      "diag_trigger()" which will save the current values of
 *      all parameters to a buffer.
 *
 * When the thread is done with all computations, the function
 * "diag_deinit()" should be called with the name of the file to
 * which all data should be written.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "diag.h"

char **diag_names, *diag_buffer;
double *diag_values;
int diag_n, diag_buffer_length;
const char *diag_filename;

#pragma omp threadprivate(diag_names, diag_buffer, diag_buffer_length, diag_values, diag_n, diag_filename)

void diag_init(const char *ofname) {
	diag_names = NULL;
	diag_values = NULL;
	diag_buffer = NULL;
	diag_buffer_length = 0;
	diag_n = 0;

	diag_filename = ofname;

	/* Let first thread delete any old 'diag' files */
	if (omp_get_thread_num() == 0) {
		remove(ofname);
	}
}

/**
 * Register a parameter to monitor. The name
 * is passed to this function and will appear
 * in the output file.
 *
 * The function returns a "parameter index" which
 * is used to refer to a parameter when communicating
 * with this utility.
 *
 * name: Name of parameter
 */
int diag_register(const char *name) {
	diag_names = realloc(diag_names, sizeof(char*)*(diag_n+1));
	diag_values = realloc(diag_values, sizeof(double)*(diag_n+1));
	diag_values[diag_n] = 0;

	long l = strlen(name)+1;
	diag_names[diag_n] = malloc(sizeof(char*)*l);
	strcpy(diag_names[diag_n], name);

	diag_n++;
	return diag_n-1;
}
double diag_get(int diag_index) {
	return diag_values[diag_index];
}
void diag_update(int diag_index, double value) {
	diag_values[diag_index] = value;
}
void diag_increase(int diag_index) {
	diag_values[diag_index]++;
}
void diag_trigger(void) {
	/* Output newline first */
	int l=0, i;
	char buffer[diag_n*20 + diag_n];
	buffer[l++] = '\n';

	for (i = 0; i < diag_n; i++) {
		l += sprintf(buffer+l, "%.12e\t", diag_values[i]);
	}

	diag_buffer = realloc(diag_buffer, sizeof(char)*(diag_buffer_length+l+1));
	strcpy(diag_buffer+diag_buffer_length, buffer);
	diag_buffer_length += l;
}

void diag_deinit(void) {
	int i;
	FILE *f;

	if (diag_buffer == NULL) return;

#pragma omp critical
{
	f = fopen(diag_filename, "a");

	if (!f) {
		perror("ERROR");
		fprintf(stderr, "Unable to open file '%s' for writing diagnostic information!\n", diag_filename);
		exit(-1);
	}

	/* If this thread is the first to access this file,
	 * also write the header! */
	if (ftell(f) == 0) {
		fprintf(f, "#");
		for (i = 0; i < diag_n; i++) {
			fprintf(f, " %s", diag_names[i]);
		}
	}

	fputs(diag_buffer, f);
	fclose(f);
}
}

