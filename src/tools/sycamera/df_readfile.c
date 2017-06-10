/* Distribution function readfile */

#include <stdio.h>
#include "distfunc.h"

char *df_readfile_file;	/* Buffered file */
char df_readfile_buffer[DISTFUNC_BUFFER_SIZE+1];
char df_readfile_newlineflag;
int df_readfile_linenr=1, df_readfile_eof_flag=0, df_readfile_ptr=0;

char df_readfile_getc(void) { return df_readfile_file[df_readfile_ptr++]; }
int df_readfile_eof(void) {
	return df_readfile_eof_flag;
}
void df_readfile_reset(void) {
	df_readfile_eof_flag = 0;
	df_readfile_linenr=1;
}
void df_readfile_to_eol(void) {
	int c;

	while ((c=df_readfile_getc())!=0 && c != '\n');
	if (c == 0) df_readfile_eof_flag = 1;
	else df_readfile_newlineflag = 1;
	df_readfile_linenr++;
}
char df_readfile_stopped_at_newline(void) {
	return df_readfile_newlineflag;
}
void df_readfile_loaddbl(double *val) {
	char *buf;
	do {
		buf = df_readfile_value();
	} while (!df_readfile_eof() && *buf==0);

	if (*buf == 0) {fprintf(stderr, "ERROR: line %d: Badly formatted file!\n", df_readfile_linenr); exit(EXIT_FAILURE);}
	if (!df_readfile_eof() && val != NULL) *val = atof(buf);
}
void df_readfile_loadint(int *val) {
	char *buf;
	do {
		buf = df_readfile_value();
	} while (!df_readfile_eof() && *buf==0);

	if (*buf == 0 && !df_readfile_eof()) {fprintf(stderr, "ERROR: line %d: Badly formatted file!\n", df_readfile_linenr); exit(EXIT_FAILURE);}
	if (!df_readfile_eof() && val != NULL) *val = atof(buf);
}
/**
 * Reads the next value from file into the buffer.
 * If buffer is empty after a call to word(),
 * an empty line has been discovered.
 *
 * f: Pointer to file to read from
 * RETURNS 0 IF EMPTY LINE WAS ENCOUNTERED
 */
char *df_readfile_value(void) {
	int c, i=0;
	df_readfile_newlineflag = 0;

	/* Get character */
	while (i<DISTFUNC_BUFFER_SIZE && (c=df_readfile_getc())!=0) {
		/* If newline... */
		if (c == '\n') {
			df_readfile_newlineflag = 1;
			df_readfile_linenr++;
			break;		/* ...stop... */
		/* If comment marker */
		} else if (c == '#') {
			df_readfile_to_eol();
			break;
		/* If space, tab or comma... */
		} else if (c == ' ' || c == '\t' || c == ',') {
			/* If buffer is empty... */
			if (i == 0) continue;	/* ...do nothing and continue */
			else break;		/* ...stop... */
		} else	/* Append char to buffer and continue */
			df_readfile_buffer[i++] = (char)c;
	}

	if (c == 0) df_readfile_eof_flag = 1;

	df_readfile_buffer[i] = 0;
	return df_readfile_buffer;
}

/**
 * Skip a given number of lines of file
 *
 * n: Number of lines to skip
 * f: Pointer to file to read from
 */
void df_readfile_skip_lines(int n) {
	int c=0;
	/* While all lines have not been skipped,
	 * and we haven't reached the end-of-file
	 */
	while (n>0 && (c=df_readfile_getc())!=0) {
		/* If c = newline, decrease counter */
		if (c == '\n') n--;
	}

	if (c == 0) df_readfile_eof_flag = 1;
}

void df_readfile_unload(void) {
	free(df_readfile_file);
}
void df_readfile_load(const char *filename) {
	FILE *f;
    f = fopen(filename, "r");
	if (!f) {
		perror("ERROR");
		fprintf(stderr, "ERROR: Unable to open distribution function file: `%s'!\n", filename);
		exit(EXIT_FAILURE);
	}

	fseek(f, 0, SEEK_END);
	size_t filesize = ftell(f);
	fseek(f, 0, SEEK_SET);

	df_readfile_file = malloc(sizeof(char)*(filesize+1));
	if (fread(df_readfile_file, sizeof(char), filesize, f) != filesize) {
		fprintf(stderr, "ERROR: Unable to read distribution function to the end!\n");
		exit(EXIT_FAILURE);
	}
	df_readfile_file[filesize] = 0;
}
