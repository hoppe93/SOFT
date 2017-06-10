/* Data file read helpers */

#include <stdio.h>

#define BUFFER_SIZE 1023
char file_buffer[BUFFER_SIZE+1];

void readfile_to_eol(FILE *f) {
	int c;
	
	while ((c=fgetc(f))!=EOF && c != '\n');
}
/**
 * Reads the next word from file into the buffer.
 * If buffer is empty after a call to word(),
 * an empty line has been discovered.
 *
 * f: Pointer to file to read from
 * RETURNS 1 IF EMPTY LINE WAS ENCOUNTERED
 */
char *readfile_word(FILE *f) {
	int c, i=0;

	/* Get character */
	while (i<BUFFER_SIZE && (c=fgetc(f))!=EOF) {
		/* If newline... */
		if (c == '\n')
			break;		/* ...stop... */
		/* If comment marker */
		else if (c == '#') {
			readfile_to_eol(f);
			break;
		/* If space, tab or equal sign... */
		} else if (c == ' ' || c == '\t' || c == '=') {
			/* If buffer is empty... */
			if (i == 0) continue;	/* ...do nothing and continue */
			else break;		/* ...stop... */
		} else	/* Append char to buffer and continue */
			file_buffer[i++] = (char)c;
	}

	file_buffer[i] = 0;
	return file_buffer;
}
/**
 * Skip a given number of lines of file
 *
 * n: Number of lines to skip
 * f: Pointer to file to read from
 */
void readfile_skip_lines(int n, FILE *f) {
	int c;
	/* While all lines have not been skipped,
	 * and we haven't reached the end-of-file
	 */
	while (n>0 && (c=fgetc(f))!=EOF) {
		/* If c = newline, decrease counter */
		if (c == '\n') n--;
	}
}
