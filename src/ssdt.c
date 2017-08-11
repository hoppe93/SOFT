/**
 * A SOFT SDT (Self-Descriptive Text) interface
 * for simplified I/O.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "sfile.h"
#include "ssdt.h"

//char *ssdt_file;
//int ssdt_keys;
//ssdt_key *ssdt_map;
//#pragma omp threadprivate(ssdt_f,ssdt_file,ssdt_keys,ssdt_map)

/**
 * Open SDT file
 *
 * filename: Name of file
 */
int ssdt_open(sFILE *s, const char *filename, enum sfile_mode mode) {
	FILE *f;
	s->mode = mode;

	switch (mode) {
		case SFILE_MODE_READ:
			f = fopen(filename, "r");
			if (f == NULL) {
				fprintf(stderr, "Unable to open SDT file: %s\n", filename);
				return 0;
			}
			s->identifier = _ssdt_load(f);
			fclose(f);
			return 1;
		case SFILE_MODE_UPDATE:
			s->identifier = fopen(filename, "w+");
			break;
		case SFILE_MODE_WRITE:
			s->identifier = fopen(filename, "w");
			break;
		default:
			fprintf(stderr, "Unrecognized option for opening SDT file: %d.\n", mode);
			s->identifier = NULL;
			return 0;
	}

	if (s->identifier == NULL) {
		fprintf(stderr, "Unable to open SDT file: %s\n", filename);
		return 0;
	}

	return 1;
}

/**
 * Close currently opened SDT file.
 */
void ssdt_close(sFILE *s) {
	int i;
	if (s->mode == SFILE_MODE_READ) {
		struct ssdt_keylist *kl = (struct ssdt_keylist*)s->identifier;
		if (kl == NULL) return;

		for (i = 0; i < kl->nkeys; i++) {
			free(kl->keys[i].name);
		}
		free(kl->keys);
		free(s->identifier);
		s->identifier = NULL;
	} else {
		if (s->identifier != NULL) {
			fclose((FILE*)s->identifier);
			s->identifier = NULL;
		}
	}

}

/******************************
 *********** INPUT ************
 ******************************/

/**
 * @INTERNAL
 * Load the entire file
 */
struct ssdt_keylist *_ssdt_load(FILE *f) {
	char *ssdt_file;
	size_t read;

	fseek(f, 0, SEEK_END);
	size_t l = ftell(f);
	fseek(f, 0, SEEK_SET);
	ssdt_file = malloc(sizeof(char)*(l+1));
	read = fread(ssdt_file, sizeof(char), l, f);

	if (read != l) {
		fprintf(stderr, "ERROR: Unable to load SDT file.\n");
		return NULL;
	}

	/* Create a map of the file */
	struct ssdt_keylist *kl;
	kl = malloc(sizeof(struct ssdt_keylist));
	kl->keys = NULL;
	kl->nkeys = 0;

	int i = 0, j, m, n;
	while (ssdt_file[i]) {
		kl->keys = realloc(kl->keys, sizeof(ssdt_key)*(kl->nkeys+1));

		/* Load m & n */
		sscanf(ssdt_file+i, "%d %d", &m, &n);
		/* Skip m & n in the string */
		while (ssdt_file[i]&& ssdt_file[i]!=' '&&ssdt_file[i]!='\t') i++;
		if (ssdt_file[i]==0) {fprintf(stderr, "ERROR[1]: Unexpected end-of-file at byte %d.\n", i); return NULL;}
		while (ssdt_file[i]&&(ssdt_file[i]==' '||ssdt_file[i]=='\t')) i++;
		if (ssdt_file[i]==0) {fprintf(stderr, "ERROR[2]: Unexpected end-of-file at byte %d.\n", i); return NULL;}
		while (ssdt_file[i]&& ssdt_file[i]!=' '&&ssdt_file[i]!='\t') i++;
		if (ssdt_file[i]==0) {fprintf(stderr, "ERROR[3]: Unexpected end-of-file at byte %d.\n", i); return NULL;}
		while (ssdt_file[i]&&(ssdt_file[i]==' '||ssdt_file[i]=='\t')) i++;
		if (ssdt_file[i]==0) {fprintf(stderr, "ERROR[4]: Unexpected end-of-file at byte %d.\n", i); return NULL;}
		kl->keys[kl->nkeys].m = m;
		kl->keys[kl->nkeys].n = n;

		/* Load name */
		j = i;
		while (ssdt_file[i] && ssdt_file[i]!='\n') i++;
		if (ssdt_file[i]==0) {fprintf(stderr, "ERROR[5]: Unexpected end-of-file at byte %d.\n", i); return NULL;}
		kl->keys[kl->nkeys].name = malloc(sizeof(char)*(i-j+1));
		strncpy(kl->keys[kl->nkeys].name, ssdt_file+j, i-j);
		kl->keys[kl->nkeys].name[i-j] = 0;

		/* Find location of data */
		while (ssdt_file[i] && ssdt_file[i]!='\n') i++;
		if (ssdt_file[i]==0) {fprintf(stderr, "ERROR[6]: Unexpected end-of-file at byte %d.\n", i); return NULL;}
		while (ssdt_file[i] && (ssdt_file[i]=='\n' || ssdt_file[i]==' ' || ssdt_file[i]=='\t')) i++;
		if (ssdt_file[i]==0) {fprintf(stderr, "ERROR[7]: Unexpected end-of-file at byte %d.\n", i); return NULL;}
		kl->keys[kl->nkeys].location = ssdt_file+i;

		kl->nkeys++;

		/* Seek to end of object or file */
		while (ssdt_file[i]&&!(ssdt_file[i]=='\n'&&ssdt_file[i+1]=='\n')) i++;
		if (ssdt_file[i] == 0) break;
		else while (ssdt_file[i]=='\n') i++;
	}

	free(ssdt_file);
	return kl;
}
/**
 * Locate the object with name "name".
 * Places the file pointer at the beginning of
 * the object data and returns object dimensions.
 */
ssdt_key *ssdt_locate(sFILE *s, const char *name) {
	int i;
	struct ssdt_keylist *kl = (struct ssdt_keylist*)s->identifier;

	for (i = 0; i < kl->nkeys; i++) {
		if (!strcmp(kl->keys[i].name, name))
			return kl->keys+i;
	}

	return NULL;
}

/**
 * Reads the string with name "name".
 *
 * name: Name of string to read.
 */
char *ssdt_get_string(sFILE *s, const char *name) {
	char *str;
	int i;
	ssdt_key *key;

	key = ssdt_locate(s, name);
	if (key == NULL) return NULL;

	/* Locate end-of-object */
	i = 0;
	while (key->location[i] && !(key->location[i]=='\n'&&key->location[i+1]=='\n')) i++;
	i--;

	str = malloc(sizeof(char)*(i+1));
	strncpy(str, key->location, i);
	str[i] = 0;

	return str;
}

/**
 * Reads an array of doubles from the object
 * with name "name". The dimensions of the
 * array are returned in "dims".
 *
 * name: Name of object
 * dims: Contains size of object on return
 *
 * Returns a 1-D array (logically 2-D) which
 * contains the data of the object. If the
 * named object does not exist, NULL is returned.
 */
double **ssdt_get_doubles(sFILE *s, const char *name, sfilesize_t *dims) {
	double **values;
	double *v;
	ssdt_key *key;
	int i,j, vindex, tot;

	key = ssdt_locate(s, name);
	if (key == NULL) return NULL;

	dims[0] = key->m;
	dims[1] = key->n;
	tot = key->m*key->n;

	v = malloc(sizeof(double)*tot);
	values = malloc(sizeof(double*)*key->m);
	values[0] = v;
	vindex = 0;

	i = j = 0;
	while (key->location[i] && j<tot && !(key->location[i]=='\n'&&key->location[i+1]=='\n')) {
		sscanf(key->location+i, "%lf", v+j); j++;

		while (key->location[i] && key->location[i]!=' ' && key->location[i]!='\t' && key->location[i]!='\n') i++;
		while (key->location[i] && (key->location[i]==' ' || key->location[i]=='\t')) i++;

		if (key->location[i]=='\n' && vindex+1 < key->m) {
			values[vindex+1] = values[vindex] + key->n;
			vindex++;
			i++;
		}
	}

	return values;
}

/******************************
 *********** OUTPUT ***********
 ******************************/
void ssdt_write_array(sFILE *s, const char *name, double **arr, int m, int n) {
	char *str;
	const int ssize = 2048;
	int i, j, k, l;

	if (s == NULL || s->mode == SFILE_MODE_READ) {
		fprintf(stderr, "ERROR: No SDT file is open, so unable to write array.\n");
		return;
	}

	fprintf((FILE*)s->identifier, "%d %d %s\n", m, n, name);

	str = malloc(sizeof(char)*(ssize+1));
	l = 0;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n-1; j++) {
			k = snprintf(str+l, ssize-l, "%.8e  ", arr[i][j]);
			if (k < ssize-l) l += k;
			else {
				fprintf((FILE*)s->identifier, "%s%.8e  ", str, arr[i][j]);
				l = 0;
			}
		}

		k = snprintf(str+l, ssize-l, "%.8e\n", arr[i][j]);
		if (k < ssize-l) l += k;
		else {
			fprintf((FILE*)s->identifier, "%s%.8e\n", str, arr[i][j]);
			l = 0;
		}
	}
	if (l > 0) fprintf((FILE*)s->identifier, "%s\n", str);
}
void ssdt_write_image(sFILE *s, const char *name, double **img, int pixels) {
	ssdt_write_array(s, name, img, pixels, pixels);
}
void ssdt_write_list(sFILE *s, const char *name, double *list, int n) {
	ssdt_write_array(s, name, &list, 1, n);
}
void ssdt_write_string(sFILE *s, const char *name, const char *str, int len) {
	if (s->identifier == NULL || s->mode == SFILE_MODE_READ) {
		fprintf(stderr, "ERROR: No SDT file is open, so unable to write array.\n");
		return;
	}

	int i;
	for (i = 0; i < len; i++) {
		if (str[i] == 0 || str[i] == '\n') {
			fprintf(stderr, "ERROR: Invalid character in string about to be written to SDT file: 0x%X\n", (unsigned int)str[i]);
			return;
		}
	}
	if (str[i] != 0) {
		fprintf(stderr, "ERROR: String about to be written to SDT file is not terminated by null character.\n");
		return;
	}

	fprintf((FILE*)s->identifier, "1 %d %s\n%s\n\n", len, name, str);
}
char *_ssdt_get_attribute_name(const char *dsetname, const char *name) {
	int l1 = strlen(dsetname), l2 = strlen(name);
	char *nname = malloc(sizeof(char)*(l1+l2+2));

	strcpy(nname, dsetname);
	nname[l1] = '_';
	strcpy(nname+l1+1, name);
	
	return nname;
}
void ssdt_write_attribute_scalar(sFILE *s, const char *dset, const char *name, double val) {
	char *fullname = _ssdt_get_attribute_name(dset, name);
	fprintf((FILE*)s->identifier, "1 1 %s\n%.8e\n\n", fullname, val);
	free(fullname);
}
void ssdt_write_attribute_string(sFILE *s, const char *dset, const char *name, const char *str, int len) {
	char *fullname = _ssdt_get_attribute_name(dset, name);
	ssdt_write_string(s, fullname, str, len);
	free(fullname);
}
