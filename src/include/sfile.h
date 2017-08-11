#ifndef _SOFT_FILETYPE
#define _SOFT_FILETYPE

enum sfile_type {
	FILETYPE_UNKNOWN,
	FILETYPE_HDF5,
	FILETYPE_MATLAB,
	FILETYPE_SDT
};
enum sfile_mode {
	SFILE_MODE_READ,
	SFILE_MODE_UPDATE,
	SFILE_MODE_WRITE
};

typedef long long unsigned int sfilesize_t;

typedef struct s_sfile {
	void *identifier;	/* Points to a value that can be anything, depending on interface */
	enum sfile_mode mode;/* Opening mode */

	/* Member functions */
	void (*close)(struct s_sfile*);
	double **(*get_doubles)(struct s_sfile*, const char*, sfilesize_t*);
	char *(*get_string)(struct s_sfile*, const char*);
	int (*open)(struct s_sfile*, const char*, enum sfile_mode);
	void (*write_array)(struct s_sfile*, const char*, double**, int, int);
	void (*write_attribute_scalar)(struct s_sfile*, const char*, const char*, double);
	void (*write_attribute_string)(struct s_sfile*, const char*, const char*, const char*, int);
	void (*write_image)(struct s_sfile*, const char*, double**, int);
	void (*write_list)(struct s_sfile*, const char*, double*, int);
	void (*write_string)(struct s_sfile*, const char*, const char*, int);
} sFILE;

/*
void sfile_close(void);
double** sfile_get_doubles(const char*, sfilesize_t*);
char* sfile_get_string(const char*);
int sfile_open(const char*, enum sfile_mode);
void sfile_write_array(const char*, double**, int, int);
void sfile_write_attribute_scalar(const char*, const char*, double);
void sfile_write_attribute_string(const char*, const char*, const char*, int);
void sfile_write_image(const char*, double**, int);
void sfile_write_list(const char*, double*, int);
void sfile_write_string(const char*, const char*, int);
*/

enum sfile_type sfile_get_filetype(const char*);
enum sfile_type sfile_name2filetype(const char*);
sFILE *sfile_init(enum sfile_type);

#endif/*_SOFT_FILETYPE*/
