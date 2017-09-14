#ifndef _SYCAMERA_H
#define _SYCAMERA_H

#include <time.h>
#include "constants.h"
#include "global.h"
#include "IO_data.h"
#include "rkf45.h"
#include "settings.h"
#include "sycout.h"
#include "tools.h"

#define SOFT_MAGIC 0x50F7

#define SYCAMERA_PCYL_CONSTANT (1/sqrt(3)*LIGHTSPEED*CHARGE*CHARGE/EPS0)
#define SYCAMERA_PCYL_LOWERBOUND (4*PI/3*LIGHTSPEED/CHARGE)

#define SYCAMERA_PAS2_CONSTANT (sqrt(3)*LIGHTSPEED*CHARGE*CHARGE/(8*PI*EPS0))

#define SYCAMERA_CHUNKSIZE 1024	/* Number of elements to allocate to orbit variable on each allocation */

typedef struct {
	unsigned short magic;			/* Magic bytes */
	/* General information */
	time_t timestamp;
	/* Detector information */
	double aperture, visang,	/* Aperture and vision angle */
		   toroidal_resolution,
		   dx, dy, dz,			/* Camera direction */
		   px, py, pz;			/* Camera position */
	/* Particle information */
	double mass, charge;
	double bounds[9];			/* Bounds are given in order:
								   r0, rend, rn
								   ppar0, pparend, pparn
								   pperp0, pperpend, pperpn
								*/
}__attribute__((packed)) mapfile_header;

#define SYCAMERA_MAP_ASCII 1
#define SYCAMERA_MAP_HDF5 2

#define ASCII_HEADER_NAME_MAX 10	/* Maximum length (bytes) of header variable names */
#define ASCII_HEADER_VALUE_MAX 21	/* Maximum length (bytes) of header variable values */
#define ASCII_PRECISION 12			/* Number of digits to keep in floating-point numbers */

enum sycamera_radiation_type {
	SYCAMERA_RADIATION_BREMSSTRAHLUNG,
	SYCAMERA_RADIATION_BREMSSTRAHLUNG_SPECTRUM,
	SYCAMERA_RADIATION_CONSTANT,
	SYCAMERA_RADIATION_SYNCHROTRON,
	SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM
};

#define MAP_PARTICLES_REALLOC_SIZE 1000
#define MAP_PIXELS_REALLOC_SIZE 10

enum sycamera_polarization_type {
	SYCAMERA_POLARIZATION_BOTH,
	SYCAMERA_POLARIZATION_PARALLEL,
	SYCAMERA_POLARIZATION_PERPENDICULAR
};

/* Shared variables */
extern vector *ddet, *e1, *e2, *Rdet;
extern double rdet, visang, sycamera_zeff;
extern enum sycamera_radiation_type radiation_type;

/* Lookup tables (pcyl) */
extern const int sycamera_pcyl_lookup_count;
extern const double sycamera_pcyl_lookup_int[];
extern const double sycamera_pcyl_lookup_lambda[];

/* Lookup tables (pdist) */
extern const int sycamera_pdist_lookup_count;
extern const double sycamera_pdist_lookup_int1[];
extern const double sycamera_pdist_lookup_int2[];
extern const double sycamera_pdist_lookup_omega[];

/* Functions required by SOFT */
void sycamera_init(struct general_settings*, struct general_settings*, int);
void sycamera_init_run(unsigned int);
ode_solution *sycamera_init_particle(particle*);
void sycamera_deinit_run(void);
void sycamera_step(ode_solution*, step_data*);
int sycamera_process_step(step_data*, double);
int sycamera_stop_condition(void);
void sycamera_output(equation*);
//void sycamera_extend_orbit(void);

/* Functions specific to the detector */
int cone_delta_can_radiation_hit(step_data*, vector*);
void cone_delta_init(enum sycamera_radiation_type, enum sycamera_polarization_type, double*, int, int);
void cone_delta_init_run(void);
void cone_delta_init_particle(particle*);
void cone_delta_init_step(step_data*);
double cone_delta_get_intensity(step_data*,vector*,vector*,vector*,vector*,vector*);
void cone_delta_get_intersections(double, double, double, double, double, double*, double*, int*, int*, double, double);
int cone_delta_process_step(step_data*, double, double);
double cone_delta_radiation_hits(step_data*,vector*,vector*,vector*,vector*,vector*);
void cone_delta_register_radiation(double, double, step_data*, double, double, double);
double *cone_delta_get_wavelengths(void);
double *cone_delta_get_spectrum(void);
int cone_delta_get_spectrum_length(void);

void cone_dist_init(enum sycamera_radiation_type, enum sycamera_polarization_type, double*, int, int);
void cone_dist_init_run(void);
void cone_dist_init_particle(particle*);
void cone_dist_init_step(step_data*);
double cone_dist_get_intensity(step_data*,vector*,vector*,vector*,vector*,vector*);
double cone_dist_Ihat(double, double, double);
double cone_dist_Ihat_spec(double, double, double);
double *cone_dist_get_wavelengths(void);
double *cone_dist_get_spectrum(void);
int cone_dist_get_spectrum_length(void);

void isotropic_init(enum sycamera_radiation_type, enum sycamera_polarization_type, double*, int, int);
void isotropic_init_run(void);
void isotropic_init_particle(particle*);
void isotropic_init_step(step_data*);
double isotropic_intensity(step_data*,vector*,vector*,vector*,vector*,vector*);
double *isotropic_get_wavelengths(void);
double *isotropic_get_spectrum(void);
int isotropic_get_spectrum_length(void);

void sphere_init(enum sycamera_radiation_type, enum sycamera_polarization_type, double*, int, int);
void sphere_init_run(void);
void sphere_init_particle(particle*);
void sphere_init_step(step_data*);
double sphere_intensity(step_data*,vector*,vector*,vector*,vector*,vector*);
double *sphere_get_wavelengths(void);
double *sphere_get_spectrum(void);
int sphere_get_spectrum_length(void);

/* Functions for spectrum weighting */
void sycamera_spectrum_init(double, double, int);
void sycamera_spectrum_init_run(void);
double sycamera_spectrum_weight(step_data*, double, double);
double *sycamera_spectrum_get_wavelengths(void);
double *sycamera_spectrum_get(void);
int sycamera_spectrum_length(void);

void sycamera_bsspec_init(double, double, int);
void sycamera_bsspec_init_run(void);
double sycamera_bsspec_int(double, double, double, double);
double *sycamera_bsspec_get_wavelengths(void);
double *sycamera_bsspec_get_spectrum(void);
int sycamera_bsspec_get_spectrum_length(void);
void sycamera_bsspec_test(void);

void sycamera_pcyl_init(double, double, int);
void sycamera_pcyl_init_run(void);
double sycamera_pcyl_int(double, double, double, double, double);
double *sycamera_pcyl_get_wavelengths(void);
double *sycamera_pcyl_get_spectrum(void);
int sycamera_pcyl_get_spectrum_length(void);

void sycamera_pas2_init(double, double, int);
int sycamera_pas2_valid(double, double, double, double);
double sycamera_pas2_int(void);

void sycamera_pdist_init(double, double, enum sycamera_polarization_type);
void sycamera_pdist_init_run(void);
void sycamera_pdist_init_particle(double);
double sycamera_pdist_int(double, double, double, double, double, double, double, double, double, double, double);
void sycamera_pdist_test(void);

double *sycamera_get_spectrum(void);
double *sycamera_get_wavelengths(void);
int sycamera_get_spectrum_length(void);

/* Output functions */
//void sycamera_image_output(const char*, camera_image*, int);
//void sycamera_asciiout(char*, camera_map*, int, int, int, int);
/* The binout format has no MPI support and so is currently disabled */
//void sycamera_binout(char*, camera_map*, int, int, int, double, double, double, vector*, vector*, double, double);
//void sycamera_hdf5out(char*, camera_map*, int, int, int, int);

/* Debugging functions */
void print_intersections(int, double*, double*, double, double, double, double, double, double, long long int, double, step_data*, double, double, double, double);
void print_particle(step_data*);
void print_particle_struct(step_data*);
void print_particle_struct_p(step_data*);

#endif/*_SYCAMERA_H*/
