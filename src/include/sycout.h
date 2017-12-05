/* Synchrotron Camera output handler */
#ifndef _SYCOUT_H
#define _SYCOUT_H

#include <stdio.h>
#include "settings.h"
#include "sfile.h"
#include "tools.h"

#define SPACE3D_CHUNKSIZE 10000

enum sycout_mpiid {
	SYCOUT_MPIID_GENERAL,
	SYCOUT_MPIID_GREEN,
	SYCOUT_MPIID_IMAGE,
	SYCOUT_MPIID_POLIMAGE,
	SYCOUT_MPIID_POLSPECTROMETER,
	SYCOUT_MPIID_SPACE3D,
	SYCOUT_MPIID_SPECTROMETER,
	SYCOUT_MPIID_TOPVIEW
};

struct sycout_data {
    step_data *sd;      /* General information about step */
    double brightness;  /* Brightness in point */
	double differential;/* Differential element */
	double RdPhi;		/* Toroidal differential element */
	double Jdtdrho;		/* Poloidal trajectory differential */
	double Jp;			/* Momentum-space jacobian */
	double particle_diffel;/* Particle differential element */
	double distribution_function;/* Distribution function evaluated in current point */
    double i, j;           /* Pixel indices */
};

typedef struct {
    char *name;
    void (*deinit_run)(void);
    void (*init)(struct general_settings*);
    void (*init_run)(void);
    void (*init_particle)(particle*);
    void (*step)(struct sycout_data*);
    void (*write)(int, int);
} sycout_type;

void sycout_deinit_run(void);
void sycout_init_handler(void);
void sycout_output_all(int, int);
void sycout_prepare_run(void);
void sycout_init_particle(particle*);
void sycout_select(const char*, struct general_settings*, int);
void sycout_step(struct sycout_data*);

/****************************
 *       GREEN SYCOUT       *
 ****************************/
enum sycout_green_dimension {
	SYCOUT_GREEN_NONE,
	SYCOUT_GREEN_RADIUS,
	SYCOUT_GREEN_IMAGEI,
	SYCOUT_GREEN_IMAGEJ,
	SYCOUT_GREEN_SPECTRUM,
	SYCOUT_GREEN_VEL1,
	SYCOUT_GREEN_VEL2
};

void sycout_green_init(struct general_settings*);
void sycout_green_init_run(void);
void sycout_green_init_particle(particle*);
void sycout_green_deinit_run(void);
void sycout_green_step(struct sycout_data*);
void sycout_green_write(int, int);

/****************************
 *       IMAGE SYCOUT       *
 ****************************/
typedef struct {
	double **canvas;	/* Camera canvas */
	int pixels;			/* Pixels per dimension */
} camera_image;

enum sycout_image_type {
	SYCOUT_IMAGE_TYPE_BW,			/* Black and white image (0's and 1's) */
	SYCOUT_IMAGE_TYPE_INTENSITY,	/* Intensity distribution */
	SYCOUT_IMAGE_TYPE_HIST			/* Histogram (pixels hold integers corresponding to number of hits) */
};

void sycout_image_deinit_run(void);
void sycout_image_init(struct general_settings*);
void sycout_image_init_run(void);
void sycout_image_init_particle(particle*);
void sycout_image_step(struct sycout_data*);
void sycout_image_combine(sFILE*, camera_image*);
void sycout_image_output(sFILE*, camera_image*);
void sycout_image_write(int, int);

/****************************
 *  POLARIZED IMAGE SYCOUT  *
 ****************************/
typedef struct {
	double **StokesI, **StokesQ,
		   **StokesU, **StokesV;
	int pixels;
} camera_polimage;

void sycout_polimage_deinit_run(void);
void sycout_polimage_init(struct general_settings*);
void sycout_polimage_init_run(void);
void sycout_polimage_init_particle(particle*);
void sycout_polimage_step(struct sycout_data*);
void sycout_polimage_combine(sFILE*, camera_polimage*);
void sycout_polimage_output(sFILE*, camera_polimage*);
void sycout_polimage_write(int, int);

/****************************
 *      TOPVIEW SYCOUT      *
 ****************************/
enum sycout_topview_type {
	SYCOUT_TOPVIEW_TYPE_BW,			/* Black and white image (0's and 1's) */
	SYCOUT_TOPVIEW_TYPE_INTENSITY,	/* Intensity distribution */
	SYCOUT_TOPVIEW_TYPE_HIST		/* Histogram (pixels hold integers corresponding to number of hits) */
};

void sycout_topview_deinit_run(void);
void sycout_topview_init(struct general_settings*);
void sycout_topview_init_run(void);
void sycout_topview_init_particle(particle*);
void sycout_topview_step(struct sycout_data*);
void sycout_topview_write(int, int);

/****************************
 *      SPACE3D SYCOUT      *
 ****************************/
typedef struct {
	int n;
	double *x, *y, *z;	/* Number of (recorded) particles */
	double *intensity;	/* Intensity in each point */
} space3d_real_t;

typedef struct {
	size_t pixels;		/* Pixels per dimension */
	double *image;		/* "Image" */
	double xmin, xmax,	/* Bounds in all dimensions */
		   ymin, ymax,
		   zmin, zmax;
} space3d_pixels_t;

#define SOFT_SPACE3D_MAGIC 0x50F8
typedef struct {
	unsigned short magic;
	time_t timestamp;
	double aperture, visang,		/* Aperture and vision angle */
		   toroidal_resolution,
		   dx, dy, dz,				/* Detector viewing direction */
		   px, py, pz;				/* Detector position */
} __attribute__((packed)) space3d_header_t;

enum space3d_model_type {
	SPACE3D_MT_PIXELS,
	SPACE3D_MT_REAL
};

void sycout_space3d_deinit_run(void);
void sycout_space3d_init(struct general_settings*);
void sycout_space3d_init_run(void);
void sycout_space3d_init_particle(particle*);
void sycout_space3d_step(struct sycout_data*);
void sycout_space3d_write(int, int);

/****************************
 *   SPECTROMETER SYCOUT    *
 ****************************/
void sycout_spectrometer_deinit_run(void);
void sycout_spectrometer_init(struct general_settings*);
void sycout_spectrometer_init_run(void);
void sycout_spectrometer_init_particle(particle*);
void sycout_spectrometer_step(struct sycout_data*);
void sycout_spectrometer_write(int, int);

/*********************************
 * POLARIZED SPECTROMETER SYCOUT *
 *********************************/
void sycout_polspectrometer_deinit_run(void);
void sycout_polspectrometer_init(struct general_settings*);
void sycout_polspectrometer_init_run(void);
void sycout_polspectrometer_init_particle(particle*);
void sycout_polspectrometer_step(struct sycout_data*);
void sycout_polspectrometer_write(int, int);

#endif/*_SYCOUT_H*/
