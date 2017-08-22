/* Particle generator */
#ifndef _PARTICLES_H
#define _PARTICLES_H

typedef struct {
	double t0;
	double tend;
	double mass;
	double charge;
	double vpar, vperp;
	double *v0;
	double *r0;
	double diffel;		/* Differential element */
	int gc_position;	/* If 1, assumes the guiding-center position (rather than particle pos.) to be given */
	int ir, iv1, iv2;	/* Radial index, velocity 1 & 2 indices */
} particle;

enum particles_inputtype {
	PARTICLES_NONE=0,
	PARTICLES_PPAR=1,
	PARTICLES_PPERP=2,
	PARTICLES_PITCH=4,
	PARTICLES_COSPITCH=8,
	PARTICLES_P=16,
	PARTICLES_E=32,
	PARTICLES_MU=128
};

struct particles_paramspec {
	enum particles_inputtype type;			/* Type of this parameter */
	char *name;								/* Parameter name */
	enum particles_inputtype compatible;	/* Parameters which this parameter can be combined with */
	double conversion;						/* Conversion factor */
};

enum particles_generation_type {
	PARTICLES_GT_EVEN,
	PARTICLES_GT_QUEUE
};
enum particles_radial_coordinate {
	PARTICLES_RC_MAJOR,
	PARTICLES_RC_DYNAMIC
};

struct particlespec {
	double t0, tend;
	double mass, charge;
	double r0, r1;
	double rdynmax;
	double ppar0, ppar1;
	double pperp0, pperp1;
	double pitch0, pitch1;
	double cospitch0, cospitch1;
	double p0, p1;
	enum particles_inputtype inputtype;		/* Bits set according to 'enum particles_inputtype' above */
	enum particles_generation_type gentype;
	int rn, pparn, pperpn, pitchn, cospitchn, pn;
	enum particles_radial_coordinate rcoord;/* Type of radial coordinate */
	int gc_position;			/* Guiding-center position (rather than particle pos.) is given */
};

void particles_init(struct particlespec*);
void particles_init_run(int, int, int, int);
int particles_gridsize(void);
void particles_gridsize3(int[3]);
void particles_indices(int[3]);
int particles_get_gridsize(struct particlespec*);
void particles_set_magnetic_axis(double);
int particles_numberof(void);
int particles_stop_condition(void);
void particles_params2p(double, double, double*, double*, double*);
particle *particles_generate(void);
particle *particles_generate_at(double, double, double);
double *particles_get_bounds(void);
double particles_get_drho(void);
double particles_get_dvel1(void);
double particles_get_dvel2(void);
double particles_find_axis_r(double, double);
char *particles_param1_name(void);
char *particles_param2_name(void);
void particles_set_progress(int);
void particles_reportprogress(int, int, int, int);

#endif/*_PARTICLES_H*/
