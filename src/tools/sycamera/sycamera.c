/* Synchrotron camera module */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "config.h"
#include "constants.h"
//#include "diag.h"
#include "distfunc.h"
#include "domain.h"
#include "equations.h"
#include "global.h"
#include "IO_data.h"
#include "particles.h"
#include "rkf45.h"
#include "settings.h"
#include "smpi.h"
#include "sycamera.h"
#include "sycout.h"
#include "tools.h"
#include "util.h"

#ifdef USE_MPI
#	include <mpi.h>
#endif

vector  *sol1=NULL, *sol2=NULL,	/* Temporary solution vectors */
		*Rdet, *ddet;			/* Synchrotron detector position and direction */
double  lasttime=0;				/* Time of last snapshot */
double  sycamera_charge=0,		/* Particle charge */
		sycamera_mass=0,		/* Particle mass */
		rdet=0,					/* Synchrotron detector aperture */
		visang=0,				/* Camera vision angle */
		tanvisang=0,			/* Tangent of half vision angle (tan(phi/2)) */
		sycamera_distfunc_weight=0,/* Distribution function weight */
		sycamera_dPhi=1;		/* Set automatically from 'toroidal_resolution' */
int toroidal_resolution=320;	/* Resolution when doing to toroidal rotations */
double *sycamera_costor, *sycamera_sintor;	/* Arrays containing cos/sin of toroidal angles */
double sycamera_lasti, sycamera_lastj,
	   sycamera_lastlx, sycamera_lastly;
double sycamera_particle_diffel;/* Particle differential element */

vector *temps=NULL, *e1, *e2;
const int NUMBER_OF_TEMPS=7;
int sycamera_has_distfunc;		/* Whether or not a distribution function is available */
enum sycamera_radiation_type radiation_type=SYCAMERA_RADIATION_SYNCHROTRON;
enum sycamera_polarization_type sycamera_polarization=SYCAMERA_POLARIZATION_BOTH;

void (*intensity_init)(enum sycamera_radiation_type, enum sycamera_polarization_type, double*, int, int)=NULL;
void (*intensity_init_run)(void)=NULL;
void (*intensity_init_particle)(particle*)=NULL;
void (*intensity_init_step)(step_data*)=NULL;
double (*intensity_function)(step_data*, vector*, vector*, vector*, vector*, vector*)=NULL;
double *(*intensity_spectrum)(void)=NULL;

#pragma omp threadprivate(sol1,sol2,temps,lasttime,sycamera_distfunc_weight,sycamera_lasti,sycamera_lastj,sycamera_lastlx,sycamera_lastly,sycamera_particle_diffel)

/**
 * Initialize the synchrotron camera.
 *
 * set: List of settings for the synchrotron camera
 * sycout_set: List of lists of settings for the sycouts to use
 * nsycout_settings: Number of lists of settings for sycouts
 *   (i.e. length of 'sycout_set')
 */
void sycamera_init(struct general_settings *set, struct general_settings *sycout_set, int nsycout_settings) {
	if (set == NULL) {
		fprintf(stderr, "ERROR: The sycamera tool requires settings to be given!\n");
		exit(-1);
	}

	/* Initialize output modules */
	sycout_init_handler();

	/* Define detector */
	Rdet = vinit(3, .0,.0,.0);
	ddet = vinit(3, .0,.0,.0);
	rdet = 0;
	visang = 0;
	sycamera_has_distfunc = 0;

	sycamera_polarization = SYCAMERA_POLARIZATION_BOTH;

	double *lambdas=NULL;
	int spectrum_resolution=50, integral_resolution=10;

	intensity_init = &cone_delta_init;
	intensity_init_particle = &cone_delta_init_particle;
	intensity_init_step = &cone_delta_init_step;
	intensity_function = &cone_delta_get_intensity;
	intensity_spectrum = &cone_delta_get_spectrum;

	/* Loop over all settings */
	int i;
	for (i = 0; i < set->n; i++) {
		if (!strcmp(set->setting[i], "aperture"))
			rdet = atof(set->value[i]);
		else if (!strcmp(set->setting[i], "cone")) {
			if (!strcmp(set->value[i], "delta")) {
				printf("Selecting cone model 'delta'...\n");
				intensity_init = &cone_delta_init;
				intensity_init_run = &cone_delta_init_run;
				intensity_init_particle = &cone_delta_init_particle;
				intensity_init_step = &cone_delta_init_step;
				intensity_function = &cone_delta_get_intensity;
				intensity_spectrum = &cone_delta_get_spectrum;
			} else if (!strcmp(set->value[i], "dist")) {
				printf("Selecting cone model 'dist'...\n");
				intensity_init = &cone_dist_init;
				intensity_init_run = &cone_dist_init_run;
				intensity_init_particle = &cone_dist_init_particle;
				intensity_init_step = &cone_dist_init_step;
				intensity_function = &cone_dist_get_intensity;
				intensity_spectrum = &cone_dist_get_spectrum;
			} else if (!strcmp(set->value[i], "isotropic")) {
				printf("Selecting cone model 'isotropic'...\n");
				intensity_init = &isotropic_init;
				intensity_init_run = &isotropic_init_run;
				intensity_init_particle = &isotropic_init_particle;
				intensity_init_step = &isotropic_init_step;
				intensity_function = &isotropic_intensity;
				intensity_spectrum = &isotropic_get_spectrum;
			} else if (!strcmp(set->value[i], "sphere")) {
				printf("Selecting cone model 'sphere'...\n");
				intensity_init = &sphere_init;
				intensity_init_run = &sphere_init_run;
				intensity_init_particle = &sphere_init_particle;
				intensity_init_step = &sphere_init_step;
				intensity_function = &sphere_intensity;
				intensity_spectrum = &sphere_get_spectrum;
			} else {
				fprintf(stderr, "Unrecognized cone type given in pi-file: %s\n", set->value[i]);
				exit(-1);
			}
		} else if (!strcmp(set->setting[i], "direction")) {
			ddet->val = atodpn(set->value[i], 3, ddet->val);
			ddet->n = 3;
		} else if (!strcmp(set->setting[i], "distribution_function")) {
			distfunc_load(set->value[i]);
			sycamera_has_distfunc = 1;
		} else if (!strcmp(set->setting[i], "integral_resolution")) {
			integral_resolution = atoi(set->value[i]);
		} else if (!strcmp(set->setting[i], "position")) {
			Rdet->val = atodpn(set->value[i], 3, Rdet->val);
			Rdet->n = 3;
		} else if (!strcmp(set->setting[i], "polarization")) {
			if (!strcmp(set->value[i], "parallel")) {
				sycamera_polarization = SYCAMERA_POLARIZATION_PARALLEL;
			} else if (!strcmp(set->value[i], "perpendicular")) {
				sycamera_polarization = SYCAMERA_POLARIZATION_PERPENDICULAR;
			} else if (!strcmp(set->value[i], "both")) {
				sycamera_polarization = SYCAMERA_POLARIZATION_BOTH;
			} else {
				fprintf(stderr, "sycamera: Unrecognized option for option 'polarization': '%s'.\n", set->value[i]);
				exit(-1);
			}
		} else if (!strcmp(set->setting[i], "product")) {
			sycout_select(set->value[i], sycout_set, nsycout_settings);
		} else if (!strcmp(set->setting[i], "radiation")) {
			if (!strcmp(set->value[i], "bremsstrahlung"))
				radiation_type = SYCAMERA_RADIATION_BREMSSTRAHLUNG;
			else if (!strcmp(set->value[i], "constant"))
				radiation_type = SYCAMERA_RADIATION_CONSTANT;
			else if (!strcmp(set->value[i], "synchrotron"))
				radiation_type = SYCAMERA_RADIATION_SYNCHROTRON;
			else if (!strcmp(set->value[i], "synchrotron_spectrum"))
				radiation_type = SYCAMERA_RADIATION_SYNCHROTRON_SPECTRUM;
			else {
				fprintf(stderr, "Invalid radiation type: %s!\n", set->value[i]);
				exit(-1);
			}
		} else if (!strcmp(set->setting[i], "spectrum")) {
			lambdas = atodpn(set->value[i], 2, NULL);
		} else if (!strcmp(set->setting[i], "spectrum_resolution")) {
			spectrum_resolution = atoi(set->value[i]);
		} else if (!strcmp(set->setting[i], "toroidal_resolution")) {
			toroidal_resolution = atoi(set->value[i]);
		} else if (!strcmp(set->setting[i], "vision_angle"))
			visang = atof(set->value[i]);
		else {
			fprintf(stderr, "Invalid sycamera setting: %s!\n", set->setting[i]);
			exit(-1);
		}
	}

	/* Initialize cone handler */
	(*intensity_init)(radiation_type, sycamera_polarization, lambdas, spectrum_resolution, integral_resolution);

	/* Compute tangent of half vision angle for future use */
	if (visang <= 0 || visang >= PI) {
		fprintf(stderr, "Invalid value of vision angle: %e\n", visang);
		exit(-1);
	}
	tanvisang = tan(visang/2);

	/* Make sure ddet is normalized */
	double ddetnorm = sqrt(ddet->val[0]*ddet->val[0] + ddet->val[1]*ddet->val[1] + ddet->val[2]*ddet->val[2]);

	if (ddetnorm == 0.) {
		fprintf(stderr, "The detector normal must be set!\n");
		exit(-1);
	}

	ddet->val[0] /= ddetnorm;
	ddet->val[1] /= ddetnorm;
	ddet->val[2] /= ddetnorm;

	/* Calculate unit vectors in canvas plane */
	e1 = vinit(3, 0.0, 0.0, 0.0); e2 = vinit(3, 0.0, 0.0, 0.0);
	if (ddet->val[1] == 0) e1->val[1] = 1;
	else {
		double dd = sqrt(ddet->val[0]*ddet->val[0] + ddet->val[1]*ddet->val[1]);
		e1->val[0] = ddet->val[1]/dd;
		e1->val[1] = -ddet->val[0]/dd;
	}

	e2->val[0] = ddet->val[2]*e1->val[1];
	e2->val[1] =-ddet->val[2]*e1->val[0];
	e2->val[2] = ddet->val[1]*e1->val[0] - ddet->val[0]*e1->val[1];

	/* Print info about image orientation */
	printf("Camera image y-axis is (%.3f, %.3f, %.3f)\n", e2->val[0], e2->val[1], e2->val[2]);

	/* Initialize toroidal angle maps, which are used when rotating
	 * particles toroidally. */
	sycamera_costor = malloc(sizeof(double)*toroidal_resolution);
	sycamera_sintor = malloc(sizeof(double)*toroidal_resolution);

	if (toroidal_resolution <= 1) {
		fprintf(stderr, "ERROR: Toroidal resolution must be greater than one.\n");
		exit(1);
	}

	/* Usually we would divide dPhi by (toroidal_resolution-1), but since
	 * the initial and final points would then be equivalent, we divide
	 * by toroidal_resolution directly to obtain the interval [-PI,PI) */
	double angle=-PI;
	sycamera_dPhi = 2.0*PI/toroidal_resolution;
	for (i = 0; i < toroidal_resolution; i++, angle += sycamera_dPhi) {
		sycamera_costor[i] = cos(angle);
		sycamera_sintor[i] = sin(angle);
	}

	sycamera_lasti = sycamera_lastj = -1;
	sycamera_lastlx = sycamera_lastly = 0;
}

void sycamera_init_run(unsigned int variables) {
	/* Initialize vector solutions */
	sol1 = vnew(variables);
	sol2 = vnew(variables);

	/* Initialize temporary vectors */
	temps = malloc(sizeof(vector)*NUMBER_OF_TEMPS);
	int i;
	for (i = 0; i < NUMBER_OF_TEMPS; i++) {
		temps[i].val = malloc(sizeof(double)*3);
		temps[i].n = 3;
	}

	/* Initialize camera image */
	sycout_prepare_run();

	if (sycamera_has_distfunc) distfunc_init_run();

	(*intensity_init_run)();
}
ode_solution *sycamera_init_particle(particle *p) {
	lasttime = p->t0;
	sycamera_charge = p->charge;
	sycamera_mass = p->mass;
	sycamera_particle_diffel = p->diffel;

	/* Initialize camera map particle */
	sycout_init_particle(p);
	(*intensity_init_particle)(p);

	/* Create ODE solution object */
	ode_solution *solver_object;
	solver_object = malloc(sizeof(ode_solution));
	solver_object->Z = sol1;
	solver_object->result = sol2;
	solver_object->step = 1e-8;

	if (sycamera_has_distfunc) {
		double r = hypot(p->r0[0], p->r0[1]);
		double v2 = p->vpar*p->vpar + p->vperp*p->vperp;
		double v = sqrt(v2);
		double gamma = 1/sqrt(1 - v2 / (LIGHTSPEED*LIGHTSPEED));
		double momentum = gamma * p->mass * v;
		double costheta = fabs(p->vpar / v);

		sycamera_distfunc_weight = distfunc_eval(r, costheta, momentum);
	}

	return solver_object;
}
void sycamera_deinit_run(void) {
	/* Assemble image (if any) */
	sycout_deinit_run();
}
/**
 * In order to utilize the toroidal symmetry correctly,
 * the orbit should always be computed for one poloidal
 * revolution. So always return FALSE here.
 */
int sycamera_stop_condition(void) {
	return 0;
}

void sycamera_step(ode_solution *solver_object, step_data *sd) {
	double x0=sd->x, y0=sd->y, vx=sd->vx, vy=sd->vy, r = hypot(x0, y0);

	(*intensity_init_step)(sd);

	int i;
	for (i = 0; i < toroidal_resolution; i++) {
		sd->x = x0*sycamera_costor[i] + y0*sycamera_sintor[i];
		sd->y =-x0*sycamera_sintor[i] + y0*sycamera_costor[i];

		sd->vx = vx*sycamera_costor[i] + vy*sycamera_sintor[i];
		sd->vy =-vx*sycamera_sintor[i] + vy*sycamera_costor[i];

		sycamera_process_step(sd, r*sycamera_dPhi);
	}

	/* Restore stepdata object */
	sd->x = x0;
	sd->y = y0;
	sd->vx = vx;
	sd->vy = vy;

	vector *t = solver_object->Z;
	solver_object->Z = solver_object->result;
	solver_object->result = t;
}

/**
 * Register radiation hitting the camera. This function
 * is called when it has been verified that some radiation
 * actually hits the camera.
 *
 * i: Pixel position in x, which is hit (continuous between 0 and 1)
 * j: Pixel position in y, which is hit (continuous between 0 and 1)
 * sd: Information about particle (mass, charge, position, momentum)
 * intensity: Fraction of radiation that hits detector
 * RdPhi: Cartesian length of toroidal step
 */
void sycamera_register_radiation(double i, double j, step_data *sd, double intensity, double RdPhi) {
	double differential;
	struct sycout_data data;
	data.sd = sd;
	//data.differential = RdPhi * sd->Jdtdrho * sycamera_particle_diffel;
	data.RdPhi = RdPhi;
	data.Jdtdrho = sd->Jdtdrho;
	data.particle_diffel = sycamera_particle_diffel;
	data.distribution_function = 1;

	if (sycamera_has_distfunc)
		data.distribution_function *= sycamera_distfunc_weight;

	data.i = i;
	data.j = j;

	differential = data.RdPhi * data.Jdtdrho * data.particle_diffel * data.distribution_function;
	data.brightness = intensity * differential;

	if (DEBUG_OUTPUT)
		print_particle_struct(sd);

	sycout_step(&data);
}

/**
 * Get the total intensity from one particle.
 */
double sycamera_get_intensity(
	step_data *sd, vector *rcp, vector *vhat, vector *empty1,
	vector *empty2, vector *empty3
) {
	return (*intensity_function)(sd, rcp, vhat, empty1, empty2, empty3);
}

/**
 * Find where radiation hits (if at all) and register
 * the radiation.
 *
 * sd: Step data object
 * dt: Length of time step
 */
int sycamera_process_step(step_data *sd, double RdPhi) {
	vector
		*l = temps, *ltilde = temps+1,
		*rcp = temps+2, *vhat = temps+3,
		*empty1 = temps+4, *empty2=temps+5,	/* These should always come last! */
		*empty3 = temps+6;

	/* Calculate r_{c-p} */
	rcp->val[0] = sd->x - Rdet->val[0];
	rcp->val[1] = sd->y - Rdet->val[1];
	rcp->val[2] = sd->z - Rdet->val[2];

	/* Calculate vhat */
	vhat->val[0] = sd->vx;
	vhat->val[1] = sd->vy;
	vhat->val[2] = sd->vz;

	double rcp_ddet  = vdot3(rcp, ddet);

	/********************************************/
	/*** CHECK WHETHER PARTICLE IS WITHIN FOV ***/
	/********************************************/
	/* Calculate vector from center of image to particle */
	ltilde->val[0] = rcp->val[0] - rcp_ddet*ddet->val[0];
	ltilde->val[1] = rcp->val[1] - rcp_ddet*ddet->val[1];
	ltilde->val[2] = rcp->val[2] - rcp_ddet*ddet->val[2];

	/* Rotate ltilde so that it has only two components */
	//l->val[0] = ddet->val[0]*ltilde->val[0] + e1->val[0]*ltilde->val[1] + e2->val[0]*ltilde->val[2];	/* Always = 0 */
	l->val[1] = ddet->val[1]*ltilde->val[0] + e1->val[1]*ltilde->val[1] + e2->val[1]*ltilde->val[2];
	l->val[2] = ddet->val[2]*ltilde->val[0] + e1->val[2]*ltilde->val[1] + e2->val[2]*ltilde->val[2];

	/* Calculate size of vision square */
	double img_side = sqrt(2) * fabs(rcp_ddet) * tanvisang,
		   img_side2 = img_side/2.0;

	/* Check if particle is within vision square */
	if (l->val[1] < -img_side2 || l->val[2] < -img_side2) return 0;
	if (l->val[1] > img_side2 || l->val[2] > img_side2) return 0;

	/* Normalize pixel indices */
	double i = (l->val[1]+img_side2)/img_side;
	double j = (l->val[2]+img_side2)/img_side;

	/* Determine whether the radiation hits, and if so register it */
	double intensity = sycamera_get_intensity(sd, rcp, vhat, empty1, empty2, empty3);

	if (intensity > 0.0) {
		/* Make sure the radiation is not blocked by any wall.
		   NOTE: This is a really heavy operation and should be done
		   after all other tests.
	 	*/
		if (domain_check3d(sd->x, sd->y, sd->z, Rdet->val[0], Rdet->val[1], Rdet->val[2]) == DOMAIN_OUTSIDE) {
			return 0;
		}

		//if (sycamera_lasti == i && sycamera_lastj == j && sycamera_lasti >= 0)
			//double r2 = rcp->val[0]*rcp->val[0] + rcp->val[1]*rcp->val[1] + rcp->val[2]*rcp->val[2];
			sycamera_register_radiation(i, j, sd, intensity, RdPhi);
		/*
		else
			sycamera_interp_pixels(
				l->val[1]+img_side2, l->val[2]+img_side2, abs(sycamera_lasti-i),
				abs(sycamera_lastj-j), img_side, sd, intensity, dt, dPhi
			);

		*/
		return 1;
	} else {
		return 0;
	}
}

void sycamera_output(equation *eq) {
	/* DEBUG 
	 * Output orbit */
	/*
	int i, j;
	fprintf(stderr, "T,X,Y,Z,~,~,~,~,~,~\n");
	for (i = 0; i < sycamera_last_orbit_index; i++) {
		j = i*sycamera_variables;
		fprintf(stderr,
			"%.12e,%.12e,%.12e,%.12e,0,0,0,0,0,0\n",
			sycamera_last_orbit_t[i],
			sycamera_last_orbit[j],
			sycamera_last_orbit[j+1],
			sycamera_last_orbit[j+2]
		);
	}
	*/

	/* If we use MPI we need to wait for an output signal
	 * before we start writing anything */
	int mpi_rank=0, nprocesses = 1;
#ifdef USE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
#endif

	sycout_output_all(mpi_rank, nprocesses);
}

/* Debugging functions */
void print_intersections(int tindex, double *x, double *y, double x0, double y0, double a, double b, double cosxi, double sinxi, long long int counter, double cms, step_data *sd, double xhit, double yhit, double ola, double cosphi) {
	int i;
	fprintf(stderr, "struct('i',%lld,'cms',%e,'time',%e,'pos',[%e,%e,%e],'vhat',[%e,%e,%e],'tanpitch',%e,'xhit',%e,'yhit',%e,'ola',%e,'cosphi',%e,'xy',[", counter, cms, sd->time,sd->x,sd->y,sd->z,sd->vx,sd->vy,sd->vz,sd->vperp/sd->vpar, xhit, yhit, ola, cosphi);
	for (i = 0; i < tindex; i++)
		fprintf(stderr, "%e,%e;", x[i], y[i]);

	double xi = acos(cosxi);
	if (sinxi < 0) xi = 2*PI - xi;
	if (tindex > 0) printf("\b");
	fprintf(stderr, "],'A',%e,'B',%e,'x0',%e,'y0',%e,'xi',%e);\n", a, b, x0, y0, xi);
}
void print_particle(step_data *p) {
	printf("r=[%e,%e,%e];", p->x, p->y, p->z);
	printf("v=[%e,%e,%e];", p->vx,p->vy,p->vz);
	printf("pitch=%e/%e;\n", p->vperp, p->vpar);
	//printf("  v0=%e, %e, %e\n", p->vx, p->vy, p->vz);
	//printf("  vpar/vperp=%e / %e\n", p->vpar, p->vperp);
	//printf("  r0=%e, %e, %e\n", p->x, p->y, p->z);
}
void print_particle_struct(step_data *p) {
	fprintf(stderr, "struct('r',[%.15e,%.15e,%.15e],'v',[%.15e,%.15e,%.15e],'pitch',%.15e/%.15e)\n",
		p->x, p->y, p->z,
		p->vx, p->vy, p->vz,
		p->vperp, p->vpar
	);
}
void print_particle_struct_p(step_data *sd) {
	double ptot = sqrt(sd->ppar2 + sd->pperp2);
	fprintf(stderr, "struct('r',[%.15e,%.15e,%.15e],'p',%e,'cospitch',%e,'B',%e)\n",
		sd->x, sd->y, sd->z,
		ptot, sqrt(sd->ppar2)/ptot,
		sd->B
	);
}
