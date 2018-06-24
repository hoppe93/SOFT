/* Particle generator */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "config.h"
#include "constants.h"
#include "domain.h"
#include "magnetic_axis.h"
#include "magnetic_field.h"
#include "global.h"
#include "particles.h"
#include "util.h"

#ifdef USE_MPI
#	include <mpi.h>
#endif

particle *particles_part;
/* Global boundary points for particles */
double particles_r0, particles_r1,
	   particles_param10, particles_param11,
	   particles_param20, particles_param21,
	   particles_dparam1, particles_dparam2, particles_dr,
       particles_effective_magnetic_axis;
enum particles_radial_coordinate particles_rcoord;

struct particles_paramspec *particles_param1spec, *particles_param2spec;

/* Variables used when particles are generated on each thread */
int particles_r=0, particles_param1=0, particles_param2=0,
	particles_rn=0,particles_param1n=0,particles_param2n=0,
	particles_done=0,
/* Local boundary points for particles, i.e. local to this thread */
	particles_rend, particles_param1end, particles_param2end,
	particles_count,
/* Variables to keep track of where we are */
	particles_curr_r, particles_curr_param1, particles_curr_param2;

/* Variables used when particles are generated centrally */
int particles_g_r=0, particles_g_param1=0, particles_g_param2=0,
	particles_g_done=0, particles_g_nextprogress=0,
	particles_g_progressteps=0, particles_g_dprogress=0,
	particles_g_progresscurr=0, particles_g_count=0,
/* Account for particle drifts with the 'rdyn' parameter */
	particles_g_drifts=0;

time_t particles_lastprogress=0;

/* Pre-computed differential element */
double particles_diffel, particles_rinner, particles_router,
	   particles_diffactor;

enum particles_generation_type particles_gentype=PARTICLES_GT_EVEN;

#define NPARAMS 7
struct particles_paramspec params[NPARAMS] = {
	{PARTICLES_COSPITCH, "xi", (PARTICLES_E|PARTICLES_P), 1.0},
	{PARTICLES_E, "E", (PARTICLES_COSPITCH|PARTICLES_MU|PARTICLES_PITCH|PARTICLES_PPAR|PARTICLES_PPERP), ENERGY},
	{PARTICLES_MU, "mu", (PARTICLES_E|PARTICLES_P|PARTICLES_PPAR), ENERGY},
	{PARTICLES_P, "p", (PARTICLES_COSPITCH|PARTICLES_MU|PARTICLES_PITCH|PARTICLES_PPAR|PARTICLES_PPERP), MOMENTUM},
	{PARTICLES_PITCH, "pitch", (PARTICLES_E|PARTICLES_P), 1.0},
	{PARTICLES_PPAR, "ppar", (PARTICLES_E|PARTICLES_P|PARTICLES_PPERP), MOMENTUM},
	{PARTICLES_PPERP, "pperp", (PARTICLES_E|PARTICLES_P|PARTICLES_PPAR), MOMENTUM}
};

#pragma omp threadprivate(particles_part,particles_r0,particles_r1,particles_param10,particles_param11, \
	particles_param20,particles_param21,particles_rend,particles_count,particles_param1end,particles_param2end, \
	particles_r,particles_param1,particles_param2,particles_dr,particles_dparam1,particles_dparam2, \
	particles_done,particles_param1spec,particles_param2spec,particles_diffel,particles_gentype, \
	particles_rcoord, particles_rinner, particles_router,particles_curr_r,particles_curr_param1, \
	particles_curr_param2, particles_diffactor,particles_effective_magnetic_axis)

void particles_set_param(double val0, double val1, int n, enum particles_inputtype type) {
	int i;
	struct particles_paramspec *pspec=NULL;
	for (i = 0; i < NPARAMS; i++) {
		if (type == params[i].type) {
			pspec = params+i;
			break;
		}
	}

	/* Set param 1? */
	if (particles_param1spec == NULL) {
		particles_param1spec = pspec;
		particles_param10 = val0*pspec->conversion;
		particles_param11 = val1*pspec->conversion;
		particles_param1n = n;
	} else {
		/* Otherwise, set param 2, right? */
		if (particles_param2spec != NULL)
			fprintf(stderr, "WARNING: More than two momentum space parameters set! The first and last ones set will be used.\n");

		particles_param2spec = pspec;
		particles_param20 = val0*pspec->conversion;
		particles_param21 = val1*pspec->conversion;
		particles_param2n = n;
	}
}
int verify_input_parameters(void) {
	/* Make sure the grid was specified correctly */
	if (particles_rcoord == PARTICLES_RC_MAJOR) {
		if (particles_rn == 1 && particles_r0 != particles_r1) {
			fprintf(stderr, "WARNING: Grid-size of 1 selected in r, but r0 != r1. "
				"The behaviour with this choice is undefined and most likely an error!\n");
		} else if (particles_rn > 1 && particles_r0 == particles_r1)
			fprintf(stderr, "WARNING: Grid-size greater than 1 selected, but r0 = r1.\n");
	}

	if (particles_param1n == 1 && particles_param10 != particles_param11) {
		fprintf(stderr, "WARNING: Grid-size of 1 selected in param1, but param10 != param11. "
			"The behaviour with this choice is undefined and most likely an error!\n");
	} else if (particles_param1n > 1 && particles_param10 == particles_param11)
		fprintf(stderr, "WARNING: Grid-size greater than 1 selected, but param10 = param11.\n");

	if (particles_param2n == 1 && particles_param20 != particles_param21) {
		fprintf(stderr, "WARNING: Grid-size of 1 selected in param2, but param20 != param21. "
			"The behaviour with this choice is undefined and most likely an error!\n");
	} else if (particles_param2n > 1 && particles_param20 == particles_param21)
		fprintf(stderr, "WARNING: Grid-size greater than 1 selected, but param20 = param21.\n");

	/* Have both coordinates been set? */
	if (particles_param1spec == NULL) {
		fprintf(stderr, "ERROR: The first momentum space parameter has not been set!\n");
		return 0;
	} else if (particles_param2spec == NULL) {
		fprintf(stderr, "ERROR: The second momentum space parameter has not been set!\n");
		return 0;
	}

	/* Are the coordinates compatible? */
	if ((particles_param1spec->compatible & particles_param2spec->type)==0) {
		fprintf(stderr,
			"ERROR: Parameter 1 is '%s' which is incompatible with parameter 2 (%s)\n",
			particles_param1spec->name,
			particles_param2spec->name
		);
		return 0;
	} else return 1;
}

/**
 * Initialize particle generator.
 *
 * t0:     Initial time used as reference
 * tend:   Final time
 * mass:   Mass of each particle (in AMU)
 * charge: Charge of each particle (in elementary charge units)
 * r0:     First radial coordinate to drop particle at
 * r1:     Final radial coordinate to drop particle at
 * rn:     Number of particles to drop between rstart and rstop
 * ppar0:  First parallel momentum point to drop particle at
 * ppar1:  Final parallel momentum point to drop particle at
 * pparn:  Number of points of parallel momentum to generate
 * pperp0: First perpendicular momentum point to drop particle at
 * pperp1: Final perpendicular momentum point to drop particle at
 * pperpn: Number of point of perpendicular momentum to generate
 */
void particles_init(struct particlespec *spec) {
	particles_gentype = spec->gentype;
	particles_part = malloc(sizeof(particle));

	if (omp_get_thread_num() == 0) {
		if (spec->gentype == PARTICLES_GT_QUEUE)
			fprintf(stdout, "Generating particles in 'queue' mode\n");
		else
			fprintf(stdout, "Generating particles evenly among threads and processes\n");
	}

	particles_param1spec = particles_param2spec = NULL;

	particles_part->t0 = spec->t0;
	particles_part->tend = spec->tend;
	particles_part->mass = spec->mass;
	particles_part->charge = spec->charge;
	particles_part->v0 = malloc(sizeof(double)*3);
	particles_part->r0 = malloc(sizeof(double)*3);
	particles_part->gc_position = spec->gc_position;
	particles_part->zeta0 = spec->zeta0;

	domain_get_bounds(&particles_router, &particles_rinner, NULL, NULL);

	particles_rcoord = spec->rcoord; particles_rn = spec->rn;
	particles_g_drifts = spec->include_drifts;
	if (particles_rcoord == PARTICLES_RC_MAJOR) {
		particles_r0 = spec->r0;
		particles_r1 = spec->r1;
	} else {
		particles_r1 = spec->rdynmax;
		if (particles_rn == 1) {
			particles_r0 = particles_r1;
			particles_rcoord = PARTICLES_RC_MAJOR;
		} else particles_r0 = magnetic_axis_r;
	}

	if (spec->inputtype & PARTICLES_PPAR)
		particles_set_param(spec->ppar0, spec->ppar1, spec->pparn, PARTICLES_PPAR);
	if (spec->inputtype & PARTICLES_PPERP)
		particles_set_param(spec->pperp0, spec->pperp1, spec->pperpn, PARTICLES_PPERP);
	if (spec->inputtype & PARTICLES_PITCH)
		particles_set_param(spec->pitch0, spec->pitch1, spec->pitchn, PARTICLES_PITCH);
	if (spec->inputtype & PARTICLES_COSPITCH)
		particles_set_param(spec->cospitch0, spec->cospitch1, spec->cospitchn, PARTICLES_COSPITCH);
	if (spec->inputtype & PARTICLES_P)
		particles_set_param(spec->p0, spec->p1, spec->pn, PARTICLES_P);

	/* Make sure the given input is valid */
	if (!verify_input_parameters())
		exit(1);

	/* Compute next point to report progress at */
	int totalpoints = particles_rn*particles_param1n*particles_param2n;
	if (particles_g_progressteps > 0) {
		particles_g_dprogress = totalpoints / particles_g_progressteps;
		particles_g_nextprogress = particles_g_dprogress;
		if (particles_gentype == PARTICLES_GT_EVEN)
			fprintf(stderr, "WARNING: Progress can only be reported when particles are launched in 'queue' mode.\n");

		particles_lastprogress = time(NULL);
	} else particles_g_nextprogress = totalpoints+1;	/* Never report progress */

	/* Compute differential element */
	particles_diffel = 1.;
	if (particles_rn <= 1 || particles_r0 == particles_r1) particles_dr = 0.;
	else {
		particles_dr = (particles_r1-particles_r0)/(particles_rn-1);
		/* This differential element enters into the factor
		 * Jdtdrho which is computed numerically in the
		 * orbit follower.
		 */
		//particles_diffel *= fabs(particles_dr);
	}

	if (particles_param1n <= 1 || particles_param10 == particles_param11) particles_dparam1 = 0.;
	else {
		particles_dparam1 = (particles_param11-particles_param10)/(particles_param1n-1);
		particles_diffel *= fabs(particles_dparam1);
	}

	if (particles_param2n <= 1 || particles_param20 == particles_param21) particles_dparam2 = 0.;
	else {
		particles_dparam2 = (particles_param21-particles_param20)/(particles_param2n-1);
		particles_diffel *= fabs(particles_dparam2);
	}

	particles_count = 0;
}

/**
 * This function can be called after initialization of the particles module.
 */
int particles_gridsize(void) {
	return particles_rn*particles_param1n*particles_param2n;
}
void particles_gridsize3(int grid[3]) {
	grid[0] = particles_rn;
	grid[1] = particles_param1n;
	grid[2] = particles_param2n;
}
void particles_local_gridsize(int grid[3]) {
	if (particles_gentype == PARTICLES_GT_QUEUE) {
		grid[0] = particles_rn;
		grid[1] = particles_param1n;
		grid[2] = particles_param2n;
	} else {
		grid[0] = particles_rn;
		grid[1] = particles_param1n;
		grid[2] = particles_param2n;
	}
}
/**
 * This function can be called at any point during execution.
 */
int particles_get_gridsize(struct particlespec *spec) {
	int gs = spec->rn;
	if (spec->inputtype & PARTICLES_PPAR) gs *= spec->pparn;
	if (spec->inputtype & PARTICLES_PPERP) gs *= spec->pperpn;
	if (spec->inputtype & PARTICLES_PITCH) gs *= spec->pitchn;
	if (spec->inputtype & PARTICLES_COSPITCH) gs *= spec->cospitchn;
	if (spec->inputtype & PARTICLES_P) gs *= spec->pn;

	return gs;
}

void particles_get_bounds6(int bounds[6]) {
	if (particles_gentype == PARTICLES_GT_QUEUE) {
		bounds[0] = particles_g_r;
		bounds[1] = particles_rend;
		bounds[2] = particles_g_param1;
		bounds[3] = particles_param1end;
		bounds[4] = particles_g_param2;
		bounds[5] = particles_param2end;
	} else {
		bounds[0] = particles_r;
		bounds[1] = particles_rend;
		bounds[2] = particles_param1;
		bounds[3] = particles_param1end;
		bounds[4] = particles_param2;
		bounds[5] = particles_param2end;
	}
}

/**
 * Compute where to start
 *
 * threadid: ID of this thread (unique among ALL running threads, not just within this process)
 * nthreads: Number of total threads ( = Number of MPI processes * threads per process)
 * mpi_rank: ID of this MPI process
 * nprocesses: Number of initiated MPI processes
 */
void particles_init_run(int threadid, int nthreads, int mpi_rank, int nprocesses) {
#ifndef USE_MPI	/* If NOT using MPI */
	if (particles_gentype == PARTICLES_GT_QUEUE) return;
#endif

	/* Compute maximum parameter index, chunk size and finally
	 * the parameter indices to be scanned by this thread */
	int maxparind = particles_gridsize();
	int parindchunk=maxparind / nthreads;
	int parindmod = maxparind % nthreads;

#ifdef USE_MPI
	if (particles_gentype == PARTICLES_GT_QUEUE) {
		maxparind = particles_gridsize();
		parindchunk=maxparind / nprocesses;
		parindmod = maxparind % nprocesses;
	}
#endif

	/* Check if the problem cannot be divided evenly among
	 * threads. If so, and if this is the first thread,
	 * warn about the resolution being ill-sized. It should
	 * be possible to distribute the remaining particles
	 * across other threads, but it requires a bit more
	 * thinking. For now, just warn the user. It shouldn't
	 * impact performance much anyway.
	 */
	if (parindmod != 0 && threadid == 0) {
		fprintf(stderr,
			"WARNING: The phase-space resolution should be a multiple of "
			"the number of CPU nodes and cores used. Currently %d particles "
			"are left-over and will be given to the last thread.\n",
			parindmod
		);
	}

	int parindstart = parindchunk*threadid;
	int parindend = parindstart + parindchunk - 1;

#ifdef USE_MPI
	if (particles_gentype == PARTICLES_GT_QUEUE)
		parindstart = parindchunk*mpi_rank;
#endif

	particles_param2    = parindstart         / (particles_param1n*particles_rn);
	particles_param2end = parindend           / (particles_param1n*particles_rn);
	parindstart        -= particles_param2    * (particles_param1n*particles_rn);
	parindend          -= particles_param2end * (particles_param1n*particles_rn);

	particles_param1    = parindstart         / particles_rn;
	particles_param1end = parindend           / particles_rn;
	parindstart        -= particles_param1    * particles_rn;
	parindend          -= particles_param1end * particles_rn;

	particles_r    = parindstart;
	particles_rend = parindend;

	/* If this is the last thread (or process), do all the remaining particles */
	if (threadid+1 == nthreads
#ifdef USE_MPI
		|| (particles_gentype == PARTICLES_GT_QUEUE && mpi_rank+1 == nprocesses)
#endif
	) {
		particles_rend = particles_rn-1;
		particles_param1end = particles_param1n-1;
		particles_param2end = particles_param2n-1;
	}

    /* Effective axis defaults to real magnetic axis.
     * If drifts are enabled, this value will be updated
     * when calculated in 'particles_generate'. */
    particles_effective_magnetic_axis = magnetic_axis_r;
}

/**
 * Returns the total number of particles
 * launched by SOFT.
 */
int particles_numberof(void) {
	if (particles_gentype == PARTICLES_GT_EVEN)
		return particles_count;
	else return particles_g_count;
}

/**
 * Returns 1 when all particles have been
 * processed. Returns 0 otherwise. */
int particles_stop_condition(void) {
	if (particles_gentype == PARTICLES_GT_EVEN)
		return particles_done;
	else return particles_g_done;
}

double particles_get_differential_factor_current(void) {
	return particles_diffactor;
}
double particles_get_differential_factor(double param1, double param2) {
	double m = particles_part->mass, c = LIGHTSPEED;

	switch (particles_param1spec->type) {
		case PARTICLES_COSPITCH:
			if (particles_param2spec->type == PARTICLES_E) return m*m*c * param1*param1 * param2*sqrt(param2*param2/(m*m*c*c*c*c)-1);
			else if (particles_param2spec->type == PARTICLES_P) return param2*param2;
			break;
		case PARTICLES_E:
			if (particles_param2spec->type == PARTICLES_COSPITCH) return m*m*c * param2*param2 * param1*sqrt(param1*param1/(m*m*c*c*c*c)-1);
			else if (particles_param2spec->type == PARTICLES_PPAR) return /* XXX */1;
			else if (particles_param2spec->type == PARTICLES_PPERP) return /* XXX */1;
			else if (particles_param2spec->type == PARTICLES_PITCH) return /* XXX */1;
			break;
		case PARTICLES_P:
			if (particles_param2spec->type == PARTICLES_COSPITCH) return param1*param1;
			else if (particles_param2spec->type == PARTICLES_PITCH) return param1*param1 * sin(param2);
			else if (particles_param2spec->type == PARTICLES_PPAR) return param1;
			else if (particles_param2spec->type == PARTICLES_PPERP) return param1*param2 / sqrt(param1*param1 - param2*param2);
			break;
		case PARTICLES_PITCH:
			if (particles_param2spec->type == PARTICLES_E) return /* XXX */1;
			else if (particles_param2spec->type == PARTICLES_P) return param2*param2 * sin(param1);
			break;
		case PARTICLES_PPAR:
			if (particles_param2spec->type == PARTICLES_E) return /* XXX */1;
			else if (particles_param2spec->type == PARTICLES_PPERP) return param2;
			else if (particles_param2spec->type == PARTICLES_P) return param2;
			break;
		case PARTICLES_PPERP:
			if (particles_param2spec->type == PARTICLES_E) return /* XXX */1;
			else if (particles_param2spec->type == PARTICLES_PPERP) return param1;
			else if (particles_param2spec->type == PARTICLES_P) return param1*param2 / sqrt(param2*param2 - param1*param1);
			break;
		default: break;
	}

	return 1;
}

/**
 * Computes the determinant of the Jacobian for the
 * current particle (i.e. the integral differential element).
 */
double particles_get_differential_element(double param1, double param2) {
	double dV = particles_diffel;
	particles_diffactor = particles_get_differential_factor(param1, param2);

	return dV * particles_diffactor;
}

/**
 * This function generates particles so that
 * they are distributed uniformly across threads
 * and processes.
 */
particle *particles_generate_distributed(void) {
	double r0 = particles_r0;
	double r = particles_r0 + particles_dr*particles_r;
	double param1 = particles_param10 + particles_dparam1*particles_param1;
	double param2 = particles_param20 + particles_dparam2*particles_param2;

	particles_curr_r = particles_r;
	particles_curr_param1 = particles_param1;
	particles_curr_param2 = particles_param2;

	if (particles_rcoord == PARTICLES_RC_DYNAMIC) {
		if (particles_g_drifts) {
			double ppar, pperp, pabs2;
			particles_params2p(param1, param2, &ppar, &pperp, &pabs2);
			r0 = particles_find_axis_r(ppar, pperp);
			/* We check if r0 > r before returning later... */
		} else {
			r0 = magnetic_axis_r;
		}
	}

	/* Last particle for this thread? */
	if (particles_r == particles_rend-1 &&
		particles_param1 == particles_param1end-1 &&
		particles_param2 == particles_param2end-1)
		particles_done = 1;

	/* Move on to next particle */
	particles_r++;
	particles_count++;

	if (particles_r >= particles_rn) {
		particles_r = 0;
		particles_param1++;

		if (particles_param1 >= particles_param1n) {
			particles_param1 = 0;
			particles_param2++;

			if (particles_param2 >= particles_param2n) {
				particles_done = 1;
				particles_param2 = 0;
			}
		}
	}

	if (r0 <= r) {
		return particles_generate_at(r, param1, param2);
	} else if (particles_done) {
		return NULL;
	} else {
		return particles_generate_distributed();
	}
}
/**
 * This function generates particles such that
 * on each call the function returns the next
 * particle which has not been processed by any
 * thread yet. It is therefore impossible to tell
 * beforehand which particles go to which threads.
 * However, this approach makes the particles more
 * evenly distributed according to the time it
 * takes to run them.
 */
particle *particles_generate_queue(void) {
	double r, param1, param2, r0 = particles_r0;

	if (particles_g_done) return NULL;

	#pragma omp critical
	{
		r      = particles_r0 + particles_dr*particles_g_r;
		param1 = particles_param10 + particles_dparam1*particles_g_param1;
		param2 = particles_param20 + particles_dparam2*particles_g_param2;

		particles_curr_r = particles_g_r;
		particles_curr_param1 = particles_g_param1;
		particles_curr_param2 = particles_g_param2;

		if (particles_rcoord == PARTICLES_RC_DYNAMIC) {
			if (particles_g_drifts) {
				double ppar, pperp, pabs2;
				particles_params2p(param1, param2, &ppar, &pperp, &pabs2);
				r0 = particles_find_axis_r(ppar, pperp);
				/* We check if r0 > r before returning later... */
			} else {
				r0 = magnetic_axis_r;
			}
		}

		/* Last particle for this thread? */
		if (particles_g_r == particles_rn-1 &&
			particles_g_param1 == particles_param1n-1 &&
			particles_g_param2 == particles_param2n-1)
			particles_g_done = 1;

		/* Move on to next particle */
		particles_g_r++;
		particles_g_count++;

		if (particles_g_count == particles_g_nextprogress) {
			particles_g_progresscurr++;
			particles_reportprogress(
				particles_g_progresscurr, particles_g_progressteps,
				particles_g_count, particles_rn*particles_param1n*particles_param2n
			);
			particles_g_nextprogress += particles_g_dprogress;
		}

		if (particles_g_r >= particles_rn) {
			particles_g_r = 0;
			particles_g_param1++;

			if (particles_g_param1 >= particles_param1n) {
				particles_g_param1 = 0;
				particles_g_param2++;

				if (particles_g_param2 >= particles_param2n) {
					particles_g_done = 1;
					particles_g_param2 = 0;
				}
			}
		}
	}

	if (r0 <= r) {
		return particles_generate_at(r, param1, param2);
	} else if (particles_g_done)
		return NULL;
	else
		return particles_generate_distributed();
}
particle *particles_generate(void) {
	if (particles_gentype == PARTICLES_GT_EVEN)
		return particles_generate_distributed();
	else
		return particles_generate_queue();
}
/**
 * Convert input parameters to ppar and pperp
 */
void particles_params2p(
	double param1, double param2,
	double *ppar, double *pperp, double *pabs2
) {
	double m2 = particles_part->mass * particles_part->mass,
		   c2 = LIGHTSPEED * LIGHTSPEED, p;

	switch (particles_param1spec->type) {
		case PARTICLES_COSPITCH:
			if (particles_param2spec->type == PARTICLES_E) {
				*pabs2 = param2*param2/c2 - m2*c2;
				p = sqrt(*pabs2);
				*ppar = p * param1;
				*pperp= p * sqrt(1 - param1*param1);
			} else if (particles_param2spec->type == PARTICLES_P) {
				*ppar = param2 * param1;
				*pperp = param2 * sqrt(1 - param1*param1);
				*pabs2 = param2*param2;
			}
			break;

		case PARTICLES_E:
			*pabs2 = param1*param1/c2 - m2*c2;
			p = sqrt(*pabs2);
			if (particles_param2spec->type == PARTICLES_COSPITCH) {
				*pperp = p * sqrt(1 - param2*param2);
				*ppar = p * param2;
			/*} else if (particles_param2spec->type == PARTICLES_MU) {
				double gmm = sqrt(1+pabs2/(m2*c2));
				*pperp = sqrt(2.0 * param2 * gmm * )
			*/
			} else if (particles_param2spec->type == PARTICLES_PITCH) {
				*pperp = p * sin(param2);
				*ppar = p * cos(param2);
			} else if (particles_param2spec->type == PARTICLES_PPAR) {
				*ppar = param2;
				*pperp = sqrt(*pabs2 - (*ppar)*(*ppar));
			} else if (particles_param2spec->type == PARTICLES_PPERP) {
				*pperp = param2;
				*pperp = sqrt(*pabs2 - (*pperp)*(*pperp));
			}
			break;
		case PARTICLES_P:
			*pabs2 = param1 * param1;
			if (particles_param2spec->type == PARTICLES_COSPITCH) {
				*pperp = param1 * sqrt(1 - param2*param2);
				*ppar = param1 * param2;
			/*} else if (particles_param2spec->type == PARTICLES_MU) {
				double gmm = sqrt(1+pabs2/(m2*c2));
				*pperp = sqrt(2.0 * param2 * gmm * )
			*/
			} else if (particles_param2spec->type == PARTICLES_PITCH) {
				*pperp = param1 * sin(param2);
				*ppar = param1 * cos(param2);
			} else if (particles_param2spec->type == PARTICLES_PPAR) {
				*ppar = param2;
				*pperp = sqrt(*pabs2 - (*ppar)*(*ppar));
			} else if (particles_param2spec->type == PARTICLES_PPERP) {
				*pperp = param2;
				*pperp = sqrt(*pabs2 - (*pperp)*(*pperp));
			}
			break;

		case PARTICLES_PITCH:
			if (particles_param2spec->type == PARTICLES_E) {
				*pabs2 = param2*param2/c2 - m2*c2;
				p = sqrt(*pabs2);
				*ppar = p * cos(param1);
				*pperp= p * sin(param1);
			} else if (particles_param2spec->type == PARTICLES_P) {
				*ppar = param2 * cos(param1);
				*pperp = param2 * sin(param1);
				*pabs2 = param2*param2;
			}
			break;

		case PARTICLES_PPAR:
			*ppar = param1;
			if (particles_param2spec->type == PARTICLES_E) {
				*pabs2 = param2*param2/c2 - m2*c2;
				*pperp = sqrt(*pabs2 - param1*param1);
			} else if (particles_param2spec->type == PARTICLES_PPERP) {
				*pperp = param2;
				*pabs2 = (*ppar)*(*ppar) + (*pperp)*(*pperp);
			} else if (particles_param2spec->type == PARTICLES_P) {
				*pabs2 = param2*param2;
				*pperp = sqrt(*pabs2 - (*ppar)*(*ppar));
			}
			break;

		case PARTICLES_PPERP:
			*pperp = param1;
			if (particles_param2spec->type == PARTICLES_PPAR) {
				*ppar = param2;
				*pabs2 = (*ppar)*(*ppar) + (*pperp)*(*pperp);
			} else if (particles_param2spec->type == PARTICLES_P) {
				*pabs2 = param2*param2;
				*ppar  = sqrt(*pabs2 - (*pperp)*(*pperp));
			}
			break;

		/* Make compiler shut up */
		default: break;
	}

	if (*pperp < 0) *pperp = fabs(*pperp);
}
particle *particles_generate_at(double r, double param1, double param2) {
	/* Compute pabs2, ppar and pperp */
	double pabs2=0., ppar=0., pperp=0.;
	particles_params2p(param1, param2, &ppar, &pperp, &pabs2);

	/* Here we should compute the most optimal location to put the particle
	 * in so that the particle's radiation can never hit the synchrotron
	 * camera.
	 */
	particles_part->r0[0] = r;
	particles_part->r0[1] = 0;
	particles_part->r0[2] = magnetic_axis_z;
	particles_part->ir    = particles_curr_r;
	particles_part->iv1   = particles_curr_param1;
	particles_part->iv2   = particles_curr_param2;

	particles_part->diffel = particles_get_differential_element(param1, param2);

	/* Compute gamma * m */
	double gm = sqrt(particles_part->mass*particles_part->mass + pabs2/(LIGHTSPEED*LIGHTSPEED));
	vpp2v(ppar, pperp, particles_part->r0, particles_part->v0);

	particles_part->v0[0] /= gm;
	particles_part->v0[1] /= gm;
	particles_part->v0[2] /= gm;
	particles_part->vpar = ppar/gm;
	particles_part->vperp = pperp/gm;

	return particles_part;
}

double *particles_get_bounds(void) {
	double *bounds = malloc(sizeof(double)*9);
	bounds[0] = particles_r0;
	bounds[1] = particles_r1;
	bounds[2] = particles_rn;
	bounds[3] = particles_param10;
	bounds[4] = particles_param11;
	bounds[5] = particles_param1n;
	bounds[6] = particles_param20;
	bounds[7] = particles_param21;
	bounds[8] = particles_param2n;

	return bounds;
}

double particles_get_drho(void) { return particles_dr; }
double particles_get_dvel1(void) { return particles_dparam1; }
double particles_get_dvel2(void) { return particles_dparam2; }

double particles_compute_X_pol(double ppar, double pperp, double r) {
	double Beffx, Beffy, Beffz, Beffpar, Xdotr, Xdotz;
	diff_data *dd;

	dd = magnetic_field_diff(r, 0.0, magnetic_axis_z);

	Beffx = dd->B->val[0] - ppar/CHARGE * dd->curlB->val[0],
	Beffy = dd->B->val[1] - ppar/CHARGE * dd->curlB->val[1],
	Beffz = dd->B->val[2] - ppar/CHARGE * dd->curlB->val[2];
	Beffpar = fabs((Beffx*dd->B->val[0] + Beffy*dd->B->val[1] + Beffz*dd->B->val[2]) / dd->Babs);

	Xdotr = ppar*Beffx + pperp*pperp / (2.0*CHARGE*dd->Babs*dd->Babs) * (dd->gradB->val[1]*dd->B->val[2] - dd->gradB->val[2]*dd->B->val[1]),
	Xdotz = ppar*Beffz + pperp*pperp / (2.0*CHARGE*dd->Babs*dd->Babs) * (dd->gradB->val[0]*dd->B->val[1] - dd->gradB->val[1]*dd->B->val[0]);

	return hypot(Xdotr, Xdotz) / Beffpar;
}
double particles_find_axis_r(double ppar, double pperp) {
	double a = particles_rinner, b = particles_router, c, d, tol = 1e-7,
		   Xpol_c=0, Xpol_d=0,
		   phi = 2 / (1+sqrt(5));

	/* Determine c and d */
	c = b + (a-b)*phi;
	d = a + (b-a)*phi;
	Xpol_c = particles_compute_X_pol(ppar, pperp, c);
	Xpol_d = particles_compute_X_pol(ppar, pperp, d);
		
	while (fabs(c-d) > tol) {
		if (Xpol_c < Xpol_d) {
			b = d;
			d = c;
			//Xpol_b = Xpol_d;
			Xpol_d = Xpol_c;
			c = b + (a-b)*phi;
			Xpol_c = particles_compute_X_pol(ppar, pperp, c);
		} else {
			a = c;
			c = d;
			//Xpol_a = Xpol_c;
			Xpol_c = Xpol_d;
			d = a + (b-a)*phi;
			Xpol_d = particles_compute_X_pol(ppar, pperp, d);
		}
	}

    particles_effective_magnetic_axis = (a+b)*0.5;
    /*
    fprintf(
        stderr, "%.16e,%.16e,%.16e\n",
        ppar/(ELECTRONMASS*LIGHTSPEED), pperp/(ELECTRONMASS*LIGHTSPEED),
        particles_effective_magnetic_axis
    );
    */

	return particles_effective_magnetic_axis;
}
/**
 * Returns the r-coordinate of the "effective" magnetic
 * axis (i.e. the point at which a guiding-center orbit
 * is just a point in the poloidal plane). If orbit
 * drifts are disabled, this is exactly the physical
 * magnetic axis.
 */
double particles_get_effective_magnetic_axis_r(void) {
    return particles_effective_magnetic_axis;
}
/**
 * Returns the distance by which the guiding-center orbit
 * is shifted away from a flux surface due to drifts.
 * (NOTE: This value can be negative!)
 */
double particles_get_orbit_drift_shift(void) {
    return (particles_get_effective_magnetic_axis_r()-magnetic_axis_r);
}

char *particles_param1_name(void) {
	return particles_param1spec->name;
}
char *particles_param2_name(void) {
	return particles_param2spec->name;
}

void particles_set_progress(int progsteps) {
	particles_g_progressteps = progsteps;
}
void particles_reportprogress(int curr, int totsteps, int parts, int ptot) {
#define SECONDS_PER_DAY (60*60*24)
#define SECONDS_PER_HOUR (60*60)
#define SECONDS_PER_MINUTE (60)
	time_t t = time(NULL);

	/* Compute time difference */
	time_t dt = t-particles_lastprogress;
	int days = dt / SECONDS_PER_DAY; dt %= SECONDS_PER_DAY;
	int hours = dt / SECONDS_PER_HOUR; dt %= SECONDS_PER_HOUR;
	int mins = dt / SECONDS_PER_MINUTE; dt %= SECONDS_PER_MINUTE;
	int secs = dt;

	if (days > 0) printf("PROGRESS (%d/%d): %d days, %d hours, %d minutes and %d seconds since last -- %d particles of %d done.\n", curr, totsteps, days, hours, mins, secs, parts, ptot);
	else if (hours > 0) printf("PROGRESS (%d/%d): %d hours, %d minutes and %d seconds since last -- %d particles of %d done.\n", curr, totsteps, hours, mins, secs, parts, ptot);
	else if (mins > 0) printf("PROGRESS (%d/%d): %d minutes and %d seconds since last -- %d particles of %d done.\n", curr, totsteps, mins, secs, parts, ptot);
	else printf("PROGRESS (%d/%d): %d seconds since last -- %d particles of %d done.\n", curr, totsteps, secs, parts, ptot);

	particles_lastprogress = t;
}

