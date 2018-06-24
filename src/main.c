/* The main program */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#include "config.h"
#include "constants.h"
#include "ctsv.h"
//#include "diag.h"
#include "domain.h"
#include "equations.h"
#include "global.h"
#include "interp2.h"
#include "IO_data.h"
#include "global.h"
#include "magnetic_axis.h"
#include "magnetic_field.h"
#include "particles.h"
#include "quantities.h"
#include "readfile.h"
#include "rkf45.h"
#include "settings.h"
#include "tools.h"
#include "util.h"

#ifdef USE_MPI
#	include <mpi.h>
#endif

#define EXIT_PARTICLE_OUTSIDE 11
#define EXIT_PARTICLE_START_OUTSIDE 12
#define EXIT_INVALID_EQUATION 13
#define EXIT_INVALID_MAGFIELD 14

unsigned int current_index=-1;
int return_value;

/* Interpolation of ODE timestep
 *   0 = Do not interpolate (linearly)
 *   n = Interpolate linearly with n points
 *       between all timesteps (n = int)
 *   x = Interpolate so that the timestep
 *       never exceeds x (x = double < 1)
 */
double interp_timestep = 0;

/* Toroidal check */
int signChanged=0;
double phi=0.;

/* Poloidal check */
double pol_lastr=0., pol_lastz=0.;
int pol_rpass=0, pol_zpass=0, pol_trapped=0, pol_rest=0;

#pragma omp threadprivate(current_index,signChanged,phi,pol_lastr,pol_rpass,pol_lastz,pol_zpass,pol_trapped,pol_rest)

/**
 * Checks whether the particle should not be simulated any more.
 * Returns 1 if simulation should stop, 0 otherwise.
 * If a negative end-time is given, this function will return
 * a termination signal if one toroidal turn has been traversed.
 *
 * t: Current time
 * tend: End time as specified in the pi-file (<0 means position condition should be used
 *       until t > -tend)
 * z0: Initial z position of the particle
 * z: Current z coordinate of the particle
 */
int stop_condition(double t, double tend, double r0, double r, double z0, double z) {
	if (tend < 0) {
		/* If we go above -tend we should cancel */
		if (t > -tend) return 1;

		/* Count number of times the particle passes
		 * the r=r0 and z=z0 lines respectively
		 * (where r0,z0 is the particle's initial position) */
		if ((pol_lastr < r0 && r0 < r) ||
			(pol_lastr > r0 && r0 > r)) {
			pol_rpass++;
		}

		if ((pol_lastz < z0 && z0 < z) ||
			(pol_lastz > z0 && z0 > z)) {
			pol_zpass++;
		}

		/* If the particle hasn't moved since last timestep... */
		if (pol_lastr == r && pol_lastz == z)
			pol_rest++;

		/* We do this to prevent particle's at rest from
		 * stopping execution */
		int ret = (pol_rest >= 2);
		//if (pol_rpass >= 2 || pol_zpass >= 2) {
		if (pol_zpass >= 2) {
			/*if (pol_trapped || fabs(r0-r) > fabs(pol_lastr-r) || fabs(z0-z) > fabs(pol_lastz-z)) {
				pol_trapped = 1;
				ret = (pol_rpass >= 4 || pol_zpass >= 4);
			} else ret = 1.;*/
			ret = 1;
		}

		pol_lastr = r;
		pol_lastz = z;
		return ret;
	} else return (t > tend);
}

/**
 * Take n linear steps between two ODE solutions.
 *
 * n: Number of steps to take (0 = y0 -> y1, 1 = y0 -> y_{1/2} -> y_1 etc.)
 * nvariables: Number of variables of equation
 * solver_object: Solver object containing current step and solution
 */
void take_n_steps(unsigned int n, unsigned int nvariables, ode_solution *solver_object, double Jdtdrho) {
	unsigned int i,j;
	double inv = 1 / (double)n,
		x = solver_object->Z->val[0],
		y = solver_object->Z->val[1],
		z = solver_object->Z->val[2],
		dx = (solver_object->result->val[0]-x) * inv,
		dy = (solver_object->result->val[1]-y) * inv,
		dz = (solver_object->result->val[2]-z) * inv,
		others[nvariables-3];

	for (i = 3; i < nvariables; i++)
		others[i-3] = solver_object->result->val[i];
	
	for (i = 0; i < n; i++) {
		solver_object->result->val[0] = solver_object->Z->val[0] + dx;
		solver_object->result->val[1] = solver_object->Z->val[1] + dy;
		solver_object->result->val[2] = solver_object->Z->val[2] + dz;

		for (j = 3; j < nvariables; j++)
			solver_object->result->val[j] = others[j-3];

		tool_update_position(
			solver_object->Z->val[0],
			solver_object->Z->val[1],
			solver_object->Z->val[2]
		);
		tool_step(solver_object, Jdtdrho);
	}
}

ode_solution *init_jacobian_solver_object(int nvars) {
	ode_solution *solver_object;
	solver_object = malloc(sizeof(ode_solution));
	solver_object->Z = vnew(nvars);
	solver_object->result = vnew(nvars);
	solver_object->step = 1e-8;

	return solver_object;
}

void deinit_jacobian_solver_object(ode_solution *so) {
	vfree(so->Z);
	vfree(so->result);
	free(so);
}

/**
 * Calculates saves the values of the particle motion with given parameters
 *
 * RETURNS the length of the average timestep taken
 */
double main_solve(particle *p, equation *eq, magnetic_handler *mh, tool *usetool) {
	double x,y,z,r;
	ode_solution *solver_object, *jacobian_so=NULL;

	/***************************/
	/* Initialize interpolator */
	/***************************/
	mh->init_particle();

	/* domain check of initial position */
	x = p->r0[0];
	y = p->r0[1];
	z = p->r0[2];

	r = hypot(x,y);
	double R[2] = {magnetic_axis_r,r};
	double Z[2] = {magnetic_axis_z,z};

	if (domain_check(R, Z) == DOMAIN_OUTSIDE) {
		printf("Particle starts outside device! (%e, %e, %e)\n", x, y, z);
		exit(EXIT_PARTICLE_START_OUTSIDE);
	}

	/* Initialize tool */
	solver_object = usetool->init_particle(p);
	/* Store initial values in solver_object */
	eq->init_run(p, solver_object);

	/* Prepare to compute the Jacobian */
	if (usetool->require_jacobian) {
		double r0[3], *tr0;
		//r0[0] = p->r0[0] + particles_get_drho();
		r0[0] = p->r0[0] + JACOBIAN_ORBIT_STEP;
		r0[1] = p->r0[1];
		r0[2] = p->r0[2];

		tr0 = p->r0;
		p->r0 = r0;

		jacobian_so = init_jacobian_solver_object(eq->variables);
		eq->init_run(p, jacobian_so);

		p->r0 = tr0;
	}

	/* Main loop. Loop until the final time has been reached */
	current_index = 0;
	double current_time = p->t0,
		dR_dt, dR_drho, dZ_dt, dZ_drho, Jdtdrho=1, drho = particles_get_drho();
	int steps = 0;
    
	/* Since the equation used may displace the solver object slightly (such as if we're
	 * following a GC orbit and specify the particle position), we should use 'solver_object->Z->val[...]'
	 * instead of 'p->r0[...]' to calculate the initial particle position */
	x = solver_object->Z->val[0];
	y = solver_object->Z->val[1];
	z = solver_object->Z->val[2];
	pol_lastr = hypot(x,y);
	pol_lastz = z;
	double r0 = pol_lastr, z0 = pol_lastz;
	while (!stop_condition(current_time, p->tend, r0, r, z0, solver_object->Z->val[2]) &&
		   !usetool->stop_condition()) {
		steps++;

		/* Take a step */
		if (usetool->require_jacobian) {
			do {
				ode_solve(eq->eq, jacobian_so, current_time);
			} while (jacobian_so->flag == REDO_STEP);

			solver_object->step = jacobian_so->actualstep;
			ode_solve(eq->eq, solver_object, current_time);

			dR_dt = hypot(solver_object->result->val[0], solver_object->result->val[1]) -
					hypot(solver_object->Z->val[0], solver_object->Z->val[1]);
			dZ_dt = solver_object->result->val[2] - solver_object->Z->val[2];
			dR_drho = hypot(jacobian_so->Z->val[0], jacobian_so->Z->val[1]) -
					hypot(solver_object->Z->val[0], solver_object->Z->val[1]);
			dZ_drho = jacobian_so->Z->val[2] - solver_object->Z->val[2];

			Jdtdrho = fabs(dR_dt*dZ_drho - dZ_dt*dR_drho)/JACOBIAN_ORBIT_STEP;
			if (drho != 0)
				Jdtdrho *= drho;
		} else {
			do {
				ode_solve(eq->eq, solver_object, current_time);
			} while (solver_object->flag == REDO_STEP);
		}

		/* Increase time */
		current_time += solver_object->actualstep;

		/* Linearly interpolate timesteps */
		if (interp_timestep >= 1) {	/* Take exactly n steps */
			unsigned int count = (unsigned int)interp_timestep;
			take_n_steps(count, eq->variables, solver_object, Jdtdrho);
		/* Take enough steps so that the stepsize < interp_timestep */
		} else if (interp_timestep > 0 && solver_object->step > interp_timestep) {
			unsigned int count = (unsigned int)(solver_object->step / interp_timestep);
			take_n_steps(count, eq->variables, solver_object, Jdtdrho);
		} else /* No interpolation */
			tool_step(solver_object, Jdtdrho);

		if (usetool->require_jacobian) {
			vector *t = jacobian_so->Z;
			jacobian_so->Z = jacobian_so->result;
			jacobian_so->result = t;
		}

		/* Check if new position is inside domain */
		x=solver_object->Z->val[0];
		y=solver_object->Z->val[1];
		z=solver_object->Z->val[2];

		R[0] = R[1]; Z[0] = Z[1];
		r = R[1] = hypot(x, y); Z[1] = z;
		if (domain_check(R, Z) == DOMAIN_OUTSIDE) {
			double v0 = hypot(p->v0[0], hypot(p->v0[1], p->v0[2])),
				   gm = 1 / sqrt(1 - v0*v0/(LIGHTSPEED*LIGHTSPEED)),
				   p0 = gm * p->mass * v0,
				   xi0= p->vpar / v0;
			printf("Particle collided with device wall!  r0 = %e, p0 = %e, xi0 = %e\n", r0, p0, xi0);
			return_value = EXIT_PARTICLE_OUTSIDE;
			break;
		}
	}

	if (usetool->require_jacobian) {
		deinit_jacobian_solver_object(jacobian_so);
	}

	pol_rpass = 0;
	pol_zpass = 0;
	pol_rest  = 0;

	free(solver_object);

	return (current_time / steps);
}

/**
 * Report progress
 */
void report_progress(int part, int total) {
	double perc = ((double)part)/((double)total)*100;
	int p = (int)(perc/2), i;

	putchar('[');
	for (i = 0; i < 50; i++) {
		if (i <= p) putchar('#');
		else putchar(' ');
	}

	printf("] %.2f%%\r", perc);
}

void run_particles(int threadid, int nthreads, int mpi_rank, int nprocesses, void *equa, void *mhndl, void *tl, void *st) {
	equation *eq = (equation*)equa;
	magnetic_handler *mh = (magnetic_handler*)mhndl;
	tool *usetool = (tool*)tl;
	settings *set = (settings*)st;

	/* Initialize internal diagnostics */
	//diag_init("diag.txt");

	/* Initialize equation and prepare tool for run */
	eq->init(set);
	tool_prepare_run();
	mh->init_run();
	usetool->init_run(eq->variables);

	/*************************/
	/* Initialize ODE solver */
	/*************************/
	ode_init(eq->variables);

	/*********************************/
	/* Initialize particle generator */
	/*********************************/
	particles_init_run(threadid, nthreads, mpi_rank, nprocesses);

	struct timeval t1, t2;
	double avtime=0, avstep = 0, partime=0;
	while (!particles_stop_condition()) {
		gettimeofday(&t1, NULL);
		/* Solve */
		particle *p = particles_generate();
		if (p == NULL) continue;
		avstep += main_solve(p, eq, mh, usetool);

		gettimeofday(&t2, NULL);

		partime = (t2.tv_sec - t1.tv_sec) * 1000.0;
		partime += (t2.tv_usec - t1.tv_usec) / 1000.0;

		avtime += partime;
	}
	
	usetool->deinit_run();
	//diag_deinit();

	int totalparts = particles_numberof();
	printf("Average time per particle on thread #%d: %.3f ms, %.3e timestep\n",
		threadid, (avtime/((double)totalparts)), (avstep/((double)totalparts)));
}

/**
 * main function: stores initial values in global variable initial of
 * type initial_data. Also stores domain, magnetic field
 * and initializes interpolation of the magnetic field.
 *
 * calls main_solve.
 * uses ctsv_write to write solution data to file.
 */
int main(int argc, char *argv[]) {
	settings *set;

	/* Make sure output is not buffered */
	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);

	return_value = 0;		/* Set default return value (SUCCESS) */

	if (argc == 1) {
		/* Load settings from stdin */
		set = load_settings(NULL);
	} else if (argc == 2) {
		/* Load settings file */
		set = load_settings(argv[1]);
	} else {
		fprintf(stderr, "Expected only name of pi file as argument!\n");
		exit(-1);
	}

	/* Load domain */
	//dom = domain_load(set->domain);

	/*************************************/
	/* Initialize magnetic field handler */
	/*************************************/
	magnetic_init();

	/* Set which type of magnetic field to use */
	magnetic_handler *mh = magnetic_handler_select(set->magfield);
	if (mh == NULL) {
		printf("Invalid magnetic field handler selected!\n");
		exit(EXIT_INVALID_MAGFIELD);
	}

	/*******************************/
	/* Initialize equation handler */
	/*******************************/
	equation_handler_init();

	/* Set which problem to use from now on, GCM or regular particle motion */
	equation *eq = select_equation(set->equation);
	if (eq == NULL) {
		printf("Invalid equation selected!\n");
		exit(EXIT_INVALID_EQUATION);
	}

	int i;
	/* Find the right magnetic field handler settings to use */
	for (i = 0; i < set->nmagnetic; i++) {
		if (!strcmp(set->magnetic[i].name, set->magfield)) {
			mh->init(set->magnetic+i);
			break;
		}
	}

	/* If no specific settings were given for the magnetic field handler,
	 * init it with empty settings (if the tool requires settings, this may
	 * (but usually won't) crash the program) */
	if (i == set->nmagnetic)
		mh->init(NULL);

	/***************************/
	/* Initialize tool handler */
	/***************************/
	tool_init_handler();

	tool *usetool = tool_select(set->tool);
	if (usetool == NULL) {
		printf("ERROR: Invalid tool selected!\n");
		exit(-1);
	}

	/* Find the right tool settings to use */
	for (i = 0; i < set->ntools; i++) {
		if (!strcmp(set->tools[i].name, set->tool)) {
			usetool->init(set, set->tools+i, set->sycouts, set->nsycouts);
			break;
		}
	}

	/* If no specific settings were given for the tool, init it with empty settings
	 * (if the tool requires settings, this may crash the program) */
	if (i == set->ntools)
		usetool->init(set, NULL, NULL, 0);

	/**************************************************/
	/* Specify tolerance of ODE Solver (if requested).*/
	/* Otherwise this is taken care of by the tool.   */
	/**************************************************/
	if (set->tolerance > 0.)
		ode_tolerance = set->tolerance;
	
	printf("Selected solver tolerance = %e\n", ode_tolerance);

	if (set->maxtimestep > 0.) {
		ode_maxtimestep = set->maxtimestep;
		printf("Maximum allowed timestep is %e s\n", ode_maxtimestep);
	}

	/**************************************************/
	/* Set interpolation scheme of timesteps          */
	/**************************************************/
	interp_timestep = set->interptimestep;

	int nprocesses = 1,		/* When not using MPI we just have 1 process */
		mpi_rank=0;

#ifdef USE_MPI
	/******************/
	/* Initialize MPI */
	/******************/
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	/************************/
	/* Initialize particles */
	/************************/
	struct particlespec *spec = set->particlespec;

	/* Should progress be reported? (set->progress
	 * contains the number of times to report progress) */
	if (set->progress > 0)
		particles_set_progress(set->progress);

	/* Determine how many threads (and MPI processes) are needed */
	/* Without MPI, we can easily run fewer particles than the
	 * number of available threads. With MPI we will have to
	 * do a bit more math, and possibly exit this process
	 */
#ifdef USE_MPI
	if (particles_get_gridsize(spec) < set->threads * nprocesses) {
		int maxproc = particles_get_gridsize(spec)/set->threads;
		if (maxproc == 0) {
			set->threads = particles_get_gridsize(spec);
			maxproc = 1;
		} else
			set->threads = particles_get_gridsize(spec) - maxproc*set->threads;

		/* Is this MPI process needed? */
		if (mpi_rank > maxproc-1) {
			MPI_Finalize();
			exit(0);
		}

		if (mpi_rank == 0)
			printf("The number of particles are fewer than the number of threads! Limiting number of threads and processes to %d * %d\n", set->threads, maxproc);
	}
#else	/* Without MPI */
	if (particles_get_gridsize(spec) < set->threads) {
		set->threads = particles_get_gridsize(spec);
		printf("The number of particles are fewer than the number of threads! Limiting number of threads to %d\n", set->threads);
	}
#endif

	#pragma omp parallel num_threads(set->threads)
	{
		particles_init(spec);

		int threadid = omp_get_thread_num() + set->threads*mpi_rank;
		run_particles(threadid, set->threads*nprocesses, mpi_rank, nprocesses, eq, mh, usetool, set);
	}

	/* Write solution data to file specified in input file */
	usetool->output(eq);

#ifdef USE_MPI
	/* Finalize MPI environment */
	MPI_Finalize();
#endif

	return return_value;
}
