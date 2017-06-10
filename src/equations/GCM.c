/* Equation of Guiding-Center motion */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "equation_GCM.h"
#include "global.h"
#include "interp2.h"
#include "IO_data.h"
#include "magnetic_field.h"
#include "readfile.h"
#include "rkf45.h"
#include "settings.h"
#include "tools.h"
#include "vector.h"


/* Indices to the quantity array. Initialized
 * in `equation_GCM_init'. */
particle *GCM_particle=NULL;
/* Pre-allocated vectors */
vector *GCM_retval=NULL, **GCM_reservoir=NULL;
int GCM_reservoir_index = 0, GCM_reservoir_count=0, GCM_nodrifts=0;

#pragma omp threadprivate(GCM_particle,GCM_retval,GCM_reservoir,GCM_reservoir_index,GCM_reservoir_count,GCM_nodrifts)

void equation_GCM_init(settings *set) {
	/* Pre-allocate vectors */
	GCM_retval = vnew(5);

	GCM_nodrifts = set->nodrifts;
}

vector *equation_GCM_vnew(void) {
	if (GCM_reservoir_index+1 >= GCM_reservoir_count) {
		GCM_reservoir = realloc(GCM_reservoir, sizeof(vector*)*(GCM_reservoir_count+1));

		GCM_reservoir[GCM_reservoir_count] = vnew(3);
		GCM_reservoir_count++;
	}

	return GCM_reservoir[GCM_reservoir_index++];
}

/**
 * Function to initialize data for the GCM simulation run:
 * Stores initial values to solution vector
 * Converts particle initial values to GC initial values
 * Input: initial_data object, contains particle initial values
 * Input: solution vector
 * 
 * RETURNS pointer to ode_solution containing initial values for GC
 */
void equation_GCM_init_run(particle *part, ode_solution *solver_object) {
	GCM_particle = part;

	double m=part->mass; // particle mass
	double e = part->charge; // particle charge

	/* initial particle position */
	double x = part->r0[0],
	y = part->r0[1],
	z = part->r0[2];

	/* initial particle velocity */
	double vx = part->v0[0],
	vy = part->v0[1],
	vz = part->v0[2];

	/* Create initial velocity vector  */
	vector *v=equation_GCM_vnew();
	vinit(3,vx,vy,vz);

	/* Get magnetic field B for initial position */
	vector *B = magnetic_field_get(x,y,z);
		   
	double Bx=B->val[0],   
	By=B->val[1],
	Bz=B->val[2];	   
	/* Absolute value of B */
	double B_abs=sqrt(Bx*Bx + By*By + Bz*Bz);
		   
	/* Calculate bhat, unit vector in B-field direction */
	vector *bhat = vmulsf(1/B_abs,B);
	double bx=bhat->val[0],
	by=bhat->val[1],
	bz=bhat->val[2];

	/* Calculate absolute value of parallel velocity */
	double vpar_abs=vdot3(bhat,v);
			   
	/* negative parallel velocity vector, to calculate
	* perpendicular velocity vector */
	vector *bhat_vpar = vmuls(-vpar_abs,bhat, equation_GCM_vnew());
	/*  perpendicular velocity */
	vector *vperp= vadd(v,bhat_vpar, equation_GCM_vnew());
	double vperp_x=vperp->val[0], vperp_y=vperp->val[1], vperp_z=vperp->val[2];

	/* Allocate and calculate Larmor vector rho */
	vector *rho=equation_GCM_vnew();
	rho->val[0]=by*vperp_z-bz*vperp_y;
	rho->val[1]=bz*vperp_x-bx*vperp_z;
	rho->val[2]=bx*vperp_y-by*vperp_x;
	double Omega = e*B_abs/m; //angular velocity Omega (no 'c' in SI-units) 
	vmulsf(-1/Omega, rho);

	/* Calculate magnetic moment mu */
	double vperp_abs = vdot3(vperp, vperp);
	double mu=m*vperp_abs/(2*B_abs);

	if (part->gc_position) {
		solver_object->Z->val[0] = x;
		solver_object->Z->val[1] = y;
		solver_object->Z->val[2] = z;
		solver_object->Z->val[3] = vpar_abs;
		solver_object->Z->val[4] = mu;
	} else {
		/* Calculate guiding center position vector X */
		double X=x+rho->val[0], Y=y+rho->val[1], Z=z+rho->val[2];

		/* store initial values to solution vector */
		solver_object->Z->val[0] = X; // Guiding center x-position
		solver_object->Z->val[1] = Y; // Guiding center y-position
		solver_object->Z->val[2] = Z; // Guiding center z-position 
		solver_object->Z->val[3] = vpar_abs; 
		solver_object->Z->val[4] = mu;        // Magnetic moment mu
	}

	/* Initialize the solver object */
	/* Save initial data in solver object */
	solver_object->step = 1e-8; /* Initial step size */
	solver_object->flag=0;
}

/**
 * Initialize output data structure
 */
solution_data *equation_GCM_init_output(solution_data *output) {
	output->nvars = 5;
	output->labels = malloc(sizeof(char*)*output->nvars);
	output->labels[0] = "X";
	output->labels[1] = "Y";
	output->labels[2] = "Z";
	output->labels[3] = "u";
	output->labels[4] = "mu";

	return output;
}

/**
 * Equation for the guiding center motion, GCM
 * 
 * T: time, not used here but needed for ode.c
 * Z: Pointer to vector containing particle position
 * coordinates in the first three values, and particle
 * velocity in the next three.
 *
 * RETURNS: vector of values of function f 
 *
 * Used as the first argument to ode_solve in ode.c if GCM is selected
 */
vector *equation_GCM_eq(double T, vector *Z) {
	particle *part = GCM_particle;
	double m = part->mass;
	double e = part->charge;

	double x=Z->val[0], y=Z->val[1], z=Z->val[2],
		   udot, Xdot1, Xdot2, Xdot3, mu;

	mu = Z->val[4];

	/***************************
	* Calculate udot          *
	***************************/
	diff_data *dd = magnetic_field_diff(x,y,z);

	/* Calculate bhat */
	double bhatval[3];
	vector bhatv;
	bhatv.val = bhatval;
	bhatv.n = 3;
	vector *bhat = vmuls(1/dd->Babs, dd->B, &bhatv);

	/* Calculate B* (effective B-field) */
	vector *B_eff = vaddf(dd->B, vmulsf(m/e*Z->val[3], dd->curlB));

	/* Get value of B* parallel to b^ */
	double Beff_par = vdot3(B_eff, bhat);

	/* == > */
	udot = -mu/(m*Beff_par) * vdot3(dd->gradB, B_eff);

	if (GCM_nodrifts) {
		Xdot1 = Z->val[3] / m * bhat->val[0];
		Xdot2 = Z->val[3] / m * bhat->val[1];
		Xdot3 = Z->val[3] / m * bhat->val[2];
	} else {
		/***************************
		* Calculate Xdot1         *
		***************************/
		Xdot1 = 1/Beff_par * (-mu/e*(dd->gradB->val[1]*bhat->val[2] - dd->gradB->val[2]*bhat->val[1]) + Z->val[3] * B_eff->val[0]);
		/***************************
		* Calculate Xdot2         *
		***************************/
		Xdot2 = 1/Beff_par * (-mu/e*(dd->gradB->val[2]*bhat->val[0] - dd->gradB->val[0]*bhat->val[2]) + Z->val[3] * B_eff->val[1]);
		/***************************
		* Calculate Xdot3         *
		***************************/
		Xdot3 = 1/Beff_par * (-mu/e*(dd->gradB->val[0]*bhat->val[1] - dd->gradB->val[1]*bhat->val[0]) + Z->val[3] * B_eff->val[2]);
	}

	vector *retval = GCM_retval;
	retval->val[0] = Xdot1;
	retval->val[1] = Xdot2;
	retval->val[2] = Xdot3;
	retval->val[3] = udot;
	retval->val[4] = 0;

	/* Calculate energy */
	double v2 = Z->val[3]*Z->val[3] + 2*dd->Babs*mu/m;
	
	tool_record_data(
		T, m*m*udot*udot, m*m*(v2-udot*udot),
		dd->Babs,
		Z->val[0], Z->val[1], Z->val[2],
		Xdot1, Xdot2, Xdot3,
		Z->val[3], sqrt(2*dd->Babs*mu/m)
	);

	return retval;
}
