/* Equation of Guiding-Center motion */

#include <math.h>
#include <stdlib.h>
#include <omp.h>

#include "equation_GCM_rel.h"
#include "global.h"
#include "interp2.h"
#include "IO_data.h"
#include "magnetic_field.h"
#include "quantities.h"
#include "readfile.h"
#include "rkf45.h"
#include "tools.h"
#include "vector.h"

particle *GCM_rel_particle=NULL;
/* Pre-allocated vectors */
vector *GCM_rel_retval=NULL, **GCM_rel_reservoir=NULL;
int GCM_rel_reservoir_index = 0, GCM_rel_reservoir_count=0, GCM_rel_nodrifts=0;

#pragma omp threadprivate(GCM_rel_particle,GCM_rel_retval,GCM_rel_reservoir,GCM_rel_reservoir_index,GCM_rel_reservoir_count,GCM_rel_nodrifts)

void equation_GCM_rel_init(settings *set) {
	/* Pre-allocate vectors */
	GCM_rel_retval = vnew(5);

	GCM_rel_nodrifts = set->nodrifts;
}

vector *equation_GCM_rel_vnew(void) {
	if (GCM_rel_reservoir_index+1 >= GCM_rel_reservoir_count) {
		GCM_rel_reservoir = realloc(GCM_rel_reservoir, sizeof(vector*)*(GCM_rel_reservoir_count+1));

		GCM_rel_reservoir[GCM_rel_reservoir_count] = vnew(3);
		GCM_rel_reservoir_count++;
	}

	return GCM_rel_reservoir[GCM_rel_reservoir_index++];
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
void equation_GCM_rel_init_run(particle *part, ode_solution *solver_object) {
	GCM_rel_particle = part;
	GCM_rel_reservoir_index = 0;

	double m = part->mass; // particle mass
	double q = part->charge; // particle charge
	double c = 299792458.0;	// speed of light

	/* initial particle position */
	double x = part->r0[0],
	y = part->r0[1],
	z = part->r0[2];
	/* Create initial particle position vector, r */

	/* initial particle momentum */
	double vx = part->v0[0],
	vy = part->v0[1],
	vz = part->v0[2];
	double v2 = vx*vx + vy*vy + vz*vz;

	/* Calculate gamma */
	double gamma = 1/sqrt(1-v2/(c*c));

	/* Calculate momentum */
	double px = vx*gamma*m,
		   py = vy*gamma*m,
		   pz = vz*gamma*m;

	/* Create initial 3-momentum vector  */
	//vector *p=vinit(3,px,py,pz);
	vector *p = equation_GCM_rel_vnew();
	p->val[0] = px; p->val[1] = py; p->val[2] = pz;

	/* Get magnetic field B for initial position */
	vector *B = magnetic_field_get(x,y,z);

	double Bx=B->val[0],
	By=B->val[1],
	Bz=B->val[2];
	/* Square of B */
	double B2 = Bx*Bx + By*By + Bz*Bz;
	/* Absolute value of B */
	double B_abs=sqrt(B2);

	/* Calculate bhat, unit vector in B-field direction */
	vector *bhat=vmulsf(1/B_abs,B);
	double bx=bhat->val[0],
	by=bhat->val[1],
	bz=bhat->val[2];

	/* Compute magnitude (with sign) of parallel momentum */
	double ppar_abs=vdot3(bhat,p);

	/* negative parallel velocity vector, to calculate
	* perpendicular momentum vector */
	vector *bhat_ppar = vmuls(-ppar_abs,bhat,equation_GCM_rel_vnew());
	/* perpendicular velocity */
	vector *pperp= vadd(p,bhat_ppar,equation_GCM_rel_vnew());
	/* Get components of p_perp */
	double pperp_x=pperp->val[0], pperp_y=pperp->val[1], pperp_z=pperp->val[2];

	/* Allocate and calculate Larmor vector rho = pperp/qB */
	vector *rho=equation_GCM_rel_vnew();
	rho->val[0]=by*pperp_z-bz*pperp_y;
	rho->val[1]=bz*pperp_x-bx*pperp_z;
	rho->val[2]=bx*pperp_y-by*pperp_x;
	vmulsf(-1/(q*B_abs), rho);

	/* Calculate magnetic moment mu */
	double mu=vdot3(pperp,pperp)/(2*m*B_abs);

	if (part->gc_position) {
		solver_object->Z->val[0] = x;
		solver_object->Z->val[1] = y;
		solver_object->Z->val[2] = z;
		solver_object->Z->val[3] = ppar_abs;
		solver_object->Z->val[4] = mu;
	} else {
		/* Calculate guiding center position vector X */
		double X=x+rho->val[0], Y=y+rho->val[1], Z=z+rho->val[2];

		/* Store initial values to solution vector */
		solver_object->Z->val[0] = X; // Guiding center x-position
		solver_object->Z->val[1] = Y; // Guiding center y-position
		solver_object->Z->val[2] = Z; // Guiding center z-position
		solver_object->Z->val[3] = ppar_abs;  // Parallel momentum
		solver_object->Z->val[4] = mu;        // Magnetic moment mu
	}

	/* Save initial data in solver object */
	solver_object->step = 1e-8; /* Initial step size */
	solver_object->flag=0;
}

/**
 * Initialize output data structure
 */
solution_data *equation_GCM_rel_init_output(solution_data *output) {
	output->nvars = 5;
	output->labels = malloc(sizeof(char*)*output->nvars);
	output->labels[0] = "X";
	output->labels[1] = "Y";
	output->labels[2] = "Z";
	output->labels[3] = "ppar";
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
vector *equation_GCM_rel_eq(double T, vector *Z) {
	particle *part = GCM_rel_particle;
	double m = part->mass;
	double q = part->charge;
	double c = 299792458;	// speed of light

	double x=Z->val[0], y=Z->val[1], z=Z->val[2],
		   Xdot1, Xdot2, Xdot3, ppardot, mu;

	mu = Z->val[4];

	/* Get magnetic field data (including gradB, curlB and Babs) */
	diff_data *dd = magnetic_field_diff(x,y,z);
	/* Calculate gamma */
	double pperp2 = 2*dd->Babs*m*mu;
	double ppar2 = Z->val[3]*Z->val[3];
	//printf("ppar = %e\n", sqrt(ppar2));
	double p2 = ppar2 + pperp2;
	double gamma = sqrt(1+p2/(c*c*m*m));
	//printf("cos(pitch) = %e\n", sqrt(ppar2/p2));

	/***************************
	* Calculate ppardot        *
	***************************/
	/* Calculate bhat */
	double bhatval[3];
	vector bhatv;
	bhatv.n = 3;
	bhatv.val = bhatval;
	vector *bhat = vmuls(1/dd->Babs, dd->B, &bhatv);

	/* Calculate B* (effective B-field) */
    // dd->curlB is destroyed here!
	//vector *B_eff = vaddf(dd->B, vmulsf(m/q*Z->val[3], dd->curlB));
	vector *B_eff = vaddf(dd->B, vmulsf(Z->val[3]/q, dd->curlB));
    /*
    vector *curlBhat = dd->curlB;
    curlBhat->val[0] = (curlBhat->val[0] + bhat->val[1]*dd->gradB->val[2] - bhat->val[2]*dd->gradB->val[1]) / dd->Babs;
    curlBhat->val[1] = (curlBhat->val[1] + bhat->val[2]*dd->gradB->val[0] - bhat->val[0]*dd->gradB->val[2]) / dd->Babs;
    curlBhat->val[2] = (curlBhat->val[2] + bhat->val[0]*dd->gradB->val[1] - bhat->val[1]*dd->gradB->val[0]) / dd->Babs;
    vector *B_eff = vaddf(dd->B, vmulsf(Z->val[3]/q, curlBhat));
    */

	/* Get value of B* parallel to b^ */
	double Beff_par = vdot3(B_eff, bhat);

	/* == > */
	if (GCM_rel_nodrifts) {
		Xdot1 = Z->val[3] / (gamma*m) * bhat->val[0];
		Xdot2 = Z->val[3] / (gamma*m) * bhat->val[1];
		Xdot3 = Z->val[3] / (gamma*m) * bhat->val[2];

        ppardot = -mu/gamma * vdot3(dd->gradB, bhat);
	} else {
		/***************************
		* Calculate Xdot1         *
		***************************/
		Xdot1 = 1/(gamma*Beff_par) * (-mu/q*(dd->gradB->val[1]*bhat->val[2] - dd->gradB->val[2]*bhat->val[1]) + Z->val[3]/m * B_eff->val[0]);
		/***************************
		* Calculate Xdot2         *
		***************************/
		Xdot2 = 1/(gamma*Beff_par) * (-mu/q*(dd->gradB->val[2]*bhat->val[0] - dd->gradB->val[0]*bhat->val[2]) + Z->val[3]/m * B_eff->val[1]);
		/***************************
		* Calculate Xdot3         *
		***************************/
		Xdot3 = 1/(gamma*Beff_par) * (-mu/q*(dd->gradB->val[0]*bhat->val[1] - dd->gradB->val[1]*bhat->val[0]) + Z->val[3]/m * B_eff->val[2]);

        ppardot = -mu/(gamma*Beff_par) * vdot3(dd->gradB, B_eff);
	}

	vector *retval = GCM_rel_retval;
	retval->val[0] = Xdot1;
	retval->val[1] = Xdot2;
	retval->val[2] = Xdot3;
	retval->val[3] = ppardot;
	retval->val[4] = 0;

	double ptov = 1/(m*gamma);
	double vpar = Z->val[3]*ptov;
	double vperp = sqrt(pperp2)*ptov;

	double Xmag = hypot(Xdot1, hypot(Xdot2, Xdot3));
	tool_record_data(
		T, ppar2, pperp2,		/* Time and momentum squared */
		dd->Babs,				/* Magnetic field strength */
		x, y, z,				/* Guiding-center position */
		Xdot1/Xmag, Xdot2/Xmag, Xdot3/Xmag,	/* vhat, direction of guiding-center velocity */
		vpar, vperp				/* Particle speed */
	);

	return retval;
}
