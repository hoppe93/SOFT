#ifndef _GLOBAL_H
#define _GLOBAL_H

/* Number of initial simulation points */
#define NUMBER_OF_SIMULATION_POINTS 1000

/* Step in rho of orbit used to calculate spatial jacobian */
#define JACOBIAN_ORBIT_STEP 1e-6

/* conversion from atomic mass units to kg */
#define AMU_TO_KG 1.66053886e-27
/* elementary charge in Coloumbs */
#define CHARGE 1.60217657e-19
/* conversion from joule to eV */
#define ENERGY 6.24150934e18

#endif/*_GLOBAL_H*/
