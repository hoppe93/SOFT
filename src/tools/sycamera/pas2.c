/* Implementation of Eq. (26) in [Pankratov, 1999]
   called P_as2(lambda). Inspired by the Matlab
   implementation done by Adam Stahl.
*/

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "global.h"
#include "magnetic_axis.h"
#include "sycamera.h"

#define IC (1/((LIGHTSPEED*LIGHTSPEED)*mass*mass))

double *sycamera_pas2_lambdas, sycamera_pas2_dl, sycamera_pas2_majorradius;
int sycamera_pas2_lambda_resolution, sycamera_pas2_haswarned;

double pperpPpara, sycamera_gamma, sycamera_gamma2, sycamera_gamma3, eta, etap1, etap12;

#pragma omp threadprivate(sycamera_pas2_haswarned,pperpPpara,sycamera_gamma,sycamera_gamma2,sycamera_gamma3,eta,etap1,etap12)

void sycamera_pas2_init(double lambda0, double lambda1, int res) {
    sycamera_pas2_haswarned = 0;

    sycamera_pas2_lambdas = malloc(sizeof(double)*res);
    sycamera_pas2_lambda_resolution = res;
    int i;
    for (i = 0; i < res; i++) {
        /* We reverse this spectrum for simpler code below */
        sycamera_pas2_lambdas[res-i-1] = lambda0 + i * (lambda1-lambda0) / (res-1);
    }

    /* Compute step length */
    sycamera_pas2_dl = fabs(sycamera_pas2_lambdas[1]-sycamera_pas2_lambdas[0]);

    /* Get device major radius (r of magnetic axis) */
    //sycamera_pas2_majorradius = magnetic_axis_r;
    sycamera_pas2_majorradius = 1.67;
}

int sycamera_pas2_valid(double ppar2, double pperp2, double Bmag, double mass) {
    pperpPpara = sqrt(pperp2/ppar2);
    sycamera_gamma2 = 1 + (ppar2 + pperp2)*IC;
    sycamera_gamma = sqrt(sycamera_gamma2);
    sycamera_gamma3 = sycamera_gamma2 * sycamera_gamma;
    eta = CHARGE * Bmag * sycamera_pas2_majorradius * pperpPpara / (sycamera_gamma * mass * LIGHTSPEED);
    etap1 = 1+eta;
    etap12 = etap1*etap1;

    double q = (4*PI/3) * sycamera_pas2_majorradius*eta / (sycamera_gamma3*etap12*etap1);
    return (sycamera_pas2_lambdas[sycamera_pas2_lambda_resolution-1]*5 < q);
}

double sycamera_pas2_int(void) {
    double prefactor = SYCAMERA_PAS2_CONSTANT * sycamera_gamma * etap12 / (sqrt(eta)*sycamera_pas2_majorradius);
    double exponent = (-4*PI/3) * sycamera_pas2_majorradius / (sycamera_gamma3 * etap1);

    int i;
    double sum=0., ilambda;
    for (i = 0; i < sycamera_pas2_lambda_resolution; i++) {
        ilambda = 1 / sycamera_pas2_lambdas[i];

        sum += prefactor * (ilambda*ilambda) * exp(exponent*ilambda);
    }

    return sum * sycamera_pas2_dl;
}

void sycamera_pas2_test(void) {
    sycamera_pas2_init(5e-7,1e-6,50);
    //double sin2theta = 0.009900990099010;
    //double Bmag = 4.923295e+00;
    double sin2theta = pow(sin(atan(0.1)),2);
    double Bmag = 2.2;

    double p2=1.603286e-20*1.603286e-20;
    sycamera_pas2_valid(p2*(1-sin2theta), p2*sin2theta, Bmag, 9.10938e-31);
    double sum = sycamera_pas2_int();
    printf("p2 = %e   ==>   %e\n", p2, sum);
}
