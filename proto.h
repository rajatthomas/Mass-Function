/*! FILE : proto.h

Contains all the function definitions required
for the CODE

Author:
Rajat Mani Thomas
IPMU, Japan
CITA, Toronto
Enquiries: thomas@cita.utoronto.ca
Private: rajatthomas@gmail.com

*/

/*---------------*/
#ifndef PROTO_H
#define PROTO_H
/*---------------*/

#ifndef ALLVARS_H
#include "allvars.h"
#endif

/* Convenience from Numerical Recipes in C, 2nd edition */
#define SQR(a) (((a)) == 0.0 ? 0.0 : (a)*(a))

#define FAIL	0
#define SUCCESS	1

/***Gravitational Constant in terms of (km/s)^2 Mpc M_sun^ -1 ***/
#define GravConst_km_sec_solarmass 4.299E-9


int readparams(char *);

int TransferFunc_Init( );
double TransferFuncCDM_BBKS(double );
double TransferFuncCDM_EH(double );

double MassFunc(double ,double );
double Hubble(double );
double omega(double );
double Delta_c(double );
double Delta_vir(double );


double GrowthFactor(double );
double GrowthFactor_aux1(double );
double GrowthFactor_aux2(double );
double GrowthFactor_aux3(double );
double GrowthFactor_aux4(double , void *);

/* Functions used to calculate the Sigma(mass) */
double Sigma(double );
double sig2tophat(double , void *);
double sig2gauss(double , void *);
double sig2sharpk(double , void *);
double sig2tophatnu(double, void *);
double sig2gaussnu(double , void *);
double sig2sharpknu(double , void *);

/* The moment function in moments.c */
double momentfunc(double , void * );
double zeta(double , double );
/*----------------*/
#endif

