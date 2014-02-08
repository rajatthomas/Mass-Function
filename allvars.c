/*! \file allvars.c
 *  \brief provides instances of all global variables.
 */

#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"

/* Cosmolgoy parameters */
double Omega_Lambda,Omega_Baryon,Omega_Matter,Omega_Neutrino;
double tilt, Sigma_8, Hubble_z0, h_z0, redshift;
double rho_crit_0;
int degen_Neutrino;
double alpha_nu,       /* The small-scale suppression */
       beta_c,         /* The correction to the log in the small-scale */
       growth_to_z0,   /* D_1(z)/D_1(0) -- the growth relative to z=0 */
       z_equality;     /*Redshift of matter-radiation equality */

/* Transfer & Mass function  */
char SigmaFilter[32], SigmaMethod[32], MassFuncType[32], TranFuncType[32];
int SigmaFirstTime;

gsl_rng *random_generator; /*!< the employed random number generator of the GSL li
brary */

/* Redshifts spanned */
double z_start,z_end,z_stride;
double mass_min,mass_max;
int mass_bins;

/* Switches */

int CalcMoment, CalcMassFunc;
