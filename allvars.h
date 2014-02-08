/*! \file allvars.h
 *  \brief provides instances of all global variables.
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>

/* Cosmolgoy parameters */
extern double Omega_Lambda,Omega_Baryon,Omega_Matter,Omega_Neutrino;
extern double tilt,Sigma_8,Hubble_z0,h_z0, redshift;
extern double rho_crit_0;
extern int degen_Neutrino;
extern double alpha_nu,/* The small-scale suppression */
       beta_c,         /* The correction to the log in the small-scale */
       growth_to_z0,   /* D_1(z)/D_1(0) -- the growth relative to z=0 */
       z_equality;     /*Redshift of matter-radiation equality */


/* Transfer & Mass function  */
extern char SigmaMethod[32], SigmaFilter[32], MassFuncType[32], TranFuncType[32];
extern int SigmaFirstTime;

extern gsl_rng *random_generator; /*!< the employed random number generator of the GSL li
brary */

/* redshifts spanned */
extern double z_start,z_end,z_stride;
extern double mass_min,mass_max;
extern int mass_bins;
#endif

/* Switches */
extern int CalcMoment, CalcMassFunc;

