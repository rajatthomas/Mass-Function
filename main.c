/*<! File: main.c

The main program for the Press-Shechter code.

Written by:
 Rajat Mani Thomas
 CITA, Toronto
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "proto.h"

int main (int argc, char** argv){

/**********************************************************************/
/* Read the parameter file which contains the  cosmology given and */
/* the type of massfunction and transferfunctions that are used    */
/**********************************************************************/
  char ParameterFile[64],filename[32];
  FILE *fp;
  double mi,zi,dm,tmp;
  double massfcn;


/***Variables for the integration of the moment function *****/
  double epsabs=0.0, abserr, epsrel=1.e-7;
  double moment_result, Mmin;

  if(argc != 2) {
    printf(" Error in usage of code... format: ./ps parameters.txt .. \n");
    printf(" Exiting \n"); exit(0);
  }

  sprintf(ParameterFile,"%s",argv[1]);

  readparams(ParameterFile); /* Defined in readparams.c */

/*********************************************************************/

  if(strcmp(MassFuncType,"PS")==0)
  printf("Generating the Press-Schechter Mass Function ... \n ");
  else if(strcmp(MassFuncType,"ST")==0)
  printf("Generating the Sheth-Tormen Mass Function ... \n ");
       else { printf(" Unknown mass function requested..exiting \n");
              exit(0);
            }

  if(strcmp(TranFuncType,"BBKS")==0)
  printf("Using BBKS type transfer Function ... \n ");
  else if(strcmp(TranFuncType,"EH")==0)
  printf("Using Eisenstein & Hu Transfer Function ... \n ");
       else { printf(" Unknown transfer function requested..exiting \n");
              exit(0);
            }

  printf(" Input parameters used are in %s \n", ParameterFile);
  printf("(Omg_b,%f) \t (Omg_m,%f) \t (Omg_L,%f) \t (Sigma_8,%f) \t (Hubble(0),%f) \n",Omega_Baryon,Omega_Matter,Omega_Lambda,Sigma_8,Hubble_z0);


/************************************************/
/**** Find the massfunction for each redshift ***/
/************************************************/

if(z_start < z_end) {
tmp = z_start;
z_start = z_end;
z_end = tmp;      /* Swap the order of redshifts if necessary */ 
}
    for(zi = z_start; zi > z_end ; zi -= z_stride){

redshift = zi;
  printf(" z = %f \n", redshift);
  TransferFunc_Init();  /* Initialize for given redshift */

  dm = (mass_max - mass_min)/mass_bins;

SigmaFirstTime = 1; /* First time always use the numerical */
                    /* method of calculating Sigma(M)      */

if(CalcMassFunc){

/****Generate a filename ********/
sprintf(filename,"MF_%s_%s_z%2.1f.dat",MassFuncType,TranFuncType,zi);
fp = fopen(filename,"w");
/********************************/

/******Evaluate the number density across the mass range*****/
  for(mi = mass_min; mi < mass_max; mi += dm){
   massfcn = MassFunc(pow(10.0,mi),zi);  
  fprintf(fp,"%le \t %le \t %le\n", pow(10.0,mi), massfcn, massfcn*SQR(pow(10.0,mi))/(rho_crit_0 * Omega_Matter));
                               
} /* End of mi loop */

fclose(fp);
}
/* Calculate the moments of the function if asked */

if(CalcMoment) /* as given in moments.c program */
{

Mmin = 1e4; /* Minimum mass from where to integrate */

gsl_integration_workspace *w = gsl_integration_workspace_alloc (10000);


         gsl_function F;
         F.function = &momentfunc;

/* integral from Mmin-to-infinity */
	 gsl_integration_qags (&F, log(Mmin), log(1e17), epsabs,epsrel, 10000, w, &moment_result, &abserr);

gsl_integration_workspace_free (w);

printf("The integral = %le at z = %f\n", moment_result, redshift);

}


} /* End of zi loop */

  return 0;
}
