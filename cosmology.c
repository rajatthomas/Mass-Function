/*<! File: cosmology.c

Contains expressions for Transfer functions, window functions
and massfunctions and so on.

Written by:
 Rajat Mani Thomas
 CITA, U of T

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "proto.h"

/*-----------------------------------------------

Initialize the cosmology, required for the 
Transfer functions. The input parameters are taken
from the parameters file.

------------------------------------------------*/

/***** GLOBAL VARIABLES for transfer functions ---*/

double	omhh,		/* Omega_matter * h_z0^2 */
        alpha_gamma,	/* sqrt(alpha_nu) */
	sound_horizon_fit,  /* The sound horizon at the drag epoch */
	p_cb,		/* The correction to the exponent after drag epoch */
	theta_cmb,	/* The temperature of the CMB, in units of 2.7 K */
	f_hdm,		/* Massive Neutrino fraction */
	f_cb,		/* Baryon + CDM fraction */
	num_degen_Neutrino,	/* Number of degenerate massive neutrino species */
	growth_k0;	/* D_1(z) -- the growth function as k->0 */

int TransferFunc_Init(){

double	f_baryon,	/* Baryon fraction */
	f_bnu,		/* Baryon + Massive Neutrino fraction */
	f_cdm,		/* CDM fraction */
	hh_z0,	        /* Need to pass Hubble constant to TFmdm_onek_hmpc() */
	k_equality,	/* The comoving wave number of the horizon at equality*/
	obhh,		/* Omega_baryon * h_z0^2 */
	omega_curv,	/* = 1 - Omega_Matter - Omega_Lambda */
	Omega_Lambda_z, /* Omega_lambda at the given redshift */
	Omega_Matter_z,	/* Omega_matter at the given redshift */
	onhh,		/* Omega_hdm * h_z0^2 */
	p_c,		/* The correction to the exponent before drag epoch */
	y_drag,		/* Ratio of z_equality to z_drag */
	z_drag;		/* Redshift of the drag epoch */


    double z_drag_b1, z_drag_b2, omega_denom;
    int qwarn;
    qwarn = 0;

    theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */

    /* Look for strange input */
    if (Omega_Baryon<0.0) {
	fprintf(stderr,
	  "TransferFunc_Init(): Negative Omega_Baryon set to trace amount.\n");
	qwarn = 1;
    }
    if (Omega_Neutrino<0.0) {
	fprintf(stderr,
	  "TransferFunc_Init(): Negative Omega_Neutrino set to trace amount.\n");
	qwarn = 1;
    }
    if (h_z0<=0.0) {
	fprintf(stderr,"TransferFunc_Init(): Negative Hubble constant illegal.\n");
	exit(1);  /* Can't recover */
    } else if (h_z0>2.0) {
	fprintf(stderr,"TransferFunc_Init(): Hubble constant should be in units of 100 km/s/Mpc.\n");
	qwarn = 1;
    }
    if (redshift<=-1.0) {
	fprintf(stderr,"TransferFunc_Init(): Redshift < -1 is illegal.\n");
	exit(1);
    } else if (redshift>99.0) {
	fprintf(stderr,
	  "TransferFunc_Init(): Large redshift entered.  TF may be inaccurate.\n");
	qwarn = 1;
    }
    if (degen_Neutrino<1) degen_Neutrino=1;
    num_degen_Neutrino = (double) degen_Neutrino;	
	/* Have to save this for TFmdm_onek_mpc() */
    /* This routine would crash if baryons or neutrinos were zero, 
	so don't allow that */
    if (Omega_Baryon<=0) Omega_Baryon=1e-5;
    if (Omega_Neutrino<=0) Omega_Neutrino=1e-5;

    omega_curv = 1.0-Omega_Matter-Omega_Lambda;
    omhh = Omega_Matter*SQR(h_z0);
    obhh = Omega_Baryon*SQR(h_z0);
    onhh = Omega_Neutrino*SQR(h_z0);
    f_baryon = Omega_Baryon/Omega_Matter;
    f_hdm = Omega_Neutrino/Omega_Matter;
    f_cdm = 1.0-f_baryon-f_hdm;
    f_cb = f_cdm+f_baryon;
    f_bnu = f_baryon+f_hdm;

    /* Compute the equality scale. */
    z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));	/* Actually 1+z_eq */
    k_equality = 0.0746*omhh/SQR(theta_cmb);

    /* Compute the drag epoch and sound horizon */
    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*
		(1.0+z_drag_b1*pow(obhh,z_drag_b2));
    y_drag = z_equality/(1.0+z_drag);

    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));

    /* Set up for the free-streaming & infall growth function */
    p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
    p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

    omega_denom = Omega_Lambda+SQR(1.0+redshift)*(omega_curv+
			Omega_Matter*(1.0+redshift));
    Omega_Lambda_z = Omega_Lambda/omega_denom;
    Omega_Matter_z = Omega_Matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
    growth_k0 = z_equality/(1.0+redshift)*2.5*Omega_Matter_z/
	    (pow(Omega_Matter_z,4.0/7.0)-Omega_Lambda_z+
	    (1.0+Omega_Matter_z/2.0)*(1.0+Omega_Lambda_z/70.0));
    growth_to_z0 = z_equality*2.5*Omega_Matter/(pow(Omega_Matter,4.0/7.0)
	    -Omega_Lambda + (1.0+Omega_Matter/2.0)*(1.0+Omega_Lambda/70.0));
    growth_to_z0 = growth_k0/growth_to_z0;	
    
    /* Compute small-scale suppression */
    alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*
	pow(1+y_drag,p_cb-p_c)*
	(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/
	(1-0.193*sqrt(f_hdm*num_degen_Neutrino)+0.169*f_hdm*pow(num_degen_Neutrino,0.2))*
	(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    alpha_gamma = sqrt(alpha_nu);
    beta_c = 1/(1-0.949*f_bnu);
    /* Done setting scalar variables */
    hh_z0 = h_z0;	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
/* Critical density of the universe */
rho_crit_0 = 3 * SQR(Hubble_z0)/(8.0 * M_PI * GravConst_km_sec_solarmass);
    return qwarn;
}

double TransferFuncCDM_BBKS(double k){

/* Verify using Eq. 15.82 in Peacock's book */

  double q;
  double t;

/* Gamma is the shape parameter Eq. 15.85 Peocock.*/
double Gamma = Omega_Baryon * h_z0;

  q = k/(Gamma);
  t = (log(1+2.34*q)/(2.34*q))*pow(
                    1.+3.89*q+(16.1*q)*(16.1*q)+
                    (5.46*q)*(5.46*q)*(5.46*q)+
                    (6.71*q)*(6.71*q)*(6.71*q)*(6.71*q),-.25);
  
  return t; 
}

double TransferFuncCDM_EH(double kk){

/* Eisenstein & Hu (1997) Better approximation when baryon */
/* content of the Universe is large */

/* Finally, the function gives its  answers as */
double   tf_cb,		/* The transfer function for density-weighted
			CDM + Baryon perturbations. */
	tf_cbnu;	/* The transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */

/* By default, these functions return tf_cb */

double	gamma_eff,	/* Effective \Gamma */
	growth_cb,	/* Growth factor for CDM+Baryon perturbations */
	growth_cbnu,	/* Growth factor for CDM+Baryon+Neutrino pert. */
	max_fs_correction,  /* Correction near maximal free streaming */
	qq,		/* Wavenumber rescaled by \Gamma */
	qq_eff,		/* Wavenumber rescaled by effective Gamma */
	qq_nu,		/* Wavenumber compared to maximal free streaming */
	tf_master,	/* Master TF */
	tf_sup,		/* Suppressed TF */
	y_freestream; 	/* The epoch of free-streaming for a given scale */

/* Given a wavenumber in Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
/* Input: kk -- Wavenumber in Mpc^-1 */
/* Output: The following are set as global variables:
	growth_cb -- the transfer function for density-weighted
			CDM + Baryon perturbations. 
 	growth_cbnu -- the transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */
/* The function returns growth_cb */

    double tf_sup_L, tf_sup_C;
    double temp1, temp2;

    qq = kk/omhh*SQR(theta_cmb);

    /* Compute the scale-dependent growth functions */
    y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*
		SQR(num_degen_Neutrino*qq/f_hdm);
    temp1 = pow(growth_k0, 1.0-p_cb);
    temp2 = pow(growth_k0/(1+y_freestream),0.7);
    growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
    growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;

    /* Compute the master function */
    gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/
		(1+SQR(SQR(kk*sound_horizon_fit*0.43))));
    qq_eff = qq*omhh/gamma_eff;

    tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
    tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
    tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));

    qq_nu = 3.92*qq*sqrt(num_degen_Neutrino/f_hdm);
    max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_Neutrino,0.3+0.6*f_hdm)/
		(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
    tf_master = tf_sup*max_fs_correction;

    /* Now compute the CDM+HDM+baryon transfer functions */
    tf_cb = tf_master*growth_cb/growth_k0;
    tf_cbnu = tf_master*growth_cbnu/growth_k0;
    return tf_cb;
}

double MassFunc(double Mass, double redshift){


double rho_average; /* Average density in the Universe at z */
double nu;          /* delta/sigma(M) */
double dlnS_dlnM;   /* dSigma(M)/Sigma(M) * M/dM */
double M1,M2,eps =.001; /* to calculate the interval dM around M */
double massfcn=0.0;

rho_average = rho_crit_0 * Omega_Matter;
nu = Delta_c(redshift) / Sigma(Mass) ;
//nu = 1.686 / Sigma(Mass) ;

//printf("rho_crit = %le \n", rho_crit_0);

/* Now we will compute the approximation to the differential */
/* dlnS/dlnM using a small dM around M */
M1 = (1.0 - eps) * Mass;
M2 = (1.0 + eps) * Mass;

dlnS_dlnM = ((Sigma(M2) - Sigma(M1))/Sigma(Mass)) * (Mass / (M2 - M1));

//printf("dlnS_dlnM = %f \n",dlnS_dlnM);
//printf("rho average = %le \n",rho_average);

/* Return the final form of the Mass functions */
  if(strcmp(MassFuncType,"PS")==0)
	massfcn =  (sqrt(2/M_PI) * (rho_average/SQR(Mass)) * nu * exp(-0.5*SQR(nu)) 
                * fabs(dlnS_dlnM));

  if(strcmp(MassFuncType,"ST")==0){
            nu = 0.84083292 * nu;
        massfcn = (sqrt(2/M_PI) * (rho_average/SQR(Mass)) * 0.3222 * nu * 
               exp(-0.5*SQR(nu)) * (1.0 + 1.0/pow(nu,0.6)) * fabs(dlnS_dlnM)); 
                             }

//printf("massfcn = %le \n",massfcn); 
//exit(0);
return massfcn;
}

/**********************************************************************

      FUNCTION Hubble(z)
----------------------------------------------------------------------
c Calculate hubble constant (in physical units; km/s/Mpc) at redshift z.
c
c----------------------------------------------------------------------*/

      
double Hubble(double z)
{
  double  fac;
  fac = Omega_Lambda + (1.0 - Omega_Lambda - Omega_Matter) * SQR(1+z) + 
    Omega_Matter * pow((1+z),3.0);
  return Hubble_z0 * sqrt(fac);	
}


/**********************************************************************

      FUNCTION omega(z)
c----------------------------------------------------------------------
c
c Calculate the density parameter omega at redshift z. Note that omega 
c is only the mater contribution.
c
c----------------------------------------------------------------------*/

double omega(double z)
{
  return (Omega_Matter * pow((1.0 + z),3.0) / pow((Hubble(z)/Hubble_z0),2.0));
}



/**********************************************************************

      FUNCTION Delta_c(z)
c----------------------------------------------------------------------
c
c Calculate the critical overdensity.
c We use the approximation in NFW97
c
c--------------------------------------------------------------------*/

double Delta_c(double z)
{

  double     ddd,dc0,omz;

  omz = omega(z);
  ddd = 1.0 / GrowthFactor(z);
  dc0 = 0.0;
 
     if(fabs(1.0-Omega_Matter-Omega_Lambda)<1.0e-4)
	 dc0 = 0.15 * pow((12.0 * M_PI),(2.0/3.0)) * pow(omz,0.0055);
     if((Omega_Matter<1.0) && (Omega_Lambda==0.0))
	 dc0 = 0.15 * pow((12.0 * M_PI),(2.0/3.0)) * pow(omz,0.0185);
     
     if(dc0==0.0)
	      {
	  printf(" Delta_c not defined for this cosmology\n");
	  printf(" Omega_Matter = %lf, Omega_Lambda: %lf \n",Omega_Matter,Omega_Lambda);
	  printf(" omega_z = %lf, GrowthFactor: %lf \n",omz,ddd);
	  printf(" at redshift: %lf \n",z);
	      }
 
//printf("dc0 = %f \t ddd= %f \n",dc0,ddd);
//return dc0 * ddd;
return dc0;
}

/**********************************************************************

     FUNCTION Delta_vir(z)
c----------------------------------------------------------------------
c
c Calculate the virial density in terms of critical density of the 
c Universe. We use the fitting formulae of Bryan & Norman, 1998, ApJ,
c 495, 80. These are accurate to better than 1% for 0.1 < Omega_Matter < 1.0
c
c---------------------------------------------------------------------*/

double Delta_vir(double z)
{

  double   x,Delta_critical;

  x = omega(z) - 1.0;

      if(fabs(1.0-Omega_Matter-Omega_Lambda)<1.0e-4)
	Delta_critical = 18.0 * M_PI*M_PI + 82.0*x - 39.0 * x*x;
      else	      
         if((Omega_Matter<1.0) && (Omega_Lambda==0.0))
	   Delta_critical = 18.0 * M_PI*M_PI + 60.0*x - 32.0 * x*x;
         else
               {
	   printf("Delta_crit not defined for this cosmology\n");
     		exit(0);
               }

      return Delta_critical;
}


/*--------------------------------------------------------------------
c
c The linear growth factor at redshift z (see NFW97)
c
c--------------------------------------------------------------------*/

double GrowthFactor(double z)
{

  double    ddtemp;
  double   ww,w_0,y_0,y;
  
   if((Omega_Matter==1.0)&&(Omega_Lambda==0.0))
     {
       ddtemp = 1.0/(1.0 + z);
       return ddtemp;
     }
      else if((Omega_Matter<1.0)&&(Omega_Lambda==0.0))
          	{
        	w_0 = 1.0/Omega_Matter - 1.0;
       		ww = w_0 / (1.0 + z);
       		ddtemp = (double)(GrowthFactor_aux1(ww)/GrowthFactor_aux1(w_0));
       		return ddtemp;
      
     		}	
    		else if(fabs(1.0-Omega_Matter-Omega_Lambda)<1.0e-4)
      		{
        	w_0 = 1.0/(double)(Omega_Matter) - 1.0;
        	y_0 = pow((2.0 * w_0),(1.0/3.0));
        	y = y_0/(1.0 + (double)(z));  
        	ddtemp = (double)((GrowthFactor_aux2(y) * GrowthFactor_aux3(y))/(GrowthFactor_aux2(y_0) * GrowthFactor_aux3(y_0)));
       		return ddtemp;
      		}
			else {
            		printf("Not the right cosmology: error in GrowthFactor \n");
                        exit(0);
 
         		}
}

/**********************************************************************

       FUNCTION GrowthFactor_aux1(x)
c--------------------------------------------------------------------
c
c Auxialiary function used in computation of growth rate
c
c--------------------------------------------------------------------*/

double GrowthFactor_aux1(double x)
{
  double   fac1,fac2;

  fac1 = sqrt(1.0+x) - sqrt(x);
  fac2 = 3.0*sqrt(1.0+x)/pow(x,1.5);

 return (1.0 + (3.0/x) + fac2 * log(fac1));
 
}

/**********************************************************************

      FUNCTION GrowthFactor_aux2(x)
c--------------------------------------------------------------------
c
c Auxialiary function used in computation of growth rate
c
c--------------------------------------------------------------------*/
       
double GrowthFactor_aux2(double x)
{
return sqrt(pow(x,3.0) + 2.0)/pow(x,1.5);
}

/**********************************************************************

      FUNCTION GrowthFactor_aux3(x)
c--------------------------------------------------------------------
c
c Auxialiary function used in computation of growth rate
c
c--------------------------------------------------------------------*/

double GrowthFactor_aux3(double x){


  double  SS;
  double epsabs=0.0, abserr, epsrel=1.4899999999999999e-08;

gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);


         gsl_function F;
         F.function = &GrowthFactor_aux4;

gsl_integration_qags(&F,0.0,x,epsabs,epsrel,1000, w, &SS,&abserr);
         
gsl_integration_workspace_free (w);

  return (double)(SS);

  }

/**********************************************************************
   
      FUNCTION GrowthFactor_aux4(x)
c--------------------------------------------------------------------
c
c Auxialiary function used in computation of growth rate
c
c--------------------------------------------------------------------*/

double GrowthFactor_aux4(double x, void * params)
{
return pow((x/(pow(x,3.0) + 2.0)),1.5);
}


