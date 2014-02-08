/* Code to calculate sigma8 for COBE-normalized CDM and MDM models, 
using the approximation of Eisenstein & Hu (1997). */

/* This code will calculate sigma_8 and other quantities for CDM and
MDM models.  Because the fitting form for the transfer function does not
include the baryon oscillations, the answer will not be accurate for
large baryon fractions.  The COBE normalization is from the papers
on Bunn & White (1997) and Hu & White (1997). */

/* Accuracy is generally better than 5%.  Recall that the COBE normalization
is only 7% accurate (1-sigma).  This program gives low answers compared to
CMBfast for lower values of Omega_0*h^2, and vice versa for higher values. 
Values of Omega_0 near 0.2 can produce answers that are 5% low. */

/* The integrator used here is a bit crude, but suffices for rapid
10^-3 accuracy.  You could get better performance by substituting an
integrator from a numerical package.  See the end of the file for details. */

/* This program uses the routines TransferFunc_Init() and TransferFuncCDM_EH(),
available from http://www.sns.ias.edu/~whu/power/power.html in the file
power.c.  To compile this program, you'll need that file as well, e.g.
	cc -O sigma8.c power.c -o sigma8 -lm 		*/

/* Lists of versions and bugs:
v1.0 -- 10/23/97 -- Original release version.
v1.1 -- 11/13/97 -- We correct a missing factor of Omega_Matter in the
		    conversion between mass and length scales [lines
		    32 and 34 of main()].  This only affects the results 
		    if you explicitly used the -mass option.
v1.2 -- 04/08/98 -- We correct a typo in the second-order tilt scaling
		    of the COBE normalization for flat models with tensors.
v1.3 -- 09/15/98 -- We correct a typo in the tilt scaling of the COBE
		    normalization for open models with tensors.
v1.4 -- 12/02/98 -- We correct the COBE normalization for tilted open models,
		    undoing the change of v1.3 and fixing what was
		    wrong there to begin with.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#define YES_NR
#define WEIRD  -2538	/* A strange number, unlikely to be entered */
#define RHOCRIT 2.78e11	/* The critical density, for h=1, in Msun/Mpc^3 */


/* Global variable */
double   scale;
double   deltaH;		/* COBE-normalization */
double   sig2_scale;	/* For passing a length scale to the integrator */
int 	qtensors, quiet, qredshift;


void crash(char *s)
{
    fprintf(stderr,"%s\n", s); exit(1); return;
}

void warn(char *s)
{
    fprintf(stderr,"%s\n", s); return;
}



int checking()
/* Look for illegal entries */
{

    scale = -1.0;
    qtensors = 0;
    quiet = 0;
    qredshift = 0;

    /* Look for illegal entries */
    if (Omega_Matter<=0.0) crash("Omega must be positive.");
    if (h_z0<=0.0) crash("Hubble constant must be positive.");
	else if (h_z0>10.0) warn("Did you get the Hubble units right?");
    if (degen_Neutrino<0) {
	warn("Negative number of degenerate species reset to 1.");
	degen_Neutrino = 1.0;
    }
    if (redshift<=-1.0) crash("Redshift must exceed -1.");

    /* Now look for fractions, if they were input. */
    if (Omega_Baryon<0.0) Omega_Baryon = Omega_Matter*(-Omega_Baryon);	
    if (Omega_Neutrino<0.0) Omega_Neutrino = Omega_Matter*(-Omega_Neutrino);	
    if (Omega_Lambda== WEIRD) Omega_Lambda = 1-Omega_Matter;	/* Set to be flat */

    if (Omega_Baryon+Omega_Neutrino>Omega_Matter) crash("Baryons + Neutrinos exceeds Omega_0.");
	else if (Omega_Baryon+Omega_Neutrino>0.6*Omega_Matter) 
	    warn("CDM fraction is below range of TF fitting formula.");
    if (tilt>1.3 || tilt < 0.7) 
	warn("Tilt out of range of Bunn and White fitting formula.");
    if (tilt!=1 && (Omega_Matter<0.2 || Omega_Matter>1.0)) 
	warn("Omega out of range of Bunn and White fitting formula.");
    if (tilt==1 && (Omega_Matter<0.2 || Omega_Matter>1.6 || Omega_Lambda<0.0 || Omega_Lambda>0.8))
	warn("Omega or Lambda out of range of Bunn and White fitting formula.");
    if (qtensors!=0&&(tilt>1.0||
	(fabs(Omega_Lambda)>1e-5&&fabs(1.0-Omega_Matter-Omega_Lambda)>1e-5))) 
	crash("Tensor formulae only hold for tilts < 1 and simple cosmologies.");
    if (Omega_Baryon/Omega_Matter>0.4) warn("Baryon oscillations will be large.  TF fitting form will be inaccurate.");
    if (Omega_Neutrino/Omega_Matter>0.4) warn("Neutrino fraction outside of range of TF fitting form.");
    if (Omega_Matter*h_z0*h_z0>0.4) warn("Omega_0 h^2 outside of range of TF fitting form.");
    if (Omega_Matter*h_z0*h_z0<0.06) warn("Omega_0 h^2 outside of range of TF fitting form.");
    if (redshift>29) warn("Redshift is higher than range of TF fitting form.");

    return SUCCESS;
}

/* ==================================== MAIN() ========================= */

double Sigma(double mass)
{
    double sigma = 0.0;
    double cobenorm();
    double epsabs=0.0, abserr, epsrel=1.489999e-08,sig2;
    int nbins = 2000;
   
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (nbins);
   	gsl_function F;
  
  
    checking(); /* Check if cosmology is sensible */

    /* Set all the background quantities */
  if(TransferFunc_Init(redshift)) printf("Some warnings in Init\n");
  
  deltaH = cobenorm();
quiet=1;
    if (!quiet) {	/* Echo cosmological information to the screen */
	printf("Omega_0  = %5.3f, Omega_baryon = %5.3f, Omega_neutrino = %5.3f\n",
		Omega_Matter, Omega_Baryon, Omega_Neutrino);
	printf("Hubble   = %5.3f, Lambda       = %5.3f, N_neutrinos    = %d\n",
		h_z0, Omega_Lambda, (int) degen_Neutrino);
	printf("Redshift = %5.3f, Growth = %6.4f relative to z=0\n",
		redshift, growth_to_z0);
	printf("Tilt     = %5.3f, COBE normalization: deltaH = %8.3e",
		tilt, deltaH);
	if (tilt!=1.0) 
	    if (qtensors) printf(" with tensors.\n");
	    else printf(" without tensors.\n");
	else printf("\n");
    }
    /* If we've been asked to compute for a particular scale, do that */
    
     
   if (scale>0 || mass>0) {
	if (mass<0)
	    mass = 4.0*M_PI/3.0*Omega_Matter*RHOCRIT*scale*scale*scale;
	else if (scale<0) 
	    scale = pow(3.0/4.0/M_PI/RHOCRIT/Omega_Matter*mass, 1.0/3.0);
	/* Integrate tophat for the given scale */
	sig2_scale = scale;

        /* Choose the right filter */
	
        if(strcmp(SigmaFilter,"top-hat")==0)
         F.function = &sig2tophat;
        
	if(strcmp(SigmaFilter,"gaussian")==0)
	 F.function = &sig2gauss;

        if(strcmp(SigmaFilter,"sharp-k")==0)
         F.function = &sig2sharpk;
	
  	 /* Performing the integral */
 	gsl_integration_qags(&F,0.0,1.0,epsabs,epsrel,nbins,w,&sig2,&abserr); 
 
	sigma = sqrt(sig2);
    }

	gsl_integration_workspace_free (w);

  return sigma;
}

/* ========================== Integration Kernals ======================= */

double W2tophat(double x)
/* The square of the tophat window */
{
    double j1onx;
    if (x<0.03) 
	j1onx = (1.0/3.0-x*x*1.0/30.0);
    else j1onx = (cos(x)-sin(x)/x)/x/x;
    return 9.0*j1onx*j1onx;
}

double sig2tophat(double x, void * params)
/* This is the integrand for the mass fluctuation integral
sig2 = \int_0^\infty (dk/k) Delta^2(k) W^2(kr)
broken and recombined so that one integrates on (0,1] in a new variable x. */
/* Here, the window is a real-space tophat of radius sig2_scale, the
latter being passed from outside in units of h^-1 Mpc */
{
    double powerspcb(double k);
    return (powerspcb(x/sig2_scale)*W2tophat(x)+
 	    powerspcb(1.0/x/sig2_scale)*W2tophat(1.0/x))/x;
}

double sig2tophatnu(double x, void * params)
/* Same as sig2tophat(), but with P_cbnu */
{
    double powerspcbnu(double k);
    return (powerspcbnu(x/sig2_scale)*W2tophat(x)+
 	    powerspcbnu(1.0/x/sig2_scale)*W2tophat(1.0/x))/x;
}

double sig2gauss(double x, void * params)
/* This is the integrand for the mass fluctuation integral
sig2 = \int_0^\infty (dk/k) Delta^2(k) W^2(kr)
broken and recombined so that one integrates on (0,1] in a new variable x. */
/* Here, the window is a gaussian of radius sig2_scale, the
latter being passed from outside in units of h^-1 Mpc */
{
    double powerspcb(double k);
    return (powerspcb(x/sig2_scale)*exp(-x*x)+
 	    powerspcb(1.0/x/sig2_scale)*exp(-1.0/x/x))/x;
}

double sig2gaussnu(double x, void * params)
/* Same as sig2gauss(), but with P_cbnu */
{
    double powerspcbnu(double k);
    return (powerspcbnu(x/sig2_scale)*exp(-x*x)+
 	    powerspcbnu(1.0/x/sig2_scale)*exp(-1.0/x/x))/x;
}

double sig2sharpk(double x, void * params)
/* This is the integrand for the mass fluctuation integral
sig2 = \int_0^\infty (dk/k) Delta^2(k) W^2(kr)
broken and recombined so that one integrates on (0,1] in a new variable x. */
/* Here, the window is a k-space tophat of radius 1/sig2_scale, the
latter being passed from outside in units of h^-1 Mpc */
{
    double powerspcb(double k);
    return (powerspcb(x/sig2_scale))/x;
}

double sig2sharpknu(double x, void * params)
/* Same as sig2sharpk(), but with P_cbnu */
{
    double powerspcbnu(double k);
    return (powerspcbnu(x/sig2_scale))/x;
}

/* ============================ Power Spectrum and COBE ================ */
double powerspcbnu(double k)
/* Returns Delta^2(k), COBE-normalized, based on the approximations of 
Eisenstein & Hu (1997).  k is in h^-1 Mpc. */
/* TransferFunc_Init() and cobenorm() must have been called before this */
/* Returns the density-weighted CDM+Baryon+Neutrino power spectrum */
{
    //extern double tf_cbnu;
    //TransferFuncCDM_EH(k);  /* This sets the value in tf_cbnu */
    return pow(2997.0*k, tilt+3.0)*SQR(deltaH*TransferFuncCDM_EH(k)*growth_to_z0);
}

double powerspcb(double k)
/* Returns Delta^2(k), COBE-normalized, based on the approximations of 
Eisenstein & Hu (1997).  k is in h^-1 Mpc. */
/* TransferFunc_Init() and cobenorm() must have been called before this */
/* Returns the density-weighted CDM+Baryon power spectrum */
{
    return pow(2997.0*k, tilt+3.0)*SQR(deltaH*TransferFuncCDM_EH(k)*growth_to_z0); 
}

double cobenorm()
/* Return the Bunn & White (1997) fit for delta_H */
/* Given Omega_Lambda, Omega_Matter, qtensors, and tilt */
/* Open model with tensors is from Hu & White */
{
    double n;
    n = tilt-1;
    if (fabs(Omega_Matter+Omega_Lambda-1.0)<1e-5) {	/* Flat universe */
	if (qtensors) 
	    return 1.94e-5*pow(Omega_Matter, -0.785-0.05*log(Omega_Matter))*
		exp(n+1.97*n*n);
	else return 1.94e-5*pow(Omega_Matter, -0.785-0.05*log(Omega_Matter))*
		exp(-0.95*n-0.169*n*n);
	/* Error -- was 0.0169 in earlier version of code */
    } else if (fabs(Omega_Lambda)<1e-5) {	/* No Omega_Lambda */
	if (qtensors)
	    return 1.95e-5*pow(Omega_Matter,-0.35-0.19*log(Omega_Matter)-0.15*n)*
		exp(+1.02*n+1.7*n*n);
	 /* Error -- was exp(-1.02n-1.7*n*n) in earlier version (v1-1.3)*/
	else return 1.95e-5*pow(Omega_Matter, -0.35-0.19*log(Omega_Matter)-0.17*n)*
		exp(-n-0.14*n*n);
	 /* Error -- was exp(-n-0.14*n*n) in earlier version (v1-1.2)*/
	 /* Error -- was exp(n+0.14*n*n) in earlier version (v1.3) */
    } else return 1e-5*(2.422-1.166*exp(Omega_Matter)+0.800*exp(Omega_Lambda)
	+3.780*Omega_Matter-2.267*Omega_Matter*exp(Omega_Lambda)+0.487*SQR(Omega_Matter)+
	0.561*Omega_Lambda+3.329*Omega_Lambda*exp(Omega_Matter)-8.568*Omega_Matter*Omega_Lambda+
	1.080*SQR(Omega_Lambda));
}

