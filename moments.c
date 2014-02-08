/*<! File: moments.c

The program basically calculates the moments of the
mass function, i.e.,
\int_lm^hm mom_func(m) dn/dm dm

Written by:
 Rajat Mani Thomas
 CITA, Toronto
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

double momentfunc(double mass, void * params){
  mass = exp(mass); /* The log integration trick */
double integrand;

 integrand = zeta(mass,3e7) * pow(mass, .6666); /* Change to whatever moment you want */
 
integrand *= MassFunc(mass,redshift) * mass; /* Extra mass for the log 
integration trick */

// printf("mass = %le, integrand = %le \n",mass,MassFunc(mass,redshift));
return integrand;
}


/* zeta as calculated in 
www.cita.utoronto.ca/~malvarez/halo_absorption
*/

double zeta(double mass, double flux){
  /*refer eq.3 of above */

  double answer = 1+70*pow((1+redshift)/21, -5.) * pow(mass/1e6, -0.33333) * (flux/3e7);

  return pow(answer,-1.0/3.0);
}
