/*! FILE : readers.c

Contains all the test functions required to read almost 
all kinds of information, i.e., parametre files, input
data etc.,

Author:
Rajat Mani Thomas

Enquiries: thomas@cita.utoronto.ca
Private: rajatthomas@gmail.com

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"
#include "allvars.h"

/*! This function parses the parameterfile in a simple way.  Each paramater
 *  is defined by a keyword (`tag'), and can be either of type double, int,
 *  or character string.  The routine makes sure that each parameter
 *  appears exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
int readparams(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag = 0;

  if(sizeof(long long) != 8)
    {
      //if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      return FAIL;//endrun(0);
    }

  if(sizeof(int) != 4)
    {
      //if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      return FAIL;//endrun(0);
    }

  if(sizeof(float) != 4)
    {
      //if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      return FAIL;//endrun(0);
    }

  if(sizeof(double) != 8)
    {
      //if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      return FAIL;//endrun(0);
    }

 

  //if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;
        
        strcpy(tag[nt], "Omega_Lambda");
	addr[nt] = &Omega_Lambda;
	id[nt++] = DOUBLE;

	strcpy(tag[nt], "Omega_Baryon");
        addr[nt] = &Omega_Baryon;
        id[nt++] = DOUBLE;

      	strcpy(tag[nt], "Omega_Matter");
      	addr[nt] = &Omega_Matter;
      	id[nt++] = DOUBLE;

	strcpy(tag[nt], "Omega_Neutrino");
        addr[nt] = &Omega_Neutrino;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "degen_Neutrino");
        addr[nt] = &degen_Neutrino;
        id[nt++] = DOUBLE;

	strcpy(tag[nt], "tilt");
        addr[nt] = &tilt;
        id[nt++] = DOUBLE;

	strcpy(tag[nt], "Sigma_8");
      	addr[nt] = &Sigma_8;
      	id[nt++] = DOUBLE;

	strcpy(tag[nt], "Hubble_z0");
      	addr[nt] = &Hubble_z0;
       	id[nt++] = DOUBLE;

	strcpy(tag[nt], "h_z0");
       	addr[nt] = &h_z0;
      	id[nt++] = DOUBLE;

	strcpy(tag[nt], "TranFuncType");
      	addr[nt] = &TranFuncType;
      	id[nt++] = STRING;

	strcpy(tag[nt], "MassFuncType");
      	addr[nt] = &MassFuncType;
      	id[nt++] = STRING;

	strcpy(tag[nt], "SigmaMethod");
        addr[nt] = &SigmaMethod;
        id[nt++] = STRING;

	strcpy(tag[nt], "SigmaFilter");
        addr[nt] = &SigmaFilter;
        id[nt++] = STRING;
	       
        strcpy(tag[nt], "z_end");
        addr[nt] = &z_end;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "z_start");
        addr[nt] = &z_start;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "z_stride");
        addr[nt] = &z_stride;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "mass_min");
        addr[nt] = &mass_min;
        id[nt++] = DOUBLE;

	strcpy(tag[nt], "mass_max");
        addr[nt] = &mass_max;
        id[nt++] = DOUBLE;

	strcpy(tag[nt], "mass_bins");
        addr[nt] = &mass_bins;
        id[nt++] = INT;

	strcpy(tag[nt], "CalcMoment");
        addr[nt] = &CalcMoment;
        id[nt++] = INT;

	strcpy(tag[nt], "CalcMassFunc");
        addr[nt] = &CalcMassFunc;
        id[nt++] = INT;


      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);
	    
	    }
	}
      else
	{
	  printf("\nParameter file %s not found.\n\n", fname);
	  errorFlag = 2;
	}

      if(errorFlag != 2)
	for(i = 0; i < nt; i++)
	  {
	    if(*tag[i])
	      {
		printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
		errorFlag = 1;
	      }
	  }

      
    }


#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS

    return SUCCESS;
}


