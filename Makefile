#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make".      
#                                                                      
#----------------------------------------------------------------------

#OPT += -DMAKETABLES

#----------------------------------------------------------------------
# Here, select compile environment for the target machine. This may need 
# Local system dependencies might be involved. Follow example below.
# Following options should work on an X86_64 machine.
#----------------------------------------------------------------------

CC       =  gcc               # sets the C-compiler

OPTIMIZE =  -O2 -Wall -g   # sets optimization and warning flags

OPTIONS =  $(OPTIMIZE) $(OPT)

EXEC   = massfunc.exe

OBJS   = main.o sigma.o cosmology.o \
	 allvars.o readparams.o moments.o

INCL   = proto.h  allvars.h Makefile

GSL_INCL =  -I/data/users/harker/gsl-1.9/include 
GSL_LIBS =  -L/data/users/harker/gsl-1.9/lib


CFLAGS = $(OPTIONS) $(GSL_INCL)


LIBS   =  -lgsl -lgslcblas -lm  $(GSL_LIBS)

# Implicit Rules
.SUFFIXES: .o .c


.c.o:
	$(CC) -c $< -o $*.o $(CFLAGS)


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f *~ $(OBJS) $(EXEC)

