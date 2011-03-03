SHELL = /bin/sh
FC = ifort
CC = icc
LINKER = $(FC)
FOPT = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN =
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
LIBS = -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
DESTDIR	= ./
EXECNAME = qtl

###############################################################################
###############################################################################

FFLAGS = $(FGEN) $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS    = reallocate.o modules.o angle.o abcd.o atpar.o dvbes1.o \
          errclr.o errflg.o \
          gtfnam2.o harmon.o inispl.o \
          inverss.o l2main.o qtlmain.o  latgen.o outerr.o \
          outwin.o ph.o psplit.o rint13.o \
          rotate.o sphbes.o symop.o \
          ylm.o cputim.o mknam.o diracout.o dergl.o\
          xsplt.o sym.o readc.o outp.o qtl.o\
          radint.o rout.o lomain.o \
         inth.o inouh.o

#..............................................................................

$(EXEC): $(OBJS)
	$(LINKER) -o $(EXEC) $(OBJS) $(LDFLAGS) $(LIBS) 

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):	Makefile

#..............................................................................
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#
.f.o:
	$(FC) $(FFLAGS) -c $<
.c.o:
	$(CC)  -c $<


.f90.o:
	$(FC) $(FFLAGS) -c $<
