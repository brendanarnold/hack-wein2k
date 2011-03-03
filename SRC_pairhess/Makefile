SHELL = /bin/sh
FC = ifort
FOPT =  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
R_LIBS = -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
DESTDIR = ./
EXECNAME = pairhess
EXECN2 = patchsymm

###############################################################################
###############################################################################

EXEC   = $(DESTDIR)$(EXECNAME)
EXEC2  = $(DESTDIR)$(EXECN2)
#..............................................................................
#
#  All object files, including $(SP_OBJS) and $(OBJSMOD) 
#
OBJS  = energy.o findneigh.o findsymm.o fixup.o gen_brav.o init.o main.o makehess.o sortag.o \
        pickneigh.o dchdc.o  errclr.o errflg.o gtfnam.o outerr.o expandset.o makehessb.o \
        findequivs.o matcon.o
OBJ2  = patch.o fixup.o gtfnam.o
#..............................................................................

$(EXEC): $(OBJS) $(OBJ2)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS) $(R_LIBS)
	$(FC) -o $(EXEC2) $(OBJ2) $(LDFLAGS) $(R_LIBS)
#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)

$(OBJS): Makefile params.inc

#..............................................................................
#  remove unnecessary files 
clean:
	rm -f *.o
#..............................................................................
#
#  define inference rules to generate object files from source files
#
.f.o:
	$(FC) $(FOPT) -c $<
