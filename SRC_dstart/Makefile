###############################################################################
#                                                                             #
#  make              ... generate executable                                  #
#  make clean        ... delete unnecessary files                             #
#  make check        ... check static semantics of the FORTRAN routines       #
#  make xref         ... generate a cross reference listing                   #
#                                                                             #
###############################################################################

#------------------------------------------------------------------------------
#
#                          IBM RS6000  AIX 3.2
#
# FC ......... compiler name
# FOPT ....... optimization flags
# FGEN ....... code generation flags (flags not related to optimization)
# LDFLAGS .... linker flags
# LIBS ....... libraries needed to build the executable
# DESTDIR .... directory where the executable should go (incl. trailing '/'!)
# EXECNAME ... name of the resulting executable (without prefixed path!)
#
#------------------------------------------------------------------------------
FC = ifort
CC = icc
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     = 
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
LIBS     =
DESTDIR  = ./
EXECNAME = dstart

###############################################################################
###############################################################################

FFLAGS = $(FGEN) $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS1 =	errclr.o  errflg.o  deter.o   kdelta.o  gtfnam.o  recpr.o  rotdef.o\
	sphbes.o  stern.o   latgen.o  rgen.o    somm1.o   stars1.o  ylm.o\
	kapp.o    rotate.o  somm2.o   dstart.o  outerr.o  gbass.o  matmm.o \
        cputim.o  sortag.o  ifflim.o

#..............................................................................
#
#  Routines which may require preprocessing before compiling on some machines
#  (Sourcecode-level optimization)
#
PREOBJS = reallocate.o struct.o rotmat.o waves.o

#..............................................................................
#
OBJS = $(PREOBJS) $(OBJS1)
PRJS = $(OBJS1:.o=.prj) $(PREOBJS:.o=.prj)

#..............................................................................
#
.SUFFIXES:	.prj

#..............................................................................
#
#  Build executable
#
$(EXEC):	$(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS) $(LIBS)

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):	param.inc Makefile 

#..............................................................................
#
#  perform static semantics checking (variables used before set, type
#  consistency, correct usage of common blocks, etc.)
#
check:	$(PRJS)
	ftnchek $(PRJS)

#..............................................................................
#
#  remove object files, preprocessed source files and files needed for
#  static semantics checking
#
clean:
	rm  -f $(OBJS) $(PREOBJS:.o=.P) $(OBJS1:.o=.P) $(PRJS) \
		$(EXECNAME).xref *.mod

#..............................................................................
#
#  build a cross reference table (file '$(EXECNAME).xref')
#
xref:
	fxref $(OBJS1:.o=.f) $(PREOBJS:.o=.f) >$(EXECNAME).xref

#..............................................................................
#
#  define rules for generating files needed to check the static semantics
#  of the FORTRAN source code and for generating object files from source files
#
.f.prj:
	ftnchek -project -f77 -declare -sixchar -novice=2 -noextern -library $<
.f.o:
	$(FC) $(FFLAGS) -c $<
.f.mod:
	$(FC) $(FFLAGS) -c $<
.f90.o:
	$(FC) $(FFLAGS) -c $<
