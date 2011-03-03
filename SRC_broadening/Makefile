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
FOPT     =  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     = 
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
LIBS     =  
DESTDIR  = ./
EXECNAME = broadening

###############################################################################
###############################################################################

FFLAGS = $(FGEN) $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS1   = broadening.o gtfnam.o outerr.o errclr.o errflg.o corebroadening.o \
          valencebroadening.o spectrobroadening.o

#..............................................................................
#
#  Routines which may require preprocessing before compiling on some machines
#  (Sourcecode-level optimization)
#
PREOBJS =

#..............................................................................
#
OBJS = $(OBJS1) $(PREOBJS)
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
$(OBJS): Makefile

#..............................................................................
#
#  remove object files, preprocessed source files and files needed for
#  static semantics checking
#
clean:
	rm  -f $(OBJS) $(PREOBJS:.o=.P) $(OBJS1:.o=.P) $(EXECNAME) 

