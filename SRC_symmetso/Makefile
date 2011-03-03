###############################################################################
#                                                                             #
#  make              ... generate executable                                  #
#  make clean        ... delete unnecessary files                             #
#                                                                             #
###############################################################################

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
FC       = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     = 
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
LIBS     =
DESTDIR  = ./
EXECNAME = symmetso

###############################################################################
###############################################################################

FFLAGS = $(FGEN) $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS1 =	symmetso.o angle.o br1dm.o checks.o class.o clmchange.o dirdef.o dirdeh.o  \
        dirlat.o dot.o gbass.o gtfnam.o inverss.o in1ch.o incch.o \
        kurki.o latgen.o latsym.o  lmsort.o locdef.o lubksb.o ludcmp.o \
        matmm.o matinv.o matrot.o minvec.o order.o ordix.o outerr.o\
	pgbsym.o pglsym.o put.o reclat.o rewr.o rotdef.o rstruc.o \
        strwri.o symso.o symho.o symop.o volnew.o volofb.o 

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
	$(FC) -o $(EXEC) $(FOPT) $(OBJS) $(LDFLAGS) $(LIBS)

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):        param.inc Makefile

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
		$(EXECNAME).xref
#..............................................................................
#
#  define rules for generating files needed to check the static semantics
#  of the FORTRAN source code and for generating object files from source files
#
.f.o:
	$(FC) $(FFLAGS) -c $<
