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
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML 
FGEN     = 
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
LIBS     =
DESTDIR  = ./
EXECNAME = clmcopy
EXECNAM1 = struct_afm_check

###############################################################################
###############################################################################

FFLAGS = $(FGEN) $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)
EXE1   = $(DESTDIR)$(EXECNAM1)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS =	clmcopy.o
OBJ1 =  modules.o  reallocate.o struct_afm_check.o readstruct.o writestruct.o outerr.o

#..............................................................................
#
#  Build executable
#
$(EXEC):        $(OBJS) $(OBJ1)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS) $(LIBS)
	$(FC) -o $(EXE1) $(OBJ1) $(LDFLAGS) $(LIBS)

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS): Makefile param.inc
$(OBJ1): Makefile modules.o reallocate.o

#..............................................................................
#
#  remove object files, preprocessed source files and files needed for
#  static semantics checking
#
clean:
	rm  -f $(OBJS) $(OBJS1) 

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
.f90.o:
	$(FC) $(FFLAGS) -c $<
