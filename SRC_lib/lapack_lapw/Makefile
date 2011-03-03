###############################################################################
#                                                                             #
#  make              ... generate LAPACK library for lapw1                    #
#  make clean        ... remove unnecessary files                             #
#  make check        ... check static semantics of the FORTRAN routines       #
#                                                                             #
###############################################################################

#------------------------------------------------------------------------------
#
#                          Generic version
#
# FC ........ compiler name
# FOPT ...... optimization flags
# FGEN ...... code generation flags (flags not related to optimization)
# DESTDIR ... directory where the library should go (including a trailing '/' !)
# LIBNAME ... name of the resulting library (without prefixed path !)
#
#------------------------------------------------------------------------------
FC = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     =
DESTDIR  = ../
LIBNAME  = liblapack_lapw.a

###############################################################################
###############################################################################

FFLAGS = $(FGEN) $(FOPT)
LIB = $(DESTDIR)$(LIBNAME)

#..............................................................................
#
#  routines which are to be included into the library
#
LIBOBJS =	$(LIB)(lapack_lapw1.o) \
		$(LIB)(lapack_lapw1_2k.o) \
		$(LIB)(ztrtri.o) \
		$(LIB)(ztrti2.o) \
		$(LIB)(zhpev.o) \
		$(LIB)(zhpevx.o) \
		$(LIB)(zladiv.o) \
		$(LIB)(zlazro.o) \
		$(LIB)(dgeqrf.o) \
		$(LIB)(dgeqr2.o) \
		$(LIB)(dlarfb.o) \
		$(LIB)(dlarft.o) \
		$(LIB)(dlapy2.o) \
		$(LIB)(dsbtrd.o) \
		$(LIB)(dlargv.o) \
		$(LIB)(dlartv.o) \
		$(LIB)(dlar2v.o) \
		$(LIB)(dsyev.o) \
		$(LIB)(dlapy3.o)
#..............................................................................
#
#  Filenames necessary for static semantics checking.
#  These filenames have to consist of the respective routine name where the
#  extension is changed to 'prj'.
#
PRJS =	lapack_lapw1.prj 

#..............................................................................
#
.SUFFIXES:	.prj

#..............................................................................
#
#  generate library
#
$(LIB):	$(LIBOBJS)
	ar rv $(LIB) $?
	ranlib $(LIB)
#	rm $?
	@echo
	@echo Library \'$(LIB)\' is up-to-date now.

#..............................................................................
#
#  perform static semantics checking (variables used before set, type
#  consistency, correct usage of common blocks, etc.)
#
check:  $(PRJS)
	ftnchek -novice=2 -f77 $(PRJS)

#..............................................................................
#
#  remove object files, listing files and files needed for
#  static semantics checking
#
clean:
	rm  -f *.o *.list *.prj ../liblapack_lapw.a

#..............................................................................
#
#  define rules for generating files needed to check the static semantics
#  of the FORTRAN source code and for generating object files to include into
#  a library
#
.f.o:
	$(FC) -c $(FFLAGS) $<
.f.a:
	$(FC) -c $(FFLAGS) $<

.f.prj:
	ftnchek -project -f77 -declare -sixchar -novice=2 \
		-portabel -noextern -library $<













