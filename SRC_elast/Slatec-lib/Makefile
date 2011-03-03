FC = pgf90
FOPT     = -O -Mfree
FGEN     =
DESTDIR  = ../
LIBNAME  = libpolfit.a

###############################################################################
###############################################################################

FFLAGS = $(FGEN)
LIB = $(DESTDIR)$(LIBNAME)

#..............................................................................
#
#  routines which are to be included into the library
#
LIBOBJS =	$(LIB)(fdump.o) \
		$(LIB)(i1mach.o) \
		$(LIB)(j4save.o) \
		$(LIB)(dpolft.o) \
		$(LIB)(dp1vlu.o) \
		$(LIB)(xercnt.o) \
		$(LIB)(xerhlt.o) \
                $(LIB)(xermsg.o) \
                $(LIB)(xerprn.o) \
                $(LIB)(xersve.o) \
                $(LIB)(xgetua.o)
#..............................................................................

#  generate library
#
$(LIB):	$(LIBOBJS)
	ar rv $(LIB) $?
	ranlib $(LIB)
	rm $?
	@echo
	@echo Library \'$(LIB)\' is up-to-date now.

#.............................................................................

.f.o:
	$(FC) $(FOPT) -c $<
