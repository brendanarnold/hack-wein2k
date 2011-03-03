###############################################################################
#                                                                             #
#  make           ... generate executable for the REAL-version                #
#  make complex   ... generate executable for the COMPLEX-version             #
#  make clean     ... delete unnecessary files                                #
#                                                                             #
#  CAUTION: This makefile must not be called via the '-f' parameter of 'make'.#
#           (Don't use 'make -f Makefile.aix32' !!)                           #
#           (Always (symbolic) link the appropriate Makefile to 'Makefile'!!) #
#                                                                             #
###############################################################################
#
# FC ........... compiler name
# FL ........... linker name (in case it is different than the compiler)
# FOPT ......... optimization flags
# FGEN ......... code generation flags (flags not related to optimization)
# LDFLAGS ...... linker flags
# R_LIBS ....... libraries needed to build the REAL executable
# C_LIBS ....... libraries needed to build the COMPLEX executable
# DESTDIR ...... dir. where the executable should go (without trailing '/'!)
# R_EXECNAME ... name of the resulting REAL executable (without prefixed path!)
# C_EXECNAME ... name of the resulting COMPLEX executable ( -"- )
#
SHELL = /bin/sh
FC = ifort
CC = icc
FOPT = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN = 
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
R_LIBS = -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
C_LIBS = $(R_LIBS)
DESTDIR	= .
R_EXECNAME = aim
C_EXECNAME = aimc

###############################################################################

TYPE='REAL'
TYPE_COMMENT='!_REAL'
FFLAGS = $(FGEN) $(FOPT)
R_EXEC = $(DESTDIR)/$(R_EXECNAME)
C_EXEC = $(DESTDIR)/$(C_EXECNAME)

#..............................................................................
#
#  Object files common to both REAL and COMPLEX type subroutines
#
OBJS =	charge.o  gbass.o  gener.o  aim.o \
        main1.o  outin.o  radial.o  reduc.o  rsur.o \
	rotat.o  rotate.o  rotato.o  rotdef1.o sternb.o \
	sum.o  vnorm.o  ylm.o  dtylm.o  sumd.o  grhosphe.o \
        dradial.o  gen_change.o  rotat_back.o  rotate_back.o \
        gen_brav.o  interst.o  sphere.o  init.o  onestep.o \
        rsurf.o  gauleg.o  critic.o  surf.o integrho.o follow.o\
        interp.o cputim.o jacobi.o ludcmp.o spline.o doit.o\
	dtdtylm.o sumdd.o hrhosphe.o matprod.o readcs.o \
        dipole.o rho.o cossin.o  intlib.o davint.o d1mach.o sortag.o \
        follown.o SplineFit.o CheckCSpline.o 

PREOBJS = moduls.o
#..............................................................................
#
#  Object files for REAL and COMPLEX versions
#
R_OBJS = $(PREOBJS) $(OBJS) 
C_OBJS = $(PREOBJS) $(OBJS)

#..............................................................................
#
#  Some source files contain both REAL and COMPLEX versions. These files
#  have the file extension '.frc' and must be treated differently.
#
.SUFFIXES:	.frc
.SUFFIXES:	.frc90

#..............................................................................
#
#  Build executable (either REAL or COMPLEX versions)
#
real: rm_c_files
	$(MAKE) TYPE='REAL' TYPE_COMMENT='!_REAL' $(R_EXEC)

complex: rm_r_files
	$(MAKE) TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' $(C_EXEC)

#..............................................................................

$(R_EXEC): $(R_OBJS)
	$(FC) -o $(R_EXEC) $(R_OBJS) $(LDFLAGS) $(R_LIBS)

$(C_EXEC): $(C_OBJS)
	$(FC) -o $(C_EXEC) $(C_OBJS) $(LDFLAGS) $(C_LIBS)
#	size $(C_EXEC)

rm_c_files:
	if [ -f .complex ]; then \
	   $(MAKE) clean; \
	fi
	touch .real

rm_r_files:
	if [ -f .real ]; then \
	   $(MAKE) clean; \
	fi
	touch .complex

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(R_OBJS):		param.inc Makefile 
$(C_OBJS):		param.inc Makefile 

#..............................................................................
#
#  some subroutines require special optimization options on HP-UX 8.05 f77
#
#hns.o:	hns.frc lohns.f rotate.f ylm.f param.inc
#	@sed -e 's/^$(TYPE_COMMENT)//'   hns.frc >_tmp_.f
#	@echo "hns.frc: $(TYPE) version extracted"
#	$(FC) +OP4 -w -WP,-unroll=6,-unroll2=24,-minvector=6 \
#		-WP,-inff='rotate.f',-inff='ylm.f',-inff='lohns.f' -c _tmp_.f
#	mv _tmp_.o hns.o
#
#hamilt.o:	hamilt.frc dvbes1.f sphbes.f rotate.f ustphx.frc \
#		warpin.f ylm.f param.inc
#	$(FC) +OP4c -WP,-unroll=6,-unroll2=24,-minvector=6 -c sphbes.f
#	mv sphbes.P _tmp1_.f
#	@sed -e 's/^$(TYPE_COMMENT)//'   ustphx.frc >_tmp2_.f
#	@echo "ustphx.frc: $(TYPE) version extracted"
#	@sed -e 's/^$(TYPE_COMMENT)//'   hamilt.frc >_tmp_.f
#	@echo "hamilt.frc: $(TYPE) version extracted"
#	$(FC) +OP4 -w -WP,-unroll=6,-unroll2=24,-minvector=6 \
#		-WP,-inff='dvbes1.f',-inff='_tmp1_.f',-inff='rotate.f' \
#		-WP,-inff='ylm.f',-inff='_tmp2_.f',-inff='warpin.f' -c _tmp_.f
#	mv _tmp_.o hamilt.o
#	rm _tmp1_.f _tmp2_.f
#
#sphbes.o:	sphbes.f param.inc
#	$(FC) +OP4 -w -WP,-unroll=6,-unroll2=24,-minvector=6 -c sphbes.f

#..............................................................................
#
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	-rm -f *.o _tmp_.* *.P .real .complex *~ *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#  (From .frc files the appropriatly typed subroutine (REAL/COMPLEX) has to
#  be extracted first.)
#
.f.o:
	$(FC) $(FFLAGS) -c $<
.frc.o:
	@sed -e 's/^$(TYPE_COMMENT)//'   $< >$*_tmp.f
	@echo "$<: $(TYPE) version extracted"
	$(FC) $(FFLAGS) -c $*_tmp.f
	mv $*_tmp.o $*.o
	rm $*_tmp.f
.f.mod:
	$(FC) $(FFLAGS) -c $<  
.c.o:
	$(CC) -c $<
.f90.o:
	$(FC) $(FFLAGS) -c $<
.frc90.o:
	@sed -e 's/^$(TYPE_COMMENT)//'   $< >_tmp_.f90
	@echo "$<: $(TYPE) version extracted"
	$(FC) $(FFLAGS) -c _tmp_.f90
	mv _tmp_.o $*.o
