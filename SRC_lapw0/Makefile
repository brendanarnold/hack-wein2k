.SUFFIXES:        .F
.SUFFIXES:        .F90
SHELL = /bin/sh
FC = ifort
MPF = ifort
CC = icc
FOPT =  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FPOPT =  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
DParallel = '-DParallel'
FGEN = $(PARALLEL)
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
RP_LIBS = -L /usr/local/SCALAPACK -L /usr/local/BLACS/LIB -lpblas -lredist -ltools -lscalapack -lfblacs -lblacs -lmpi

S_EXECNAME = lapw0
P_EXECNAME = lapw0_mpi

###############################################################################
###############################################################################

S_EXEC   = $(DESTDIR)$(S_EXECNAME)
P_EXEC   = $(DESTDIR)$(P_EXECNAME)
#..............................................................................
#
#  All object files, including $(SP_OBJS) and $(OBJSMOD) 
#
OBJS  = cputim.o modules.o reallocate.o ainv.o blyp.o charg2.o  charg3.o \
        charge.o chfac.o chslv.o  corgga.o  cub_xc_back.o corlsd.o \
        dlu.o drho.o  dfxhpbe.o dylm.o  efg.o  energy.o errclr.o errflg.o \
        ev92.o    ev92ex.o exch.o exch17.o exrevpbe.o  fithi.o fxhpbe.o \
        gbass.o gcor.o  gea.o geaex.o  getfft.o  getff1.o gpoint.o  grans.o \
        gtfnam.o hcth.o ifflim.o kcis.o lapw0.o  latgen.o  multfc.o multsu.o \
        outerr.o  pbe1.o poissn.o potfac.o pwxad4.o pwxad5.o  qranf.o \
        readstruct.o rean0.o rean1.o rean3.o  rean4.o  rotate.o  rotdef.o \
        setff0.o  setff1.o  setfft.o setff2.o seval.o sevald.o sevaldd.o \
        sevali.o sevalin.o sicpbe.o sphbes.o  spline.o srolyl.o stern.o \
        sumfac.o  suml.o th1.o th2.o vpw91.o   vresp.o vs98.o vxc15.o \
        vxc16.o vxc17.o vxc24.o  vxc26.o  vxclm2.o vxcpw2.o vxi35.o vxi35a.o \
        workf1.o xcener.o  xcpot.o  xcpot1.o xcpot3.o ykav.o  ylm.o zfft3d.o \
        sicpbe_tpss.o corpbe_tpss.o fx_tpss.o dfxtpss.o rhopw.o b88.o wc05.o \
        eramps.o epot1.o rpbe.o

# object file for parallel and sequential
SP_OBJS = modules.o reallocate.o energy.o gtfnam.o lapw0.o 

# Object files using modules
OBJSMOD = drho.o energy.o getff1.o gtfnam.o latgen.o lapw0.o multsu.o \
          poissn.o pwxad4.o readstruct.o rean0.o rean4.o rotdef.o setff1.o \
          setff2.o setfft.o stern.o xcpot1.o xcpot3.o ykav.o  

#........................................................................
default: seq

all:    
	make seq para

last: 
	if [ ! -f .sequential -a ! -f .parallel ]; then make seq; fi
	if [ -f .sequential ]; then make seq; fi
	if [ -f .parallel ]; then make para; fi
        
seq: keep_s_files
	$(MAKE) $(S_EXEC) FORT=$(FC) FFLAGS='$(FGEN) $(FOPT)'

para: keep_p_files
	$(MAKE) PARALLEL=$(DParallel) $(P_EXEC) \
	  FORT=$(MPF) FFLAGS='$(FGEN) $(FPOPT) $(DParallel)'

#..............................................................................

$(S_EXEC): $(OBJS)
	$(FORT) -o $(S_EXEC) $(OBJS) $(LDFLAGS) 

$(P_EXEC): $(OBJS)
	$(FORT) -o $(P_EXEC) $(OBJS) $(LDFLAGS) $(RP_LIBS) 

# sequential
keep_s_files:
	if [ -f .parallel ]; then \
	   rm -f .parallel $(SP_OBJS) *.mod; \
	fi
	touch .sequential

# parallel
keep_p_files:
	if [ -f .sequential ]; then \
	   rm -f .sequential $(SP_OBJS) *.mod; \
	fi
	touch .parallel

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):		param.inc Makefile
$(OBJSMOD):		modules.o Makefile


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
	$(FORT) $(FFLAGS) -c $<
.F.o:
	$(FORT) $(FFLAGS) -c $<
.c.o:
	$(CC) -c $<

.f90.o:
	$(FORT) $(FFLAGS) -c $<
.F90.o:
	$(FORT) $(FFLAGS) -c $<
