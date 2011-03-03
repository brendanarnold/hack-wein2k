FC = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     =
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread     
DESTDIR  = ./
EXECN1 = structgen

EXEC   = $(DESTDIR)$(EXECN1)

OBJS =	module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o \
        structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o \
        scan_octave.o

$(EXEC):	$(OBJS) 
	$(FC) -o $(EXEC) $(FFLAGS) $(LDFLAGS) $(OBJS)

$(OBJS):	Makefile

clean:	  
	rm -f $(OBJS) structgen
	@echo clean

.f.o:
	$(FC) -c $(FOPT) $<
.f90.o:
	$(FC) $(FOPT) -c $<
