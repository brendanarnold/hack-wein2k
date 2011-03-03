FC = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     =
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread     
DESTDIR  = ./
EXECN1 = readwrite

EXEC   = $(DESTDIR)$(EXECN1)

OBJS =	module.o readwrite.o convert.o rwoctave.o

$(EXEC):	$(OBJS) 
	$(FC) -o $(EXEC) $(FFLAGS) $(LDFLAGS) $(OBJS)

$(OBJS):	Makefile

clean:	  
	rm -f $(OBJS) *.mod readwrite
	@echo clean

.f.o:
	$(FC) -c $(FOPT) $<
.f90.o:
	$(FC) $(FOPT) -c $<
