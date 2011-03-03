FC = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     = 
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
DESTDIR  = ./
EXECNAME = eosfit

EXEC   = $(DESTDIR)$(EXECNAME)

OBJS =	eosfit.o gtfnam.o lmdif.o

$(EXEC):	$(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(FOPT) $(LDFLAGS) 

$(OBJS):	Makefile

clean:	  
	-rm eosfit *.o
	@echo clean

.f.o:
	$(FC) $(FOPT) -c $<	
.f90.o:
	$(FC) $(FOPT) -c $<
