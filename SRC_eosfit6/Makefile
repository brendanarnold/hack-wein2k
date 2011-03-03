FC = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
FGEN     = 
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
DESTDIR  = ./
EXECNAME = eosfit6

EXEC   = $(DESTDIR)$(EXECNAME)

OBJS =	eosfit6.o gtfnam.o all.o

$(EXEC):	$(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(FOPT) $(LDFLAGS) 

$(OBJS):	Makefile

clean:	  
	-rm eosfit6 *.o
	@echo clean

.f.o:
	$(FC) $(FOPT) -c $<	
.f90.o:
	$(FC) $(FOPT) -c $<
