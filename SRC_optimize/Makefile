FC = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
LDFLAGS  = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
DESTDIR  = ./
EXECNAME = optimize

EXEC   = $(DESTDIR)$(EXECNAME)

OBJS =	optimize.o

$(EXEC):	$(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS)

$(OBJS):	Makefile

clean:	  
	-rm optimize $(OBJS)
	@echo clean

.f.o:
	$(FC) $(FOPT) -c $<


.f90.o:
	$(FC) $(FOPT) -c $<
