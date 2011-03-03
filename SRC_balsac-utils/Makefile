FC = ifort
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML 
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
DESTDIR  = ../
EXECN1 = outnn2plt
EXECN2 = str2lat
EXECN3 = str2plt
EXECN4 = str2xyz

EXEC1   = $(DESTDIR)$(EXECN1)
EXEC2   = $(DESTDIR)$(EXECN2)
EXEC3   = $(DESTDIR)$(EXECN3)
EXEC4   = $(DESTDIR)$(EXECN4)

OBJS1 =	outnn2plt.o
OBJS2 =	str2lat.o
OBJS3 =	str2plt.o
OBJS4 =	str2xyz.o

$(EXEC1):	$(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)
	$(FC)  -o $(EXEC1) $(OBJS1) $(LDFLAGS)
	$(FC)  -o $(EXEC2) $(OBJS2) $(LDFLAGS)
	$(FC)  -o $(EXEC3) $(OBJS3) $(LDFLAGS)
	$(FC)  -o $(EXEC4) $(OBJS4) $(LDFLAGS)

#$(EXEC2):	$(OBJS) $(OBJS2)
#	$(FC) -o $(EXEC2) $(OBJS2)

clean:	  
	-rm  *.o
	@echo clean

.f.o:
	$(FC) $(FOPT) -c $<  
.f90.o:
	$(FC) $(FOPT) -c $<
