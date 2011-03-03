FC = ifort
FOPT = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
DESTDIR  = ../
EXECN1 = rhomb_in5
EXECN2 = hex2rhomb
EXECN3 = rhomb2hex
EXECN4 = xyz2struct
EXECN5 = convert_in2
EXECN6 = plane
EXECN7 = clminter

EXEC   = $(DESTDIR)$(EXECN1)
EXEC2   = $(DESTDIR)$(EXECN2)
EXEC3   = $(DESTDIR)$(EXECN3)
EXEC4   = $(DESTDIR)$(EXECN4)
EXEC5   = $(DESTDIR)$(EXECN5)
EXEC6   = $(DESTDIR)$(EXECN6)
EXEC7   = $(DESTDIR)$(EXECN7)

OBJS =	rhomb_in5.o
OBJS2 =	hex2rhomb.o
OBJS3 =	rhomb2hex.o
OBJS4 =	xyz2struct.o
OBJS5 =	convert_in2.o
OBJS6 =	plane.o
OBJS7 =       clminter.o

$(EXEC):	$(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(OBJS7)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS)
	$(FC) -o $(EXEC2) $(OBJS2) $(LDFLAGS)
	$(FC) -o $(EXEC3) $(OBJS3) $(LDFLAGS)
	$(FC) -o $(EXEC4) $(OBJS4) $(LDFLAGS)
	$(FC) -o $(EXEC5) $(OBJS5) $(LDFLAGS)
	$(FC) -o $(EXEC6) $(OBJS6) $(LDFLAGS)
	$(FC) -o $(EXEC7) $(OBJS7) $(LDFLAGS)

#$(EXEC2):	$(OBJS) $(OBJS2)
#	$(FC) -o $(EXEC2) $(OBJS2)

clean:	  
	-rm  *.o
	@echo clean

.f.o:
	$(FC) $(FOPT) -c $<   
.f90.o:
	$(FC) $(FOPT) -c $<
