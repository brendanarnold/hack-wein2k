FC = ifort
FOPT = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
DESTDIR  = ../
EXECN1 = phonon2wien2k_sym
EXECN2 = wien2k2phonon_sym

EXEC   = $(DESTDIR)$(EXECN1)
EXEC2   = $(DESTDIR)$(EXECN2)

OBJS =	phonon2wien2k_sym.o
OBJS2 =	wien2k2phonon_sym.o

$(EXEC):	$(OBJS) $(OBJS2)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS)
	$(FC) -o $(EXEC2) $(OBJS2) $(LDFLAGS)

#$(EXEC2):	$(OBJS) $(OBJS2)
#	$(FC) -o $(EXEC2) $(OBJS2)

clean:	  
	-rm  *.o
	@echo clean

.f.o: 
	$(FC) $(FOPT) -c $<

.f90.o:
	$(FC) $(FOPT) -c $<
