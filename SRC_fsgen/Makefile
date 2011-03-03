FC = ifort
FOPT = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
LDFLAGS = $(FOPT) -L/opt/intel/mkl/9.0/lib/em64t -lpthread
DESTDIR  = ../
EXECN1 = bcc_fs_mesh
EXECN2 = fcc_fs_mesh
EXECN3 = hex_fs_mesh
EXECN4 = sc_fs_mesh
EXECN5 = spagh2rho

EXEC   = $(DESTDIR)$(EXECN1)
EXEC2   = $(DESTDIR)$(EXECN2)
EXEC3   = $(DESTDIR)$(EXECN3)
EXEC4   = $(DESTDIR)$(EXECN4)
EXEC5   = $(DESTDIR)$(EXECN5)

OBJS =	bcc_fs_mesh.o
OBJS2 =	fcc_fs_mesh.o
OBJS3 =	hex_fs_mesh.o
OBJS4 =	sc_fs_mesh.o
OBJS5 =	spagh2rho.o

$(EXEC):	$(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS)
	$(FC) -o $(EXEC2) $(OBJS2) $(LDFLAGS)
	$(FC) -o $(EXEC3) $(OBJS3) $(LDFLAGS)
	$(FC) -o $(EXEC4) $(OBJS4) $(LDFLAGS)
	$(FC) -o $(EXEC5) $(OBJS5) $(LDFLAGS)

#$(EXEC2):	$(OBJS) $(OBJS2)
#	$(FC) -o $(EXEC2) $(OBJS2)

clean:	  
	-rm  *.o
	@echo clean

.f.o: 
	$(FC) $(FOPT) -c $<

.f90.o:
	$(FC) $(FOPT) -c $<
