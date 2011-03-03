rm  -f reallocate.o struct.o rotmat.o waves.o errclr.o  errflg.o  deter.o   kdelta.o  gtfnam.o  recpr.o  rotdef.o sphbes.o  stern.o   latgen.o  rgen.o    somm1.o   stars1.o  ylm.o kapp.o    rotate.o  somm2.o   dstart.o  outerr.o  gbass.o  matmm.o cputim.o  sortag.o  ifflim.o reallocate.P struct.P rotmat.P waves.P errclr.P  errflg.P  deter.P   kdelta.P  gtfnam.P  recpr.P  rotdef.P sphbes.P  stern.P   latgen.P  rgen.P    somm1.P   stars1.P  ylm.P kapp.P    rotate.P  somm2.P   dstart.P  outerr.P  gbass.P  matmm.P cputim.P  sortag.P  ifflim.P errclr.prj  errflg.prj  deter.prj   kdelta.prj  gtfnam.prj  recpr.prj  rotdef.prj sphbes.prj  stern.prj   latgen.prj  rgen.prj    somm1.prj   stars1.prj  ylm.prj kapp.prj    rotate.prj  somm2.prj   dstart.prj  outerr.prj  gbass.prj  matmm.prj cputim.prj  sortag.prj  ifflim.prj reallocate.prj struct.prj rotmat.prj waves.prj \
dstart.xref *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c struct.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c waves.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c deter.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kdelta.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c recpr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c somm1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stars1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kapp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c somm2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dstart.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matmm.f
icc    -c -o cputim.o cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sortag.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ifflim.f
ifort -o ./dstart reallocate.o struct.o rotmat.o waves.o errclr.o  errflg.o  deter.o   kdelta.o  gtfnam.o  recpr.o  rotdef.o sphbes.o  stern.o   latgen.o  rgen.o    somm1.o   stars1.o  ylm.o kapp.o    rotate.o  somm2.o   dstart.o  outerr.o  gbass.o  matmm.o cputim.o  sortag.o  ifflim.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
