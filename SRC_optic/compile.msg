rm -f *.o _tmp_.* *.P .real .complex *.mod
if [ -f .complex ]; then \
   make clean; \
fi
touch .real
cp -p param.inc_r param.inc
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./optic
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_optic'
modules.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp.f
mv modules_tmp.o modules.o
rm modules_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar_op.f
icc -c cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c harmon.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
lomain_op.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lomain_op_tmp.f
mv lomain_op_tmp.o lomain_op.o
rm lomain_op_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mmatrix.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c opmain.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
planew.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c planew_tmp.f
mv planew_tmp.o planew.o
rm planew_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rintuout.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
sph-UP.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sph-UP_tmp.f
mv sph-UP_tmp.o sph-UP.o
rm sph-UP_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort -o ./optic modules.o abc.o  atpar_op.o  cputim.o  dvbes1.o errclr.o errflg.o gtfnam.o harmon.o  latgen.o  lomain_op.o  mmatrix.o  opmain.o  outerr.o outwin.o  planew.o  outmat.o radint.o  rint13.o  rintuout.o rotate.o sph-UP.o  sphbes.o  symop.o  sym.o  ylm.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_optic'
rm param.inc
if [ -f .real ]; then \
   make clean; \
fi
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_optic'
rm -f *.o _tmp_.* *.P .real .complex *.mod
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_optic'
touch .complex
cp -p param.inc_c param.inc
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./opticc
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_optic'
modules.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp.f
mv modules_tmp.o modules.o
rm modules_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar_op.f
icc -c cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c harmon.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
lomain_op.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lomain_op_tmp.f
mv lomain_op_tmp.o lomain_op.o
rm lomain_op_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mmatrix.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c opmain.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
planew.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c planew_tmp.f
mv planew_tmp.o planew.o
rm planew_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rintuout.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
sph-UP.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sph-UP_tmp.f
mv sph-UP_tmp.o sph-UP.o
rm sph-UP_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort -o ./opticc modules.o abc.o  atpar_op.o  cputim.o  dvbes1.o errclr.o errflg.o gtfnam.o harmon.o  latgen.o  lomain_op.o  mmatrix.o  opmain.o  outerr.o outwin.o  planew.o  outmat.o radint.o  rint13.o  rintuout.o rotate.o sph-UP.o  sphbes.o  symop.o  sym.o  ylm.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_optic'
rm param.inc
