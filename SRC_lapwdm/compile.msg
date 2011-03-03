rm -f *.o _tmp_.* *.P .real .complex *.mod
if [ -f .complex ]; then \
   make clean; \
fi
touch .real
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./lapwdm
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapwdm'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
modules.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp.f
mv modules_tmp.o modules.o
rm modules_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c angle.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c euler.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c diracout.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dergl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c harmon.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inispl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inverss.f
l2main.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c l2main_tmp.f
mv l2main_tmp.o l2main.o
rm l2main_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lmtmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapwdm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c orb.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ph.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inouh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inth.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c t3j.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symmetrization.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
icc    -c -o cputim.o cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mknam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xsplt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c locmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c timeinv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c readw.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rout.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lomain.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c couplx.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c printx.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xoper.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xqproduct.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symmso.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c addtinvso.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c addtinv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c output.f
ifort -o ./lapwdm reallocate.o modules.o abc.o angle.o atpar.o dvbes1.o dmat.o errclr.o errflg.o euler.o diracout.o dergl.o gtfnam2.o harmon.o inispl.o inverss.o l2main.o lmtmat.o lapwdm.o  latgen.o outerr.o outwin.o orb.o ph.o rint13.o inouh.o inth.o rotate.o rotdef.o sphbes.o  symop.o t3j.o symmetrization.o ylm.o cputim.o mknam.o xsplt.o locmat.o timeinv.o sym.o readw.o radint.o rout.o lomain.o spin.o couplx.o printx.o xoper.o xqproduct.o symmso.o addtinvso.o addtinv.o output.o   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapwdm'
if [ -f .real ]; then \
   make clean; \
fi
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapwdm'
rm -f *.o _tmp_.* *.P .real .complex *.mod
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapwdm'
touch .complex
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./lapwdmc
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapwdm'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
modules.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp.f
mv modules_tmp.o modules.o
rm modules_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c angle.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c euler.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c diracout.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dergl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c harmon.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inispl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inverss.f
l2main.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c l2main_tmp.f
mv l2main_tmp.o l2main.o
rm l2main_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lmtmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapwdm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c orb.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ph.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inouh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inth.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c t3j.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symmetrization.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
icc    -c -o cputim.o cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mknam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xsplt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c locmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c timeinv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c readw.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rout.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lomain.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c couplx.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c printx.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xoper.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xqproduct.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symmso.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c addtinvso.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c addtinv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c output.f
ifort -o ./lapwdmc reallocate.o modules.o abc.o angle.o atpar.o dvbes1.o dmat.o errclr.o errflg.o euler.o diracout.o dergl.o gtfnam2.o harmon.o inispl.o inverss.o l2main.o lmtmat.o lapwdm.o  latgen.o outerr.o outwin.o orb.o ph.o rint13.o inouh.o inth.o rotate.o rotdef.o sphbes.o  symop.o t3j.o symmetrization.o ylm.o cputim.o mknam.o xsplt.o locmat.o timeinv.o sym.o readw.o radint.o rout.o lomain.o spin.o couplx.o printx.o xoper.o xqproduct.o symmso.o addtinvso.o addtinv.o output.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapwdm'
