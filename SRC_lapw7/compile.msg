rm -f *.o _tmp_.* *.P .real .complex
if [ -f .complex ]; then \
   make clean; \
fi
touch .real
cp -p param.inc_r param.inc
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./lapw7
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw7'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c moduls.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw7.f
main.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main.frc_tmp_.f
mv main.frc_tmp_.o main.o
rm main.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c locdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c orth.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c grdgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spcgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c findmt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c auggen.f
augpw.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c augpw.frc_tmp_.f
mv augpw.frc_tmp_.o augpw.o
rm augpw.frc_tmp_.f
auglo.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c auglo.frc_tmp_.f
mv auglo.frc_tmp_.o auglo.o
rm auglo.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wavsph.f
wavint.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wavint.frc_tmp_.f
mv wavint.frc_tmp_.o wavint.o
rm wavint.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c trans.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c bessel.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort -o ./lapw7 moduls.o lapw7.o  main.o   latgen.o rotdef2.o locdef.o orth.o   grdgen.o spcgen.o findmt.o auggen.o augpw.o  auglo.o  wavsph.o wavint.o trans.o  bessel.o gbass.o  outwin.o rint13.o ylm.o    sphbes.o dvbes1.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw7'
rm param.inc
if [ -f .real ]; then \
   make clean; \
fi
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw7'
rm -f *.o _tmp_.* *.P .real .complex
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw7'
touch .complex
cp -p param.inc_c param.inc
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./lapw7c
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw7'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c moduls.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw7.f
main.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main.frc_tmp_.f
mv main.frc_tmp_.o main.o
rm main.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c locdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c orth.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c grdgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spcgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c findmt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c auggen.f
augpw.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c augpw.frc_tmp_.f
mv augpw.frc_tmp_.o augpw.o
rm augpw.frc_tmp_.f
auglo.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c auglo.frc_tmp_.f
mv auglo.frc_tmp_.o auglo.o
rm auglo.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wavsph.f
wavint.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wavint.frc_tmp_.f
mv wavint.frc_tmp_.o wavint.o
rm wavint.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c trans.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c bessel.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort -o ./lapw7c moduls.o lapw7.o  main.o   latgen.o rotdef2.o locdef.o orth.o   grdgen.o spcgen.o findmt.o auggen.o augpw.o  auglo.o  wavsph.o wavint.o trans.o  bessel.o gbass.o  outwin.o rint13.o ylm.o    sphbes.o dvbes1.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw7'
rm param.inc
