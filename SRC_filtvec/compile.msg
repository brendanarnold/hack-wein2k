rm -f *.o _tmp_.* *.P .real .complex
if [ -f .complex ]; then \
   make clean; \
fi
touch .real
cp -p param.inc_r param.inc
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./filtvec
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_filtvec'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c filtvec.f
main.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main.frc_tmp_.f
mv main.frc_tmp_.o main.o
rm main.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c info.f
ifort -o ./filtvec filtvec.o main.o info.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_filtvec'
rm param.inc
if [ -f .real ]; then \
   make clean; \
fi
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_filtvec'
rm -f *.o _tmp_.* *.P .real .complex
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_filtvec'
touch .complex
cp -p param.inc_c param.inc
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./filtvecc
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_filtvec'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c filtvec.f
main.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main.frc_tmp_.f
mv main.frc_tmp_.o main.o
rm main.frc_tmp_.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c info.f
ifort -o ./filtvecc filtvec.o main.o info.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_filtvec'
rm param.inc
