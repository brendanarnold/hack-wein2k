rm -f *.o _tmp_.* *.P .real .complex *.mod
if [ -f .complex ]; then \
   make clean; \
fi
touch .real
cp -p param.inc_r param.inc
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./lapw3
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw3'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atomgrid.f
kgrid.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kgrid_tmp.f
mv kgrid_tmp.o kgrid.o
rm kgrid_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ainter.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dfac.f
fourir.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fourir_tmp.f
mv fourir_tmp.o fourir.o
rm fourir_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ftrans.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kdelta.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw3.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c recpr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef5.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sinint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern.f
strfac.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c strfac_tmp.f
mv strfac_tmp.o strfac.o
rm strfac_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c suml.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumlm.f
ustphx.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ustphx_tmp.f
mv ustphx_tmp.o ustphx.o
rm ustphx_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c locdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matmm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c deter.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort -o ./lapw3 atomgrid.o kgrid.o reallocate.o ainter.o atpar.o dfac.o fourir.o ftrans.o kdelta.o lapw3.o recpr.o rotate.o rotdef5.o sinint.o stern.o strfac.o  suml.o sumlm.o ustphx.o ylm.o gbass.o locdef.o matmm.o deter.o latgen.o outerr.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw3'
rm param.inc
if [ -f .real ]; then \
   make clean; \
fi
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw3'
rm -f *.o _tmp_.* *.P .real .complex *.mod
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw3'
touch .complex
cp -p param.inc_c param.inc
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./lapw3c
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw3'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atomgrid.f
kgrid.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kgrid_tmp.f
mv kgrid_tmp.o kgrid.o
rm kgrid_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ainter.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dfac.f
fourir.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fourir_tmp.f
mv fourir_tmp.o fourir.o
rm fourir_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ftrans.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kdelta.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw3.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c recpr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef5.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sinint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern.f
strfac.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c strfac_tmp.f
mv strfac_tmp.o strfac.o
rm strfac_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c suml.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumlm.f
ustphx.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ustphx_tmp.f
mv ustphx_tmp.o ustphx.o
rm ustphx_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c locdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matmm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c deter.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort -o ./lapw3c atomgrid.o kgrid.o reallocate.o ainter.o atpar.o dfac.o fourir.o ftrans.o kdelta.o lapw3.o recpr.o rotate.o rotdef5.o sinint.o stern.o strfac.o  suml.o sumlm.o ustphx.o ylm.o gbass.o locdef.o matmm.o deter.o latgen.o outerr.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw3'
rm param.inc
