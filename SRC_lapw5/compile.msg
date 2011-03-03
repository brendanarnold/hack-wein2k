rm -f *.o _tmp_.* *.P .real .complex *.mod
if [ -f .complex ]; then \
   make clean; \
fi
touch .real
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./lapw5
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw5'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atomgrid.f
wavegrid.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wavegrid_tmp.f
mv wavegrid_tmp.o wavegrid.o
rm wavegrid_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpos.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gener.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hslgr1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hslgrs.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c interp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw5.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main1.f
outin.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outin_tmp.f
mv outin_tmp.o outin.o
rm outin_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c prtrh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radial.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reduc.f
rhoout.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rhoout_tmp.f
mv rhoout_tmp.o rhoout.o
rm rhoout_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotato.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef1.f
sternb.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sternb_tmp.f
mv sternb_tmp.o sternb.o
rm sternb_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sum.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vdif.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vnorm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort -o ./lapw5 atomgrid.o wavegrid.o reallocate.o atpos.o charg.o  charge.o gbass.o  gener.o hslgr1.o hslgrs.o interp.o lapw5.o  main1.o outin.o prtrh.o  radial.o reduc.o  rhoout.o rotat.o  rotate.o rotato.o rotdef1.o sternb.o sum.o    vdif.o   vnorm.o  ylm.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw5'
if [ -f .real ]; then \
   make clean; \
fi
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw5'
rm -f *.o _tmp_.* *.P .real .complex *.mod
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw5'
touch .complex
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./lapw5c
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw5'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atomgrid.f
wavegrid.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wavegrid_tmp.f
mv wavegrid_tmp.o wavegrid.o
rm wavegrid_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpos.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gener.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hslgr1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hslgrs.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c interp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw5.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main1.f
outin.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outin_tmp.f
mv outin_tmp.o outin.o
rm outin_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c prtrh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radial.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reduc.f
rhoout.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rhoout_tmp.f
mv rhoout_tmp.o rhoout.o
rm rhoout_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotato.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef1.f
sternb.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sternb_tmp.f
mv sternb_tmp.o sternb.o
rm sternb_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sum.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vdif.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vnorm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort -o ./lapw5c atomgrid.o wavegrid.o reallocate.o atpos.o charg.o  charge.o gbass.o  gener.o hslgr1.o hslgrs.o interp.o lapw5.o  main1.o outin.o prtrh.o  radial.o reduc.o  rhoout.o rotat.o  rotate.o rotato.o rotdef1.o sternb.o sum.o    vdif.o   vnorm.o  ylm.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw5'
