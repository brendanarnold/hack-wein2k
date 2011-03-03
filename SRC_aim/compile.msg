rm -f *.o _tmp_.* *.P .real .complex *~ *.mod
if [ -f .complex ]; then \
   make clean; \
fi
touch .real
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./aim
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_aim'
moduls.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c moduls_tmp.f
mv moduls_tmp.o moduls.o
rm moduls_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gener.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c aim.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main1.f
outin.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outin_tmp.f
mv outin_tmp.o outin.o
rm outin_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radial.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reduc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rsur.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotato.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef1.f
sternb.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sternb_tmp.f
mv sternb_tmp.o sternb.o
rm sternb_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sum.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vnorm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dtylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c grhosphe.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dradial.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gen_change.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotat_back.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate_back.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gen_brav.f
interst.frc: REAL version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c interst_tmp.f
mv interst_tmp.o interst.o
rm interst_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphere.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c init.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c onestep.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rsurf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gauleg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c critic.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c surf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c integrho.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c follow.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c interp.f
icc -c cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c jacobi.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ludcmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spline.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c doit.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dtdtylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumdd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hrhosphe.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matprod.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c readcs.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dipole.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rho.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cossin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c intlib.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c davint.f
icc -c d1mach.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sortag.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c follown.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c SplineFit.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c CheckCSpline.f
ifort -o ./aim moduls.o charge.o  gbass.o  gener.o  aim.o main1.o  outin.o  radial.o  reduc.o  rsur.o rotat.o  rotate.o  rotato.o  rotdef1.o sternb.o sum.o  vnorm.o  ylm.o  dtylm.o  sumd.o  grhosphe.o dradial.o  gen_change.o  rotat_back.o  rotate_back.o gen_brav.o  interst.o  sphere.o  init.o  onestep.o rsurf.o  gauleg.o  critic.o  surf.o integrho.o follow.o interp.o cputim.o jacobi.o ludcmp.o spline.o doit.o dtdtylm.o sumdd.o hrhosphe.o matprod.o readcs.o dipole.o rho.o cossin.o  intlib.o davint.o d1mach.o sortag.o follown.o SplineFit.o CheckCSpline.o   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_aim'
if [ -f .real ]; then \
   make clean; \
fi
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_aim'
rm -f *.o _tmp_.* *.P .real .complex *~ *.mod
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_aim'
touch .complex
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./aimc
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_aim'
moduls.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c moduls_tmp.f
mv moduls_tmp.o moduls.o
rm moduls_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gener.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c aim.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main1.f
outin.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outin_tmp.f
mv outin_tmp.o outin.o
rm outin_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c radial.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reduc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rsur.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotato.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef1.f
sternb.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sternb_tmp.f
mv sternb_tmp.o sternb.o
rm sternb_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sum.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vnorm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dtylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c grhosphe.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dradial.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gen_change.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotat_back.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate_back.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gen_brav.f
interst.frc: COMPLEX version extracted
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c interst_tmp.f
mv interst_tmp.o interst.o
rm interst_tmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphere.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c init.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c onestep.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rsurf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gauleg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c critic.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c surf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c integrho.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c follow.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c interp.f
icc -c cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c jacobi.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ludcmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spline.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c doit.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dtdtylm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumdd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hrhosphe.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matprod.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c readcs.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dipole.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rho.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cossin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c intlib.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c davint.f
icc -c d1mach.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sortag.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c follown.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c SplineFit.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c CheckCSpline.f
ifort -o ./aimc moduls.o charge.o  gbass.o  gener.o  aim.o main1.o  outin.o  radial.o  reduc.o  rsur.o rotat.o  rotate.o  rotato.o  rotdef1.o sternb.o sum.o  vnorm.o  ylm.o  dtylm.o  sumd.o  grhosphe.o dradial.o  gen_change.o  rotat_back.o  rotate_back.o gen_brav.o  interst.o  sphere.o  init.o  onestep.o rsurf.o  gauleg.o  critic.o  surf.o integrho.o follow.o interp.o cputim.o jacobi.o ludcmp.o spline.o doit.o dtdtylm.o sumdd.o hrhosphe.o matprod.o readcs.o dipole.o rho.o cossin.o  intlib.o davint.o d1mach.o sortag.o follown.o SplineFit.o CheckCSpline.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_aim'
