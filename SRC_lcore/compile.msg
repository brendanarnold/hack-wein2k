rm -f *.o *_tmp_.* *.P *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cdsld.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dalp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dfradial.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fcore.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fpot.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hfsd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inouh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c insld.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inth.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c potsl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c resld.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sevald.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c somm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spline.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ykdir.f
ifort -o ./lcore cdsld.o  charge.o  dalp.o   dfradial.o errclr.o errflg.o  fcore.o fpot.o   gtfnam.o  hfsd.o   inouh.o    insld.o  inth.o    outerr.o potsl.o  resld.o   sevald.o somm.o     spline.o ykdir.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
