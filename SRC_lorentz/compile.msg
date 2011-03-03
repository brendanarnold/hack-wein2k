rm  -f lorentz.o gtfnam.o outerr.o errclr.o errflg.o   lorentz.P gtfnam.P outerr.P errclr.P errflg.P lorentz.prj gtfnam.prj outerr.prj errclr.prj errflg.prj  \
	lorentz.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lorentz.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort -o ./lorentz lorentz.o gtfnam.o outerr.o errclr.o errflg.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
