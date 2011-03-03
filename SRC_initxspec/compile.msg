rm  -f initxspec.o  gtfnam.o  strln.o errclr.o      errflg.o  outerr.o   initxspec.P  gtfnam.P  strln.P errclr.P      errflg.P  outerr.P initxspec.prj  gtfnam.prj  strln.prj errclr.prj      errflg.prj  outerr.prj  \
	initxspec.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c initxspec.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c strln.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort -o ./initxspec initxspec.o  gtfnam.o  strln.o errclr.o      errflg.o  outerr.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
