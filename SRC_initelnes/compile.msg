rm  -f initelnes.o  gtfnam.o  strln.o errclr.o      errflg.o  outerr.o getsplit.o   initelnes.P  gtfnam.P  strln.P errclr.P      errflg.P  outerr.P getsplit.P initelnes.prj  gtfnam.prj  strln.prj errclr.prj      errflg.prj  outerr.prj getsplit.prj  \
	initelnes.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c initelnes.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c strln.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c getsplit.f
ifort -o ./initelnes initelnes.o  gtfnam.o  strln.o errclr.o      errflg.o  outerr.o getsplit.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
