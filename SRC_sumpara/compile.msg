rm  -f errclr.o errflg.o gtfnam2.o dmatscf.o latgen.o sumpara.o  outerr.o scfsum.o mknam.o strln.o   errclr.P errflg.P gtfnam2.P dmatscf.P latgen.P sumpara.P  outerr.P scfsum.P mknam.P strln.P errclr.prj errflg.prj gtfnam2.prj dmatscf.prj latgen.prj sumpara.prj  outerr.prj scfsum.prj mknam.prj strln.prj  \
	sumpara.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dmatscf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumpara.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c scfsum.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mknam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c strln.f
ifort -o ./sumpara errclr.o errflg.o gtfnam2.o dmatscf.o latgen.o sumpara.o  outerr.o scfsum.o mknam.o strln.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
