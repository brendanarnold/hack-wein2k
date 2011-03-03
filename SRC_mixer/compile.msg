rm -f *.o *_tmp_.* *.P *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c distan.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c eqival.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hyper.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c iffpar.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c integr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mixer.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c nn.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c normal.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c qmix5.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c read_denmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rean0a.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rean3.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c scfana.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setffx.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c forcint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c invmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c write_denmat.f
ifort -o ./mixer charge.o distan.o eqival.o errclr.o errflg.o gtfnam.o hyper.o iffpar.o integr.o latgen.o mixer.o  nn.o     normal.o outerr.o qmix5.o  read_denmat.o rean0a.o  rean3.o  scfana.o setffx.o stern.o forcint.o invmat.o write_denmat.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
