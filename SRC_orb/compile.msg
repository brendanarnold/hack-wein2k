rm  -f modules.o main.o init.o angle.o opnfs.o  errclr.o  errflg.o gtfnam.o latgen2.o orbmom.o  outwin.o parop.o rint13.o charg2.o charg3.o rotate.o vldau.o Vcalc.o t3j.o vorbmix.o       \
orb    .xref *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c modules.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c main.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c init.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c angle.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c opnfs.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c latgen2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c orbmom.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c outwin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c parop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c charg2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c charg3.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c vldau.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c Vcalc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c t3j.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -c vorbmix.f
ifort -o ./orb     modules.o main.o init.o angle.o opnfs.o  errclr.o  errflg.o gtfnam.o latgen2.o orbmom.o  outwin.o parop.o rint13.o charg2.o charg3.o rotate.o vldau.o Vcalc.o t3j.o vorbmix.o    -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML      -L/opt/intel/mkl/9.0/lib/em64t -lpthread  
