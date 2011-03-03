rm -f *.o
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c energy.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c findneigh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c findsymm.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fixup.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gen_brav.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c init.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c makehess.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sortag.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pickneigh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dchdc.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c expandset.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c makehessb.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c findequivs.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matcon.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c patch.f
ifort -o ./pairhess energy.o findneigh.o findsymm.o fixup.o gen_brav.o init.o main.o makehess.o sortag.o pickneigh.o dchdc.o  errclr.o errflg.o gtfnam.o outerr.o expandset.o makehessb.o findequivs.o matcon.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
ifort -o ./patchsymm patch.o fixup.o gtfnam.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
