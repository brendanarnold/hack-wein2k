rm  -f broadening.o gtfnam.o outerr.o errclr.o errflg.o corebroadening.o valencebroadening.o spectrobroadening.o   broadening.P gtfnam.P outerr.P errclr.P errflg.P corebroadening.P valencebroadening.P spectrobroadening.P broadening 
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o broadening.o broadening.f
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o gtfnam.o gtfnam.f
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o outerr.o outerr.f
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o errclr.o errclr.f
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o errflg.o errflg.f
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o corebroadening.o corebroadening.f
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o valencebroadening.o valencebroadening.f
ifort   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c -o spectrobroadening.o spectrobroadening.f
ifort  -o ./broadening broadening.o gtfnam.o outerr.o errclr.o errflg.o corebroadening.o valencebroadening.o spectrobroadening.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread  
