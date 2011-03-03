rm eosfit *.o
rm: cannot remove `eosfit': No such file or directory
rm: cannot remove `*.o': No such file or directory
make: [clean] Error 1 (ignored)
clean
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c eosfit.f	
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f	
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lmdif.f	
ifort -o ./eosfit eosfit.o gtfnam.o lmdif.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
