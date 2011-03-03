rm optimize optimize.o
rm: cannot remove `optimize': No such file or directory
rm: cannot remove `optimize.o': No such file or directory
make: [clean] Error 1 (ignored)
clean
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c optimize.f
ifort -o ./optimize optimize.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
