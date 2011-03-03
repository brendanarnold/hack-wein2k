rm -f spacegroup.o 
clean
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML spacegroup.f
ifort -o ./spacegroup  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread      spacegroup.o 
