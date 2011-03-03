rm  *.o
rm: cannot remove `*.o': No such file or directory
make: [clean] Error 1 (ignored)
clean
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rhomb_in5.f   
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hex2rhomb.f   
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rhomb2hex.f   
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xyz2struct.f   
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c convert_in2.f   
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c plane.f   
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c clminter.f   
ifort -o ../rhomb_in5 rhomb_in5.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../hex2rhomb hex2rhomb.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../rhomb2hex rhomb2hex.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../xyz2struct xyz2struct.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../convert_in2 convert_in2.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../plane plane.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../clminter clminter.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
