rm  *.o
rm: cannot remove `*.o': No such file or directory
make: [clean] Error 1 (ignored)
clean
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c outnn2plt.f  
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c str2lat.f  
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c str2plt.f  
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c str2xyz.f  
ifort  -o ../outnn2plt outnn2plt.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort  -o ../str2lat str2lat.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort  -o ../str2plt str2plt.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort  -o ../str2xyz str2xyz.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread
