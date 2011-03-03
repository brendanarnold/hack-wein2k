rm  *.o
rm: cannot remove `*.o': No such file or directory
make: [clean] Error 1 (ignored)
clean
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c phonon2wien2k_sym.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wien2k2phonon_sym.f
ifort -o ../phonon2wien2k_sym phonon2wien2k_sym.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../wien2k2phonon_sym wien2k2phonon_sym.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
