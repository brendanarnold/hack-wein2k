rm  *.o
rm: cannot remove `*.o': No such file or directory
make: [clean] Error 1 (ignored)
clean
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c bcc_fs_mesh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fcc_fs_mesh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hex_fs_mesh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sc_fs_mesh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spagh2rho.f
ifort -o ../bcc_fs_mesh bcc_fs_mesh.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../fcc_fs_mesh fcc_fs_mesh.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../hex_fs_mesh hex_fs_mesh.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../sc_fs_mesh sc_fs_mesh.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -o ../spagh2rho spagh2rho.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
