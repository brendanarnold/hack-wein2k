rm  *.o
rm: cannot remove `*.o': No such file or directory
make: [clean] Error 1 (ignored)
clean
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c genetempl.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c c2rh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c c2te.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xc2rh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setelast.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c anaelast.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lmdif.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fdump.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c i1mach.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c j4save.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dpolft.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dp1vlu.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xercnt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xerhlt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xermsg.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xerprn.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xersve.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xgetua.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -o ../genetempl genetempl.o c2rh.o c2te.o xc2rh.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -o ../setelast setelast.o c2rh.o c2te.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -o ../anaelast anaelast.o lmdif.o fdump.o i1mach.o j4save.o dpolft.o dp1vlu.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o xgetua.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread
