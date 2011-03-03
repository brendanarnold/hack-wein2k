rm -f module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o 
clean
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML module.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML spacegroup.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML getlattype.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML getsgnum.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML getsgname.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML getzz.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML structgen.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML scan_cif.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML ciftbx.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML hash_funcs.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML test_sgname.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML scan_in.f
ifort -o ./cif2struct  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread      module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o
