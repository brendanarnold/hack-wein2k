rm -f module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o 
clean
ifort -c -FR -mp -w module.f
ifort -c -FR -mp -w spacegroup.f
ifort -c -FR -mp -w getlattype.f
ifort -c -FR -mp -w getsgnum.f
ifort -c -FR -mp -w getsgname.f
ifort -c -FR -mp -w getzz.f
ifort -c -FR -mp -w structgen.f
ifort -c -FR -mp -w scan_cif.f
ifort -c -FR -mp -w ciftbx.f
ifort -c -FR -mp -w hash_funcs.f
ifort -c -FR -mp -w test_sgname.f
ifort -c -FR -mp -w scan_in.f
ifort -o ./structgen  -L/area51/opt/intel/mkl721/lib/32 -static-libcxa      module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o
