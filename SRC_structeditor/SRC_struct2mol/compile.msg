rm -f *.o _tmp_.* *.P .real .complex *.mod
ifort  -FR -mp -w -c reallocate.f
ifort  -FR -mp -w -c module.f
ifort  -FR -mp -w -c gtfnam.f
ifort  -FR -mp -w -c errflg.f
ifort  -FR -mp -w -c errclr.f
ifort  -FR -mp -w -c outerr.f
ifort  -FR -mp -w -c struct2mol.f
ifort  -FR -mp -w -c pdb.f
ifort  -FR -mp -w -c xtl.f
ifort  -FR -mp -w -c xyz.f
ifort  -FR -mp -w -c pov.f
ifort  -FR -mp -w -c dataexp.f
ifort  -FR -mp -w -c make_list_P.f
ifort  -FR -mp -w -c make_list_F.f
ifort  -FR -mp -w -c make_list_B.f
ifort  -FR -mp -w -c make_list_H.f
ifort  -FR -mp -w -c make_list_R.f
ifort  -FR -mp -w -c make_list_CXY.f
ifort  -FR -mp -w -c make_list_CXZ.f
ifort  -FR -mp -w -c make_list_CYZ.f
ifort  -FR -mp -w -c add2list.f
ifort  -FR -mp -w -c create_in.f
ifort  -FR -mp -w -c getelem.f
ifort  -FR -mp -w -c create_cf.f
ifort -o ./struct2mol reallocate.o module.o gtfnam.o errflg.o errclr.o outerr.o struct2mol.o pdb.o xtl.o xyz.o pov.o dataexp.o make_list_P.o make_list_F.o make_list_B.o make_list_H.o make_list_R.o make_list_CXY.o make_list_CXZ.o make_list_CYZ.o add2list.o create_in.o getelem.o create_cf.o -L/area51/opt/intel/mkl721/lib/32 -static-libcxa 
size ./struct2mol
   text	   data	    bss	    dec	    hex	filename
 532499	  33644	  70652	 636795	  9b77b	./struct2mol
