make clean -C SRC_ncmsymmetry 
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_ncmsymmetry'
rm -f *.o *_tmp_.* *.P *.mod ncmsymmetry
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_ncmsymmetry'
make clean -C SRC_readwrite
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_readwrite'
rm -f module.o readwrite.o convert.o rwoctave.o *.mod readwrite
clean
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_readwrite'
make clean -C SRC_struct2mol
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_struct2mol'
rm -f *.o _tmp_.* *.P .real .complex *.mod struct2mol
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_struct2mol'
make clean -C SRC_structgen 
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_structgen'
rm -f module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o scan_octave.o structgen
clean
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_structgen'
make -C SRC_ncmsymmetry 
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_ncmsymmetry'
make[1]: Circular module.o <- module.o dependency dropped.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c module.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ncmsymmetry.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c determinant.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c euler.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inversa.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_struct_klist.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c test_symetry.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c def_Ylm_rot_mat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c def_spin_rot_mat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_point_groups.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_lm_list.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fix_nonmagatom.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapack.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fix_rotloc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fiteta.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c class.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dirdeh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kurki.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lmsort.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matrot.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dirdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapack2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c order_symoper.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pglsym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c put.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matinv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ludcmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lubksb.f
ifort -o ./ncmsymmetry module.o gtfnam.o errflg.o errclr.o outerr.o ncmsymmetry.o determinant.o euler.o inversa.o make_struct_klist.o test_symetry.o def_Ylm_rot_mat.o def_spin_rot_mat.o make_point_groups.o make_lm_list.o fix_nonmagatom.o lapack.o fix_rotloc.o fiteta.o symop.o class.o dirdeh.o kurki.o lmsort.o matrot.o dirdef.o lapack2.o symgen.o  order_symoper.o pglsym.o put.o matinv.o ludcmp.o lubksb.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_ncmsymmetry'
cp SRC_ncmsymmetry/ncmsymmetry bin/
make -C SRC_readwrite
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_readwrite'
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML module.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML readwrite.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML convert.f
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML rwoctave.f
ifort -o ./readwrite  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread      module.o readwrite.o convert.o rwoctave.o
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_readwrite'
cp SRC_readwrite/readwrite bin/
make -C SRC_struct2mol
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_struct2mol'
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c module.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c struct2mol.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pdb.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xtl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xyz.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pov.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dataexp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_P.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_F.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_B.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_H.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_R.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_CXY.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_CXZ.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c make_list_CYZ.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c add2list.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c create_in.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c getelem.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c create_cf.f
ifort -o ./struct2mol reallocate.o module.o gtfnam.o errflg.o errclr.o outerr.o struct2mol.o pdb.o xtl.o xyz.o pov.o dataexp.o make_list_P.o make_list_F.o make_list_B.o make_list_H.o make_list_R.o make_list_CXY.o make_list_CXZ.o make_list_CYZ.o add2list.o create_in.o getelem.o create_cf.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
size ./struct2mol
   text	   data	    bss	    dec	    hex	filename
 744463	  29528	  81776	 855767	  d0ed7	./struct2mol
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_struct2mol'
cp SRC_struct2mol/struct2mol bin/
make -C SRC_structgen 
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_structgen'
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
ifort -c -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML scan_octave.f
ifort -o ./structgen  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread      module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o scan_octave.o
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_structeditor/SRC_structgen'
cp SRC_structgen/structgen bin/
