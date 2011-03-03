rm  -f main.o sdefl.o matmm.o sdef.o groups.o bravai.o arbmsh.o basdiv.o reduz.o zuord.o tetdiv.o tetcnt.o ord1.o indexx.o gbass.o def0.o drval.o initi.o initdr.o divisi.o addinv.o   main.P sdefl.P matmm.P sdef.P groups.P bravai.P arbmsh.P basdiv.P reduz.P zuord.P tetdiv.P tetcnt.P ord1.P indexx.P gbass.P def0.P drval.P initi.P initdr.P divisi.P addinv.P main.prj sdefl.prj matmm.prj sdef.prj groups.prj bravai.prj arbmsh.prj basdiv.prj reduz.prj zuord.prj tetdiv.prj tetcnt.prj ord1.prj indexx.prj gbass.prj def0.prj drval.prj initi.prj initdr.prj divisi.prj addinv.prj  \
	kgen.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c main.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sdefl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matmm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c groups.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c bravai.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c arbmsh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c basdiv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reduz.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c zuord.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c tetdiv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c tetcnt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ord1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c indexx.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c def0.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drval.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c initi.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c initdr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c divisi.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c addinv.f
ifort -o ./kgen main.o sdefl.o matmm.o sdef.o groups.o bravai.o arbmsh.o basdiv.o reduz.o zuord.o tetdiv.o tetcnt.o ord1.o indexx.o gbass.o def0.o drval.o initi.o initdr.o divisi.o addinv.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
