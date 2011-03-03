rm  -f felder.o arbdos.o  noculc.o  readop.o joint.o errclr.o errflg.o filename.o gtfnam.o outerr.o cputim.o rstru.o latgen.o invmat.o rotdef.o drude.o opt1.o dumm.o read_diag.o  felder.P arbdos.P  noculc.P  readop.P joint.P errclr.P errflg.P filename.P gtfnam.P outerr.P cputim.P rstru.P latgen.P invmat.P rotdef.P drude.P opt1.P dumm.P read_diag.P arbdos.prj  noculc.prj  readop.prj joint.prj errclr.prj errflg.prj filename.prj gtfnam.prj outerr.prj cputim.prj rstru.prj latgen.prj invmat.prj rotdef.prj drude.prj opt1.prj dumm.prj read_diag.prj felder.prj \
joint .xref *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c felder.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c arbdos.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c noculc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c readop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c joint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c filename.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
icc -c cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rstru.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c invmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drude.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c opt1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dumm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c read_diag.f
ifort -o ./joint  felder.o arbdos.o  noculc.o  readop.o joint.o errclr.o errflg.o filename.o gtfnam.o outerr.o cputim.o rstru.o latgen.o invmat.o rotdef.o drude.o opt1.o dumm.o read_diag.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
