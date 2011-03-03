rm  -f atomgrid.o reallocate.o atpar.o  cdsld.o   charge.o  dalp.o errclr.o errflg.o  fpot.o gtfnam.o hfsd.o    inilpw.o  inouh.o insld.o  inth.o    outerr.o  outwin.o potsl.o  resld.o   rint13.o somm.o   txspec.o  ykdir.o atomgrid.P reallocate.P atpar.P  cdsld.P   charge.P  dalp.P errclr.P errflg.P  fpot.P gtfnam.P hfsd.P    inilpw.P  inouh.P insld.P  inth.P    outerr.P  outwin.P potsl.P  resld.P   rint13.P somm.P   txspec.P  ykdir.P atpar.prj  cdsld.prj   charge.prj  dalp.prj errclr.prj errflg.prj  fpot.prj gtfnam.prj hfsd.prj    inilpw.prj  inouh.prj insld.prj  inth.prj    outerr.prj  outwin.prj potsl.prj  resld.prj   rint13.prj somm.prj   txspec.prj  ykdir.prj atomgrid.prj reallocate.prj \
txspec.xref *.mod
make: Circular atomgrid.o <- atomgrid.o dependency dropped.
make: Circular reallocate.o <- atomgrid.o dependency dropped.
make: Circular reallocate.o <- reallocate.o dependency dropped.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atomgrid.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cdsld.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dalp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fpot.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hfsd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inilpw.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inouh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c insld.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inth.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c potsl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c resld.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c somm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c txspec.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ykdir.f
ifort -o ./txspec atomgrid.o reallocate.o atpar.o  cdsld.o   charge.o  dalp.o errclr.o errflg.o  fpot.o gtfnam.o hfsd.o    inilpw.o  inouh.o insld.o  inth.o    outerr.o  outwin.o potsl.o  resld.o   rint13.o somm.o   txspec.o  ykdir.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
