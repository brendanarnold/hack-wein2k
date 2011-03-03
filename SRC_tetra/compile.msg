rm  -f reallocate.o tetra.o arbdos.o noculc.o ados.o gtfnam.o errflg.o errclr.o broad.o reallocate.P tetra.P arbdos.P noculc.P ados.P gtfnam.P errflg.P errclr.P broad.P tetra.prj arbdos.prj noculc.prj ados.prj gtfnam.prj errflg.prj errclr.prj broad.prj reallocate.prj \
	tetra.xref
make: Circular reallocate.o <- reallocate.o dependency dropped.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c tetra.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c arbdos.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c noculc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c ados.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c broad.f
ifort -o ./tetra reallocate.o tetra.o arbdos.o noculc.o ados.o gtfnam.o errflg.o errclr.o broad.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
