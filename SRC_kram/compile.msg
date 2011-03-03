rm  -f drude.o kram.o eshift.o green.o green_reg.o gtfnam.o losmo.o losmo1.o losmo2.o   drude.P kram.P eshift.P green.P green_reg.P gtfnam.P losmo.P losmo1.P losmo2.P drude.prj kram.prj eshift.prj green.prj green_reg.prj gtfnam.prj losmo.prj losmo1.prj losmo2.prj  \
	kram .xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drude.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kram.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c eshift.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c green.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c green_reg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c losmo.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c losmo1.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c losmo2.f
ifort -o ./kram  drude.o kram.o eshift.o green.o green_reg.o gtfnam.o losmo.o losmo1.o losmo2.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
