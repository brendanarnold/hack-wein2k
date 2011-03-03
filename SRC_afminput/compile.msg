rm  -f afminput.o   afminput.P afminput.prj  \
	afminput.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c afminput.f
ifort -o ./afminput afminput.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
