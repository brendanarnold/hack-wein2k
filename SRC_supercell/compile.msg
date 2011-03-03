rm  -f  supercell.o    supercell.P  supercell.prj   \
supercell.xref *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c supercell.f
ifort -o ./supercell  supercell.o   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
