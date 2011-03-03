rm  -f clmcopy.o  
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c clmcopy.f
make: Circular modules.o <- modules.o dependency dropped.
make: Circular reallocate.o <- modules.o dependency dropped.
make: Circular reallocate.o <- reallocate.o dependency dropped.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c modules.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c struct_afm_check.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c readstruct.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c writestruct.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c outerr.f
ifort -o ./clmcopy clmcopy.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread  
ifort -o ./struct_afm_check modules.o  reallocate.o struct_afm_check.o readstruct.o writestruct.o outerr.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread  
