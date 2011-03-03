rm  -f struk.o variable_fields.o nn.o dirlat.o ord2.o reduce_alloc.o  struk.P variable_fields.P nn.P dirlat.P ord2.P reduce_alloc.P nn.prj dirlat.prj ord2.prj reduce_alloc.prj struk.prj variable_fields.prj \
nn.xref *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c struk.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c variable_fields.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c nn.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c dirlat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c ord2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -c reduce_alloc.f
ifort -o ./nn struk.o variable_fields.o nn.o dirlat.o ord2.o reduce_alloc.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML  -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
