ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML blas_lapw1.f
make: ifort: Command not found
make: *** [../libblas_lapw.a(blas_lapw1.o)] Error 127
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML blas_lapw1.f
ar rv ../libblas_lapw.a blas_lapw1.o
ar: creating ../libblas_lapw.a
a - blas_lapw1.o
ranlib ../libblas_lapw.a
rm blas_lapw1.o

Library '../libblas_lapw.a' is up-to-date now.
