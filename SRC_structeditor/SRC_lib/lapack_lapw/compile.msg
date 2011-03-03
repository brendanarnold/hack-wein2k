pgf77 -c  -O lapack_lapw1.f
pgf77 -c  -O zhpev.f
pgf77 -c  -O zhpevx.f
pgf77 -c  -O zladiv.f
pgf77 -c  -O zlazro.f
pgf77 -c  -O dlapy2.f
pgf77 -c  -O dlapy3.f
ar rv ../liblapack_lapw.a lapack_lapw1.o zhpev.o zhpevx.o zladiv.o zlazro.o dlapy2.o dlapy3.o
a - lapack_lapw1.o
a - zhpev.o
a - zhpevx.o
a - zladiv.o
a - zlazro.o
a - dlapy2.o
a - dlapy3.o
ranlib ../liblapack_lapw.a
rm lapack_lapw1.o zhpev.o zhpevx.o zladiv.o zlazro.o dlapy2.o dlapy3.o

Library '../liblapack_lapw.a' is up-to-date now.
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML lapack_lapw1.f
make: ifort: Command not found
make: *** [../liblapack_lapw.a(lapack_lapw1.o)] Error 127
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML lapack_lapw1.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML lapack_lapw1_2k.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML ztrtri.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML ztrti2.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML zhpev.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML zhpevx.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML zladiv.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML zlazro.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dgeqrf.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dgeqr2.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dlarfb.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dlarft.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dlapy2.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dsbtrd.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dlargv.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dlartv.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dlar2v.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dsyev.f
ifort -c  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML dlapy3.f
ar rv ../liblapack_lapw.a lapack_lapw1.o lapack_lapw1_2k.o ztrtri.o ztrti2.o zhpev.o zhpevx.o zladiv.o zlazro.o dgeqrf.o dgeqr2.o dlarfb.o dlarft.o dlapy2.o dsbtrd.o dlargv.o dlartv.o dlar2v.o dsyev.o dlapy3.o
ar: creating ../liblapack_lapw.a
a - lapack_lapw1.o
a - lapack_lapw1_2k.o
a - ztrtri.o
a - ztrti2.o
a - zhpev.o
a - zhpevx.o
a - zladiv.o
a - zlazro.o
a - dgeqrf.o
a - dgeqr2.o
a - dlarfb.o
a - dlarft.o
a - dlapy2.o
a - dsbtrd.o
a - dlargv.o
a - dlartv.o
a - dlar2v.o
a - dsyev.o
a - dlapy3.o
ranlib ../liblapack_lapw.a

Library '../liblapack_lapw.a' is up-to-date now.
