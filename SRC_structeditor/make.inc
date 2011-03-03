FC = ifc
CC = cc
FOPT     = -FR -mp1 -w -prec_div -pc80 -pad -ip
LDFLAGS  =  -L/opt/intel/mkl70/lib/32 -static -pthread  -lguide -Vaxlib
LIBS     = -lmkl_lapack -lmkl_ia32 -lguide
