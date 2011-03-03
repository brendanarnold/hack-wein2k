rm  -f felder.o chrct.o   classe.o  cputim.o  crwcnd.o  dprodci.o  errclr.o  errflg.o eulang.o  gtfnam.o init.o    invmat.o  invmati.o  kgroup.o  kptin.o matmm2.o  mdfpg.o   latgen2.o opnfs.o   pntgrp.o   raxis.o   rotkv.o rmprop.o  su2op.o   trsyma.o   trsymb.o   wrtir.o   irrep.o  wrtdate.o felder.P chrct.P   classe.P  cputim.P  crwcnd.P  dprodci.P  errclr.P  errflg.P eulang.P  gtfnam.P init.P    invmat.P  invmati.P  kgroup.P  kptin.P matmm2.P  mdfpg.P   latgen2.P opnfs.P   pntgrp.P   raxis.P   rotkv.P rmprop.P  su2op.P   trsyma.P   trsymb.P   wrtir.P   irrep.P  wrtdate.P chrct.prj   classe.prj  cputim.prj  crwcnd.prj  dprodci.prj  errclr.prj  errflg.prj eulang.prj  gtfnam.prj init.prj    invmat.prj  invmati.prj  kgroup.prj  kptin.prj matmm2.prj  mdfpg.prj   latgen2.prj opnfs.prj   pntgrp.prj   raxis.prj   rotkv.prj rmprop.prj  su2op.prj   trsyma.prj   trsymb.prj   wrtir.prj   irrep.prj  wrtdate.prj felder.prj \
irrep .xref *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c felder.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c chrct.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c classe.f
icc -c cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c crwcnd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dprodci.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c eulang.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c init.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c invmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c invmati.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kgroup.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kptin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matmm2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mdfpg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c opnfs.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pntgrp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c raxis.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotkv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rmprop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c su2op.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c trsyma.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c trsymb.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wrtir.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c irrep.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wrtdate.f
ifort  -o ./irrep  felder.o chrct.o   classe.o  cputim.o  crwcnd.o  dprodci.o  errclr.o  errflg.o eulang.o  gtfnam.o init.o    invmat.o  invmati.o  kgroup.o  kptin.o matmm2.o  mdfpg.o   latgen2.o opnfs.o   pntgrp.o   raxis.o   rotkv.o rmprop.o  su2op.o   trsyma.o   trsymb.o   wrtir.o   irrep.o  wrtdate.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread     
