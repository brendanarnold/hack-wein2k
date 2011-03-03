rm  -f symmetry.o symop.o class.o kurki.o latsym.o checks.o matmm.o dirlat.o pgbsym.o pglsym.o reclat.o minvec.o matinv.o put.o order.o ordix.o lmsort.o volnew.o lubksb.o ludcmp.o volofb.o dot.o rstruc.o gtfnam.o dirdef.o dirdeh.o matrot.o strwri.o   symmetry.P symop.P class.P kurki.P latsym.P checks.P matmm.P dirlat.P pgbsym.P pglsym.P reclat.P minvec.P matinv.P put.P order.P ordix.P lmsort.P volnew.P lubksb.P ludcmp.P volofb.P dot.P rstruc.P gtfnam.P dirdef.P dirdeh.P matrot.P strwri.P symmetry.prj symop.prj class.prj kurki.prj latsym.prj checks.prj matmm.prj dirlat.prj pgbsym.prj pglsym.prj reclat.prj minvec.prj matinv.prj put.prj order.prj ordix.prj lmsort.prj volnew.prj lubksb.prj ludcmp.prj volofb.prj dot.prj rstruc.prj gtfnam.prj dirdef.prj dirdeh.prj matrot.prj strwri.prj  \
	symmetry.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symmetry.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c symop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c class.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kurki.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latsym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c checks.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matmm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dirlat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pgbsym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pglsym.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reclat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c minvec.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matinv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c put.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c order.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ordix.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lmsort.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c volnew.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lubksb.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ludcmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c volofb.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dot.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rstruc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dirdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dirdeh.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matrot.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c strwri.f
ifort -o ./symmetry symmetry.o symop.o class.o kurki.o latsym.o checks.o matmm.o dirlat.o pgbsym.o pglsym.o reclat.o minvec.o matinv.o put.o order.o ordix.o lmsort.o volnew.o lubksb.o ludcmp.o volofb.o dot.o rstruc.o gtfnam.o dirdef.o dirdeh.o matrot.o strwri.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
