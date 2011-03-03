rm  -f reallocate.o modules.o bndind.o bz_lin.o cartco.o clipin.o  comprel.o con_ev.o defins.o drawt.o  get_ei.o get_k.o inview.o movet.o  pgrpnr.o pointi.o  psend.o  psinit.o spag.o seppt.o transt.o writln.o writs.o  writz.o wrtdate.o reallocate.P modules.P bndind.P bz_lin.P cartco.P clipin.P  comprel.P con_ev.P defins.P drawt.P  get_ei.P get_k.P inview.P movet.P  pgrpnr.P pointi.P  psend.P  psinit.P spag.P seppt.P transt.P writln.P writs.P  writz.P wrtdate.P modules.prj bndind.prj bz_lin.prj cartco.prj clipin.prj  comprel.prj con_ev.prj defins.prj drawt.prj  get_ei.prj get_k.prj inview.prj movet.prj  pgrpnr.prj pointi.prj  psend.prj  psinit.prj spag.prj seppt.prj transt.prj writln.prj writs.prj  writz.prj wrtdate.prj reallocate.prj \
	spaghetti.xref
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c bndind.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c bz_lin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cartco.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c clipin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c comprel.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c con_ev.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c defins.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drawt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c get_ei.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c get_k.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inview.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c movet.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pgrpnr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pointi.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c psend.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c psinit.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spag.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seppt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c transt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c writln.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c writs.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c writz.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wrtdate.f
ifort -o ./spaghetti reallocate.o modules.o bndind.o bz_lin.o cartco.o clipin.o  comprel.o con_ev.o defins.o drawt.o  get_ei.o get_k.o inview.o movet.o  pgrpnr.o pointi.o  psend.o  psinit.o spag.o seppt.o transt.o writln.o writs.o  writz.o wrtdate.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
