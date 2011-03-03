rm  -f atomgrid.o mxpmgrid.o dfpmin.o  errclr.o  errflg.o  finish.o  func.o    gbass.o gtfnam.o  haupt.o   ininos.o  inv.o     latgen.o  lnsrch.o maxstp.o  mini.o    mold0.o   mold.o    nn.o      nose.o nose0.o   nwtmin.o  outerr.o  pairdis.o rotate.o  rotdef.o wrtscf.o  latmix.o  distxy.o mdble.o d1mach.o ditsum.o dr7mdc.o  dv7cpy.o    i7mdcn.o divset.o drldst.o  dv7dfl.o  dl7itv.o drmng.o   dv7scp.o stopx.o dl7ivm.o dv7vmp.o dl7tvm.o drmngwien.o dw7zbf.o da7sst.o dl7upd.o func2.o  dd7dog.o dl7vml.o  dv2axy.o    traptouch.o dd7tpr.o dparck.o  dv2nrm.o  constr.o atomgrid.P mxpmgrid.P dfpmin.P  errclr.P  errflg.P  finish.P  func.P    gbass.P gtfnam.P  haupt.P   ininos.P  inv.P     latgen.P  lnsrch.P maxstp.P  mini.P    mold0.P   mold.P    nn.P      nose.P nose0.P   nwtmin.P  outerr.P  pairdis.P rotate.P  rotdef.P wrtscf.P  latmix.P  distxy.P mdble.P d1mach.P ditsum.P dr7mdc.P  dv7cpy.P    i7mdcn.P divset.P drldst.P  dv7dfl.P  dl7itv.P drmng.P   dv7scp.P stopx.P dl7ivm.P dv7vmp.P dl7tvm.P drmngwien.P dw7zbf.P da7sst.P dl7upd.P func2.P  dd7dog.P dl7vml.P  dv2axy.P    traptouch.P dd7tpr.P dparck.P  dv2nrm.P  constr.P dfpmin.prj  errclr.prj  errflg.prj  finish.prj  func.prj    gbass.prj gtfnam.prj  haupt.prj   ininos.prj  inv.prj     latgen.prj  lnsrch.prj maxstp.prj  mini.prj    mold0.prj   mold.prj    nn.prj      nose.prj nose0.prj   nwtmin.prj  outerr.prj  pairdis.prj rotate.prj  rotdef.prj wrtscf.prj  latmix.prj  distxy.prj mdble.prj d1mach.prj ditsum.prj dr7mdc.prj  dv7cpy.prj    i7mdcn.prj divset.prj drldst.prj  dv7dfl.prj  dl7itv.prj drmng.prj   dv7scp.prj stopx.prj dl7ivm.prj dv7vmp.prj dl7tvm.prj drmngwien.prj dw7zbf.prj da7sst.prj dl7upd.prj func2.prj  dd7dog.prj dl7vml.prj  dv2axy.prj    traptouch.prj dd7tpr.prj dparck.prj  dv2nrm.prj  constr.prj atomgrid.prj mxpmgrid.prj \
	mini.xref
make: Circular atomgrid.o <- atomgrid.o dependency dropped.
make: Circular mxpmgrid.o <- atomgrid.o dependency dropped.
make: Circular mxpmgrid.o <- mxpmgrid.o dependency dropped.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mxpmgrid.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atomgrid.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dfpmin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c finish.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c func.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c haupt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ininos.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lnsrch.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c maxstp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mini.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mold0.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mold.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c nn.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c nose.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c nose0.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c nwtmin.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pairdis.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wrtscf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latmix.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c distxy.f
icc -c mdble.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c d1mach.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ditsum.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dr7mdc.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dv7cpy.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c i7mdcn.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c divset.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drldst.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dv7dfl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dl7itv.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drmng.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dv7scp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stopx.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dl7ivm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dv7vmp.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dl7tvm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drmngwien.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dw7zbf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c da7sst.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dl7upd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c func2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dd7dog.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dl7vml.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dv2axy.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c traptouch.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dd7tpr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dparck.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dv2nrm.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c constr.f
ifort -o ./mini atomgrid.o mxpmgrid.o dfpmin.o  errclr.o  errflg.o  finish.o  func.o    gbass.o gtfnam.o  haupt.o   ininos.o  inv.o     latgen.o  lnsrch.o maxstp.o  mini.o    mold0.o   mold.o    nn.o      nose.o nose0.o   nwtmin.o  outerr.o  pairdis.o rotate.o  rotdef.o wrtscf.o  latmix.o  distxy.o mdble.o d1mach.o ditsum.o dr7mdc.o  dv7cpy.o    i7mdcn.o divset.o drldst.o  dv7dfl.o  dl7itv.o drmng.o   dv7scp.o stopx.o dl7ivm.o dv7vmp.o dl7tvm.o drmngwien.o dw7zbf.o da7sst.o dl7upd.o func2.o  dd7dog.o dl7vml.o  dv2axy.o    traptouch.o dd7tpr.o dparck.o  dv2nrm.o  constr.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread
