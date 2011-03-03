rm -f *.o *_tmp_.* *.P .real .complex .sequential .parallel *.mod
if [ -f .complex ]; then \
   rm -f .complex calkpt.o hamilt.o hns.o horb.o modules.o seclit.o seclr4.o seclr5.o tapewf.o ustphx.o wfpnt.o wfpnt1.o param.inc *.mod; \
fi
touch .real
cp -p param.inc_r param.inc
if [ -f .parallel ]; then \
   rm -f .parallel calkpt.o gtfnam.o hamilt.o hldau.o hns.o modules.o pdsyevx16.o prtkpt.o pzheevx16.o seclit.o seclr4.o seclr5.o inikpt.o inilpw.o lapw1.o pdsyr2m.o pzher2m.o setkpt.o *.mod; \
fi
touch .sequential
make TYPE='REAL' TYPE_COMMENT='\!_REAL' \
  ./lapw1 FORT=ifort FFLAGS=' -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML'
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw1'
modules.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp_.F
mv modules_tmp_.o modules.o
rm modules_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c bandv1.f
calkpt.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c calkpt_tmp_.F
mv calkpt_tmp_.o calkpt.o
rm calkpt_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cbcomb.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c coors.f
icc -c cputim.c
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dblr2k.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dgeqrl.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dgewy.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dgewyg.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dlbrfg.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsbein1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dscgst.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dstebz2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsyevx2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsymm2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsyr2m.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsyrb4.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsyrb5l.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsyrdt4.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsytrd2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsywyv.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dsyxev4.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c eisps.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c forfhs.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gaunt1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gaunt2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
gtfnam.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam_tmp_.F
mv gtfnam_tmp_.o gtfnam.o
rm gtfnam_tmp_.F
hamilt.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hamilt_tmp_.F
mv hamilt_tmp_.o hamilt.o
rm hamilt_tmp_.F
hns.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hns_tmp_.F
mv hns_tmp_.o hns.o
rm hns_tmp_.F
horb.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c horb_tmp_.F
mv horb_tmp_.o horb.o
rm horb_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inikpt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inilpw.f
lapw1.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw1_tmp_.F
mv lapw1_tmp_.o lapw1.o
rm lapw1_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lmsort.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c locdef.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lohns.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lopw.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c matmm.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c nn.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
pdsyevx16.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pdsyevx16_tmp_.F
mv pdsyevx16_tmp_.o pdsyevx16.o
rm pdsyevx16_tmp_.F
prtkpt.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c prtkpt_tmp_.F
mv prtkpt_tmp_.o prtkpt.o
rm prtkpt_tmp_.F
prtres.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c prtres_tmp_.F
mv prtres_tmp_.o prtres.o
rm prtres_tmp_.F
pzheevx16.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pzheevx16_tmp_.F
mv pzheevx16_tmp_.o pzheevx16.o
rm pzheevx16_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rdswar.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
seclit.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seclit_tmp_.F
mv seclit_tmp_.o seclit.o
rm seclit_tmp_.F
seclr4.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seclr4_tmp_.F
mv seclr4_tmp_.o seclr4.o
rm seclr4_tmp_.F
seclr5.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seclr5_tmp_.F
mv seclr5_tmp_.o seclr5.o
rm seclr5_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c select.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c service.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setkpt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setwar.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
stern.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern_tmp_.F
mv stern_tmp_.o stern.o
rm stern_tmp_.F
tapewf.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c tapewf_tmp_.F
mv tapewf_tmp_.o tapewf.o
rm tapewf_tmp_.F
ustphx.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ustphx_tmp_.F
mv ustphx_tmp_.o ustphx.o
rm ustphx_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vectf.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c warpin.f
wfpnt.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wfpnt_tmp_.F
mv wfpnt_tmp_.o wfpnt.o
rm wfpnt_tmp_.F
wfpnt1.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wfpnt1_tmp_.F
mv wfpnt1_tmp_.o wfpnt1.o
rm wfpnt1_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c zhcgst.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c zheevx2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c zhemm2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c zher2m.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c zhetrd2.f
pdsyr2m.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pdsyr2m_tmp_.F
mv pdsyr2m_tmp_.o pdsyr2m.o
rm pdsyr2m_tmp_.F
pzher2m.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pzher2m_tmp_.F
mv pzher2m_tmp_.o pzher2m.o
rm pzher2m_tmp_.F
ifort  -o ./lapw1 abc.o atpar.o bandv1.o calkpt.o cbcomb.o coors.o cputim.o dblr2k.o dgeqrl.o dgewy.o dgewyg.o dlbrfg.o dsbein1.o dscgst.o dstebz2.o dsyevx2.o dsymm2.o dsyr2m.o dsyrb4.o dsyrb5l.o dsyrdt4.o dsytrd2.o dsywyv.o dsyxev4.o dvbes1.o eisps.o errclr.o errflg.o forfhs.o gaunt1.o gaunt2.o gbass.o gtfnam.o hamilt.o hns.o horb.o inikpt.o inilpw.o lapw1.o latgen.o lmsort.o locdef.o lohns.o lopw.o matmm.o modules.o nn.o outerr.o outwin.o pdsyevx16.o prtkpt.o prtres.o pzheevx16.o rdswar.o rint13.o rotate.o rotdef.o seclit.o seclr4.o seclr5.o select.o service.o setkpt.o setwar.o sphbes.o stern.o tapewf.o ustphx.o vectf.o warpin.o wfpnt.o wfpnt1.o ylm.o zhcgst.o zheevx2.o zhemm2.o zher2m.o zhetrd2.o pdsyr2m.o pzher2m.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw1'
rm param.inc
if [ -f .real ]; then \
   rm -f .real calkpt.o hamilt.o hns.o horb.o modules.o seclit.o seclr4.o seclr5.o tapewf.o ustphx.o wfpnt.o wfpnt1.o param.inc *.mod; \
fi
touch .complex
cp -p param.inc_c param.inc
if [ -f .parallel ]; then \
   rm -f .parallel calkpt.o gtfnam.o hamilt.o hldau.o hns.o modules.o pdsyevx16.o prtkpt.o pzheevx16.o seclit.o seclr4.o seclr5.o inikpt.o inilpw.o lapw1.o pdsyr2m.o pzher2m.o setkpt.o *.mod; \
fi
touch .sequential
make TYPE='COMPLEX' TYPE_COMMENT='\!_COMPLEX' \
  ./lapw1c FORT=ifort FFLAGS=' -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML'
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw1'
modules.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp_.F
mv modules_tmp_.o modules.o
rm modules_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar.f
calkpt.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c calkpt_tmp_.F
mv calkpt_tmp_.o calkpt.o
rm calkpt_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cbcomb.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c coors.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c forfhs.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gaunt1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gaunt2.f
gtfnam.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam_tmp_.F
mv gtfnam_tmp_.o gtfnam.o
rm gtfnam_tmp_.F
hamilt.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hamilt_tmp_.F
mv hamilt_tmp_.o hamilt.o
rm hamilt_tmp_.F
hns.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hns_tmp_.F
mv hns_tmp_.o hns.o
rm hns_tmp_.F
horb.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c horb_tmp_.F
mv horb_tmp_.o horb.o
rm horb_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inikpt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inilpw.f
lapw1.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw1_tmp_.F
mv lapw1_tmp_.o lapw1.o
rm lapw1_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lohns.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lopw.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c nn.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
prtkpt.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c prtkpt_tmp_.F
mv prtkpt_tmp_.o prtkpt.o
rm prtkpt_tmp_.F
prtres.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c prtres_tmp_.F
mv prtres_tmp_.o prtres.o
rm prtres_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rdswar.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
seclit.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seclit_tmp_.F
mv seclit_tmp_.o seclit.o
rm seclit_tmp_.F
seclr4.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seclr4_tmp_.F
mv seclr4_tmp_.o seclr4.o
rm seclr4_tmp_.F
seclr5.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seclr5_tmp_.F
mv seclr5_tmp_.o seclr5.o
rm seclr5_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setkpt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setwar.f
tapewf.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c tapewf_tmp_.F
mv tapewf_tmp_.o tapewf.o
rm tapewf_tmp_.F
ustphx.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ustphx_tmp_.F
mv ustphx_tmp_.o ustphx.o
rm ustphx_tmp_.F
wfpnt.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wfpnt_tmp_.F
mv wfpnt_tmp_.o wfpnt.o
rm wfpnt_tmp_.F
wfpnt1.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wfpnt1_tmp_.F
mv wfpnt1_tmp_.o wfpnt1.o
rm wfpnt1_tmp_.F
ifort  -o ./lapw1c abc.o atpar.o bandv1.o calkpt.o cbcomb.o coors.o cputim.o dblr2k.o dgeqrl.o dgewy.o dgewyg.o dlbrfg.o dsbein1.o dscgst.o dstebz2.o dsyevx2.o dsymm2.o dsyr2m.o dsyrb4.o dsyrb5l.o dsyrdt4.o dsytrd2.o dsywyv.o dsyxev4.o dvbes1.o eisps.o errclr.o errflg.o forfhs.o gaunt1.o gaunt2.o gbass.o gtfnam.o hamilt.o hns.o horb.o inikpt.o inilpw.o lapw1.o latgen.o lmsort.o locdef.o lohns.o lopw.o matmm.o modules.o nn.o outerr.o outwin.o pdsyevx16.o prtkpt.o prtres.o pzheevx16.o rdswar.o rint13.o rotate.o rotdef.o seclit.o seclr4.o seclr5.o select.o service.o setkpt.o setwar.o sphbes.o stern.o tapewf.o ustphx.o vectf.o warpin.o wfpnt.o wfpnt1.o ylm.o zhcgst.o zheevx2.o zhemm2.o zher2m.o zhetrd2.o pdsyr2m.o pzher2m.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw1'
rm param.inc
