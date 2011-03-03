rm -f *.o *_tmp_.* *.P .real .complex .sequential .parallel *.mod
if [ -f .complex ]; then \
   rm -f .complex modules.o reallocate.o essl.o fourir.o fsumai1.o l2main.o lapw2.o read_vec.o setfft1.o setfft2.o  *.mod; \
fi
touch .real
if [ -f .parallel ]; then \
   rm -f .parallel reallocate.o atpar.o fourir.o gtfnam.o l2main.o lapw2.o modules.o read_vec.o *.mod; \
fi
touch .sequential
make TYPE='REAL' TYPE_COMMENT='\!_REAL' \
  ./lapw2 FORT=ifort FFLAGS=' -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML'
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw2'
icc -c cputim.c
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
modules.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp_.F
mv modules_tmp_.o modules.o
rm modules_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
atpar.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar_tmp_.F
mv atpar_tmp_.o atpar.o
rm atpar_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c c3fft.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c csplit.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dfrad.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c d5splt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dergl.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c deter.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c diracout.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dvbes1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
essl.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c essl_tmp_.F
mv essl_tmp_.o essl.o
rm essl_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c f7splt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fermi.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fermi5.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fomai1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fomai2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fomai3.f
fourir.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fourir_tmp_.F
mv fourir_tmp_.o fourir.o
rm fourir_tmp_.F
fsumai1.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fsumai1_tmp_.F
mv fsumai1_tmp_.o fsumai1.o
rm fsumai1_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gaunt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c getfft.f
gtfnam.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam_tmp_.F
mv gtfnam_tmp_.o gtfnam.o
rm gtfnam_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c harmon.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ifflim.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inouh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inth.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kdelta.f
l2main.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c l2main_tmp_.F
mv l2main_tmp_.o l2main.o
rm l2main_tmp_.F
lapw2.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw2_tmp_.F
mv lapw2_tmp_.o lapw2.o
rm lapw2_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lomain.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mag.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c mknam.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c notri.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c open_helpfile.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outp.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c p3splt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ph.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c psplit.f
read_vec.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c read_vec_tmp_.F
mv read_vec_tmp_.o read_vec.o
rm read_vec_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c recfil.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c recpr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
setfft1.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setfft1_tmp_.F
mv setfft1_tmp_.o setfft1.o
rm setfft1_tmp_.F
setfft2.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setfft2_tmp_.F
mv setfft2_tmp_.o setfft2.o
rm setfft2_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spline.f
stern.F: REAL version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern_tmp_.F
mv stern_tmp_.o stern.o
rm stern_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumupfft.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c t3j.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c t3j0.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vdrho.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xsplt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sortag.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c efgsplit.f
ifort  -o ./lapw2 cputim.o reallocate.o modules.o abc.o  atpar.o c3fft.o charge.o csplit.o dfrad.o d5splt.o dergl.o deter.o diracout.o dvbes1.o errclr.o errflg.o essl.o f7splt.o fermi.o fermi5.o fomai1.o fomai2.o fomai3.o fourir.o fsumai1.o gaunt.o getfft.o gtfnam.o harmon.o ifflim.o inouh.o inth.o kdelta.o l2main.o lapw2.o latgen.o lomain.o mag.o mknam.o notri.o open_helpfile.o outerr.o outp.o outwin.o p3splt.o ph.o psplit.o read_vec.o recfil.o recpr.o rint13.o rotate.o rotdef.o setfft1.o setfft2.o sphbes.o spline.o stern.o sumupfft.o t3j.o t3j0.o vdrho.o xsplt.o ylm.o  sortag.o efgsplit.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw2'
if [ -f .real ]; then \
   rm -f .real modules.o reallocate.o essl.o fourir.o fsumai1.o l2main.o lapw2.o read_vec.o setfft1.o setfft2.o  *.mod; \
fi
touch .complex
if [ -f .parallel ]; then \
   rm -f .parallel reallocate.o atpar.o fourir.o gtfnam.o l2main.o lapw2.o modules.o read_vec.o *.mod; \
fi
touch .sequential
make TYPE='COMPLEX' TYPE_COMMENT='\!_COMPLEX' \
  ./lapw2c FORT=ifort FFLAGS=' -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML'
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw2'
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
modules.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules_tmp_.F
mv modules_tmp_.o modules.o
rm modules_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c abc.f
atpar.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c atpar_tmp_.F
mv atpar_tmp_.o atpar.o
rm atpar_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c csplit.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c d5splt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dergl.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c diracout.f
essl.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c essl_tmp_.F
mv essl_tmp_.o essl.o
rm essl_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c f7splt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fermi.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fermi5.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fomai1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fomai2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fomai3.f
fourir.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fourir_tmp_.F
mv fourir_tmp_.o fourir.o
rm fourir_tmp_.F
fsumai1.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fsumai1_tmp_.F
mv fsumai1_tmp_.o fsumai1.o
rm fsumai1_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c getfft.f
gtfnam.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam_tmp_.F
mv gtfnam_tmp_.o gtfnam.o
rm gtfnam_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c harmon.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inouh.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c inth.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kdelta.f
l2main.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c l2main_tmp_.F
mv l2main_tmp_.o l2main.o
rm l2main_tmp_.F
lapw2.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw2_tmp_.F
mv lapw2_tmp_.o lapw2.o
rm lapw2_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lomain.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outp.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outwin.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c p3splt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c psplit.f
read_vec.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c read_vec_tmp_.F
mv read_vec_tmp_.o read_vec.o
rm read_vec_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c recfil.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c recpr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rint13.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
setfft1.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setfft1_tmp_.F
mv setfft1_tmp_.o setfft1.o
rm setfft1_tmp_.F
setfft2.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setfft2_tmp_.F
mv setfft2_tmp_.o setfft2.o
rm setfft2_tmp_.F
stern.F: COMPLEX version extracted
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern_tmp_.F
mv stern_tmp_.o stern.o
rm stern_tmp_.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vdrho.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xsplt.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort  -o ./lapw2c cputim.o reallocate.o modules.o abc.o  atpar.o c3fft.o charge.o csplit.o dfrad.o d5splt.o dergl.o deter.o diracout.o dvbes1.o errclr.o errflg.o essl.o f7splt.o fermi.o fermi5.o fomai1.o fomai2.o fomai3.o fourir.o fsumai1.o gaunt.o getfft.o gtfnam.o harmon.o ifflim.o inouh.o inth.o kdelta.o l2main.o lapw2.o latgen.o lomain.o mag.o mknam.o notri.o open_helpfile.o outerr.o outp.o outwin.o p3splt.o ph.o psplit.o read_vec.o recfil.o recpr.o rint13.o rotate.o rotdef.o setfft1.o setfft2.o sphbes.o spline.o stern.o sumupfft.o t3j.o t3j0.o vdrho.o xsplt.o ylm.o  sortag.o efgsplit.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread -lmkl_lapack64 -lmkl_em64t -lguide -lvml -pthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw2'
