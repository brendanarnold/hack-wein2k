rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod
if [ -f .parallel ]; then \
   rm -f .parallel modules.o reallocate.o energy.o gtfnam.o lapw0.o  *.mod; \
fi
touch .sequential
make lapw0 FORT=ifort FFLAGS=' -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML'
make[1]: Entering directory `/usr/local/WIEN2K_07_1/SRC_lapw0'
icc -c cputim.c
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c modules.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c reallocate.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ainv.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c blyp.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charg2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charg3.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c charge.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c chfac.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c chslv.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c corgga.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c cub_xc_back.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c corlsd.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dlu.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c drho.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dfxhpbe.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dylm.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c efg.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c energy.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errclr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c errflg.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ev92.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ev92ex.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c exch.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c exch17.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c exrevpbe.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fithi.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fxhpbe.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gbass.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gcor.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gea.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c geaex.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c getfft.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c getff1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gpoint.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c grans.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c gtfnam.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c hcth.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ifflim.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c kcis.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c lapw0.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c latgen.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c multfc.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c multsu.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c outerr.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pbe1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c poissn.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c potfac.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pwxad4.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c pwxad5.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c qranf.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c readstruct.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rean0.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rean1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rean3.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rean4.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotate.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rotdef.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setff0.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setff1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setfft.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c setff2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c seval.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sevald.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sevaldd.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sevali.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sevalin.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sicpbe.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sphbes.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c spline.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c srolyl.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c stern.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sumfac.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c suml.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c th1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c th2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vpw91.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vresp.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vs98.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxc15.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxc16.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxc17.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxc24.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxc26.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxclm2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxcpw2.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxi35.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c vxi35a.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c workf1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xcener.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xcpot.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xcpot1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c xcpot3.F
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ykav.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c ylm.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c zfft3d.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c sicpbe_tpss.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c corpbe_tpss.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c fx_tpss.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c dfxtpss.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rhopw.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c b88.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c wc05.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c eramps.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c epot1.f
ifort -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -c rpbe.f
ifort -o lapw0 cputim.o modules.o reallocate.o ainv.o blyp.o charg2.o  charg3.o charge.o chfac.o chslv.o  corgga.o  cub_xc_back.o corlsd.o dlu.o drho.o  dfxhpbe.o dylm.o  efg.o  energy.o errclr.o errflg.o ev92.o    ev92ex.o exch.o exch17.o exrevpbe.o  fithi.o fxhpbe.o gbass.o gcor.o  gea.o geaex.o  getfft.o  getff1.o gpoint.o  grans.o gtfnam.o hcth.o ifflim.o kcis.o lapw0.o  latgen.o  multfc.o multsu.o outerr.o  pbe1.o poissn.o potfac.o pwxad4.o pwxad5.o  qranf.o readstruct.o rean0.o rean1.o rean3.o  rean4.o  rotate.o  rotdef.o setff0.o  setff1.o  setfft.o setff2.o seval.o sevald.o sevaldd.o sevali.o sevalin.o sicpbe.o sphbes.o  spline.o srolyl.o stern.o sumfac.o  suml.o th1.o th2.o vpw91.o   vresp.o vs98.o vxc15.o vxc16.o vxc17.o vxc24.o  vxc26.o  vxclm2.o vxcpw2.o vxi35.o vxi35a.o workf1.o xcener.o  xcpot.o  xcpot1.o xcpot3.o ykav.o  ylm.o zfft3d.o sicpbe_tpss.o corpbe_tpss.o fx_tpss.o dfxtpss.o rhopw.o b88.o wc05.o eramps.o epot1.o rpbe.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -L/opt/intel/mkl/9.0/lib/em64t -lpthread 
make[1]: Leaving directory `/usr/local/WIEN2K_07_1/SRC_lapw0'
make: *** No rule to make target `complex'.  Stop.
