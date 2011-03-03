      subroutine sicpbe_tpss(up,dn,grup,grdn,gr,tauup,taudn,grupdn,ectpss)

      implicit real*8 (a-h,o-z)

!c grupdn = del n_up dot del n_down
!c ectpss = TPSS meta-GGA correlation energy per electron
!c grup = sqrt( grada(1)**2 + grada(2)**2 + grada(3)**2 )
!c grdn = sqrt( gradb(1)**2 + gradb(2)**2 + gradb(3)**2 )
!c gr   = sqrt( grad(1)**2 + grad(2)**2 + grad(3)**2 )
!c grupdn = grada(1)*gradb(1) + grada(2)*gradb(2) + grada(3)*gradb(3)

      parameter ( PI = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( TWOTHRD = 2.0d0*thrd)
      parameter ( FAC = 3.09366772628013593097d0 )
!c      FAC = (3*PI**2)**(1/3)
      parameter ( ALPHA = 1.91915829267751300662482032624669d0 )

      parameter ( GGPOL = 0.793700526d0)
!c     GGPOL = 0.5*(2**TWOTHRD)

      aa = 0.87d0
      bb = 0.5d0
      cc = 2.26d0
      dd = 2.8d0
      dnu = 1.d0

      dens = up + dn
      if (dens.gt.1.0d-18) then
         zeta = (up - dn)/dens
         CCC = 0.53d0 + aa*zeta**2.d0 + bb*zeta**4.d0 + cc*zeta**6.d0
         gg = 0.5d0*((1.0d0+zeta)**TWOTHRD + (1.0d0-zeta)**TWOTHRD)
         fer = FAC*dens**THRD
         rs = ALPHA/fer
         sk = dsqrt(4.0d0*fer/PI)
         twoksg = 2.0d0*sk*gg
         tt = gr/(twoksg*dens)
         call CORPBE_tpss(RS,ZETA,TT,ECC,HH,ECPBE)
         if (up.gt.1.0d-18) then
            tauwup = 0.125d0*grup*grup/up
         else
            tauwup = 0.0d0
         endif
         if (dn.gt.1.0d-18) then
            tauwdn = 0.125d0*grdn*grdn/dn
         else
            tauwdn = 0.0d0
         endif
         tauw = 0.125d0*gr**2.d0/dens
         tau = tauup + taudn
       if(tauw.gt.tau.and.up.lt.1.d-1) write(6,*) 'tauw.gt.tau',tauw,tau
             if(tauw.gt.tau) tau=tauw

         xx2 = (tauw/tau)**2.d0
         xx3 = (tauw/tau)**3.d0
!c----------------------------
         ter1 = (1.d0 - zeta)**2.d0*grup**2.d0
         ter2 = (1.d0 + zeta)**2.d0*grdn**2.d0
         ter3 = 2.d0*(1.d0 - zeta**2.d0)*grupdn
         sumz = ter1 + ter2 - ter3
         sumz = dabs(sumz)
         grz = dsqrt(sumz)/dens
         szeta2 = (grz/(2.d0*fer))**2.d0

         if (dabs(zeta).lt.0.9999d0) then
            term1 = 1.d0/(1.d0 + zeta)**(4.d0/3.d0)
            term2 = 1.d0/(1.d0 - zeta)**(4.d0/3.d0) 
            yyz = szeta2*(term1 + term2)/2.d0
         else
            yyz = 0.d0
         end if
         CCC = CCC/(1.d0 + dnu*yyz)**4.d0
         pbe = ecc + hh
         en = pbe*(1.d0 + CCC*xx2)

         delc = 0.0d0
!c        spin up correction
         if (up.gt.1.0d-18) then
            zeta = 1.0d0
            fer = FAC*up**THRD
            rs = ALPHA/fer
            sk = dsqrt(4.0d0*fer/PI)
            twoksg = 2.0d0*sk*GGPOL
            tt = dabs(grup)/(twoksg*up)
            call CORPBE_tpss(RS,ZETA,TT,ECC,HH,ECPBE)
            tfac = (tauw/tau)**2.d0
            pbeup = ecc + hh
            if(pbeup.lt.pbe) then
              delc = delc + tfac*up*pbe
            else
              delc = delc + tfac*up*pbeup
            end if
         endif

!c        spin down correction
         if (dn.gt.1.0d-18) then
            zeta = -1.0d0
            fer = FAC*dn**THRD
            rs = ALPHA/fer
            sk = dsqrt(4.0d0*fer/PI)
            twoksg = 2.0d0*sk*GGPOL
            tt = dabs(grdn)/(twoksg*dn)
            call CORPBE_tpss(RS,ZETA,TT,ECC,HH,ECPBE)
            tfac = (tauw/tau)**2.d0
            pbedn = ecc + hh
            if(pbedn.lt.pbe) then
              delc = delc + tfac*dn*pbe
            else
              delc = delc + tfac*dn*pbedn
            end if
         endif
         delc = -(1.d0 + CCC)*delc/dens
         revsic = en + delc
         ectpss = revsic*(1.d0 + dd*revsic*xx3)
      else
         ectpss = 0.0d0
      endif


      return
      end

