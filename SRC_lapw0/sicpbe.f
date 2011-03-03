      subroutine sicpbe(up,dn,grup,grdn,gr,tauup,taudn,ecsicpbe)

!     up: spin-up density
!     dn: spin-dn density
!     grup: grad up   grdn: grad dn 
!     gr: grad (up+dn)
!     tauup: kinetic energy density for spin up
!     taudn: kinetic energy density for spin dn
!     ecsicpbe: correlation energy per particle of the self-interaction 
!               corrected PBE correlation energy (or hyper-PBE)
     
      implicit real*8 (a-h,o-z)

      parameter ( PI = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( TWOTHRD = 2.0d0*thrd)
      parameter ( FAC = 3.09366772628013593097d0 )
!     FAC = (3*PI**2)**(1/3)
      parameter ( ALPHA = 1.91915829267751300662482032624669d0 )
      parameter ( ZERO = 0.0d0 )
      parameter ( GGPOL = 0.793700526d0)
!     GGPOL = 0.5*(2**TWOTHRD)

      parameter ( CCC = 0.53d0 )

      dens = up + dn 
      if (dens.gt.1.0d-18) then 
         zeta = (up - dn)/dens
         gg = 0.5d0*((1.0d0+zeta)**TWOTHRD + (1.0d0-zeta)**TWOTHRD)
         fer = FAC*dens**THRD
         rs = ALPHA/fer
         sk = dsqrt(4.0d0*fer/PI)
         twoksg = 2.0d0*sk*gg
         tt = gr/(twoksg*dens)
         uu = ZERO
         vv = ZERO
         ww = ZERO
         lgga = 1
         lpot = 0 
         call corpbe(rs,zeta,tt,uu,vv,ww,lgga,lpot,ecc,vcupl,vcdnl, &
                     hh,dvcup,dvcdn)
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
         tauw = tauwup + tauwdn
         tau = tauup + taudn
!             if(tauw.lt.1.d-7.or.tau.lt.1.d-7) then
!                tauw=1.d-7
!                tau=1.d-7
!             write(6,*) 'tauw.lt.0',tauw,tau
!             endif
       if(tauw.gt.tau.and.up.lt.1.d-1) write(6,*) 'tauw.gt.tau',tauw,tau
             if(tauw.gt.tau) tau=tauw
         en = (ecc + hh)*(1.0d0 + ccc*(tauw/tau)**2.0d0)
             enpbe=ecc+hh
             efak1=(1.0d0 + ccc*(tauw/tau)**2.0d0)
         delc = 0.0d0
!        spin up correction
         if (up.gt.1.0d-18) then
            zeta = 1.0d0
            fer = FAC*up**THRD
            rs = ALPHA/fer
            sk = dsqrt(4.0d0*fer/PI)
            twoksg = 2.0d0*sk*GGPOL
            tt = dabs(grup)/(twoksg*up)
            uu = ZERO 
            vv = ZERO
            ww = ZERO
            lgga = 1
            lpot = 0
            call corpbe(rs,zeta,tt,uu,vv,ww,lgga,lpot,ecc,vcupl,vcdnl, &
                        hh,dvcup,dvcdn)
!             if(tauwup.lt.1.d-7.or.tauup.lt.1.d-7) then
!                tauwup=1.d-7
!                tauup=1.d-7
!             endif
             if(tauwup.gt.tauup) tauup=tauwup

            tfac = (tauwup/tauup)**2.0d0
            delc = delc + tfac*up*(ecc + hh)
         endif

!        spin down correction
         if (dn.gt.1.0d-18) then
            zeta = -1.0d0
            fer = FAC*dn**THRD
            rs = ALPHA/fer
            sk = dsqrt(4.0d0*fer/PI)
            twoksg = 2.0d0*sk*GGPOL
            tt = dabs(grdn)/(twoksg*dn)
            uu = ZERO 
            vv = ZERO
            ww = ZERO
            lgga = 1
            lpot = 0
            call corpbe(rs,zeta,tt,uu,vv,ww,lgga,lpot,ecc,vcupl,vcdnl, &
                        hh,dvcup,dvcdn)
!             if(tauwdn.lt.1.d-7.or.taudn.lt.1.d-7) then
!                tauwdn=1.d-7
!                taudn=1.d-7
!             endif
             if(tauwdn.gt.taudn) taudn=tauwdn
            tfac = (tauwdn/taudn)**2.0d0
            delc = delc + tfac*dn*(ecc + hh)
         endif
         delc = -(1.0d0 + CCC)*delc/dens
         ecsicpbe = en + delc
!      if(tauw.gt.tau.and.up.lt.1.d-1)
!     *    write(6,*) 'tauw.gt.tau',tauw,tau,up,ecsicpbe,enpbe,en,delc
      else 
         ecsicpbe = 0.0d0
      endif
      return
      end

!-------------------------------------------------------------------
