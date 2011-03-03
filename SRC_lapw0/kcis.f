!-------------------------------------------------------------------
!     correlation energy functional of Krieger, Chen, Iafrate, 
!                                      and Savin (KCIS)
!
! Reference: J.B. Krieger, J. Chen, G.I. Iafrate, A. Savin
!            ``Construction of an accurate self-interaction-
!              corrected correlation energy functional based on 
!              an electron gas with a gap'', preprint
!
! Input:  up: spin-up density
!         dn: spin-down density
!         grup: | grad up |
!         grdn: | grad dn |
!         tauup: 0.5*sum_{i,occ} | grad psi_{i,up} |^2
!                kinetic energy density of up-electrons
!         taudn: 0.5*sum_{i,occ} | grad psi_{i,down} |^2
!                kinetic energy density of down-electrons
! Output: ecsic: KCIS correlation energy per particle
!
! Requires subroutine corlsd to calculate the correlation energy 
! per particle of the uniform electron gas using the parametrization 
! of Perdew and Wang 
!
!-------------------------------------------------------------------
      subroutine kciscor(up,dn,grup,grdn,gr,tauup,taudn,ecsic)
     
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
       if(tauw.gt.tau.and.up.lt.1.d-1) write(6,*) 'tauw.gt.tau',tauw,tau
             if(tauwup.gt.tauup) tauup=tauwup
             if(tauwdn.gt.taudn) taudn=tauwdn

      dens = up + dn 
      if (dens.gt.1.0d-18) then 
         call ggagap(up,dn,gr,ecgapgga)
         ecsic = ecgapgga

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

         delc = 0.0d0
!        spin up correction
         if (up.gt.1.0d-18) then
            dd = ZERO
            call ggagap(up,dd,grup,ecgapup)
            tfac = tauwup/tauup
            delc = delc + tfac*up*ecgapup
         endif

!        spin down correction
         if (dn.gt.1.0d-18) then
            uu = ZERO
            call ggagap(uu,dn,grdn,ecgapdn)
            tfac = tauwdn/taudn
            delc = delc + tfac*dn*ecgapdn
         endif
         ecsic = ecsic - delc/dens
      else 
         ecsic = 0.0d0
      endif
      return
      end

!-------------------------------------------------------------------

      subroutine ggagap(up,dn,grad,ecggagap)
      implicit real*8 (a-h,o-z)

      parameter ( PI = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( TWOTHRD = 2.0d0*THRD )
      parameter ( FRTHRD = 4.0d0*THRD )
      parameter ( FAC = 3.09366772628013593097d0 )
!     FAC = (3*PI**2)**(1/3)
      parameter ( FAC2 = 0.519842100d0)
!     FAC2 = 2*( 2^(1/3) - 1 )
      parameter ( ALPHA = 1.91915829267751300662482032624669d0 )
      parameter ( BETA = 0.066725d0 )

      gap = 0.0d0
      dens = up + dn
      if (dens.gt.1.0d-18) then
         zeta = (up-dn)/dens
         fer = FAC*dens**THRD
         rs = ALPHA/fer
         sk = dsqrt(4.0d0*fer/PI)

         gap = 0.125*grad*grad/(dens*dens)
       
         call ccunp(rs,cc1,cc2,cc3,eclda,eclsd)
 
         fzet = (1.0d0 + zeta)**FRTHRD + (1.0d0 - zeta)**FRTHRD
         fzet = (fzet - 2.0d0)/FAC2

         ecgg1 = 0.0d0
         ecgg2 = 0.0d0

         gg = 1.0d0
         twoksg = 2.0d0*sk*gg
         tt = dabs(grad)/(twoksg*dens)
         fgga = 1.0d0 + tt*tt*(gg**3.0d0)/dabs(eclda)
         fgga = BETA*dlog(fgga)
         ecgg = eclda/(1.0d0 + fgga)
         ecgg1 = (ecgg + cc1*gap)/(1.0d0 + cc2*gap &
                  + cc3*gap*gap)
         gg = 0.5d0**THRD
         twoksg = 2.0d0*sk*gg
         tt = dabs(grad)/(twoksg*dens)
         fgga = 1.0d0 + tt*tt*(gg**3.0d0)/dabs(eclsd)
         fgga = BETA*dlog(fgga)
         ecgg = eclsd/(1.0d0 + fgga)
         ecgg2 = (ecgg + 0.70d0*cc1*gap)/(1.0d0 &
                  + 1.5d0*cc2*gap + 2.59d0*cc3*gap*gap)
         ecggagap = ecgg1 + fzet*(ecgg2 - ecgg1)

      else
         ecggagap = 0.0d0
      endif
      
      return
      end

!-------------------------------------------------------------------

      subroutine ccunp(rs,cc1unp,cc2unp,cc3unp,eclda,eclsd)
      implicit real*8 (a-h,o-z)

      parameter ( A1 = 0.04953D0)
      parameter ( A2 = 1.07924D0)
      parameter ( A3 = 0.07928D0)

      parameter ( B3 = -0.02504D0)
      parameter ( B4 = 0.007026D0)
      parameter ( B5 = -0.001268D0)
      parameter ( B6 = 0.0001136D0)
      parameter ( B7 = -0.000003841D0)

      parameter ( ALPHA = 1.91915829267751300662482032624669d0 )
      parameter ( CCFAC = 0.06483D0)

      rsth = rs**1.5d0
      epsp = 1.0d0 + A2*dsqrt(rs) + A3*rs + A1*rsth
      epsp = A1*rsth/epsp

      epspp = B3*rs**3.0d0 + B4*rs**4.0d0 + B5*rs**5d0
      epspp = epspp + B6*rs**6.0d0 + B7*rs**7.0d0
    
      fer = ALPHA/rs
      cceps = CCFAC*fer*fer

      zeta = 0.0d0
      call corlsd(rs,zeta,eclda,vcup,vcdn,ecrs,eczet,alfc)

      denom = 2.0d0*(cceps*epsp - eclda*eclda)
      cc3unp = - (2.0d0*epsp*epsp - eclda*epspp)/denom
      cc2unp = (2.0d0*eclda*epsp - cceps*epspp)/denom
      cc1unp = - cceps*cc3unp

      zeta = 1.0d0
      call corlsd(rs,zeta,eclsd,vcup,vcdn,ecrs,eczet,alfc)
      return 
      end

!-------------------------------------------------------------------
!     end of subroutines for KCIS correlation functional 
!-------------------------------------------------------------------
