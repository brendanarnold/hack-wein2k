!-------------------------------------------------------------------
!      revised PBE exchange of Hammer, Hansen, and Norskov
!      (PRB 59, Vol.11 (1999))
!
! Input:  up: spin-up density
!         dn: spin-down density
!         grup: | grad up |
!         grdn: | grad dn |
!
! Output: exrpbe:  RPBE exchange energy per particle
!
!-------------------------------------------------------------------
      subroutine exrevpbe(up,dn,grup,grdn,exrpbe)
      implicit real*8 (a-h,o-z)

      exrup = 0.0d0
      exrdn = 0.0d0

      if (up.gt.1.0d-18) then
         up2 = 2.0d0*up
         grup2 = 2.0d0*grup
         call exrrpbe(up2,grup2,exrup)
         exrup = up2*exrup/(up+dn)
      endif

      if (dn.gt.1.0d-18) then 
         dn2 = 2.0d0*dn
         grdn2 = 2.0d0*grdn
         call exrrpbe(dn2,grdn2,exrdn)
         exrdn = dn2*exrdn/(up+dn)
      endif

      exrpbe = 0.5d0*(exrup + exrdn)

      return
      end

!-----------------------------------------------------------------
      subroutine exrrpbe(rho,grad,exrr)
      implicit real*8 (a-h,o-z)
      parameter(THRD=1.d0/3.d0)
      parameter(PI=3.14159265358979323846264338327950d0)
      parameter(AX=-0.738558766382022405884230032680836d0)
      parameter(UM=0.2195149727645171d0,UK=0.8040d0,ul=um/uk)
      parameter (FAC = 3.09366772628013593097D0)
!     FAC = (3*PI**2)**(1/3)
      if (rho.gt.1.0d-18) then
         fer = FAC*(rho**THRD)
         denom = 2.0d0*fer*rho
         ss = dabs(grad)/denom
         exunif = AX*rho**THRD
         s2 = ss*ss
         fxrpbe = 1.0d0 + UK*(1.0d0 - dexp(-UM*s2/UK) )
         exrr = exunif*fxrpbe
      else
         exrr = 0.0d0
      endif

      return
      end

!-------------------------------------------------------------------
