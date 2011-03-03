!-------------------------------------------------------------------
!     exchange-correlation functional of Becke, Lee, Yang, and Parr
!
! Reference:  Exchange: A. Becke, PRA 38, p.3098 (1988)
!             Correlation: C. Lee, W. Yang, R.G. Parr, PRB 37, p.785 
!                          (1988); the implemented formula is obtained 
!                          by partial integration from the original 
!                          work and is taken from B. Miehlich, 
!                          A. Savin, H. Stoll, and H. Preuss, 
!                          Chem. Phys. Lett. 157, p.200 (1989)
!
! Input:  up: spin-up density
!         dn: spin-down density
!         grup: | grad up |
!         grdn: | grad dn |
!         gr  : | grad (up + dn) |
! Output: exb88: Becke 88 exchange energy per particle
!         eclyp: Lee, Yang, Parr 88 correlation energy per particle
!
!-------------------------------------------------------------------
      subroutine excblyp(up,dn,grup,grdn,gr,exb88,eclyp)
      implicit real*8 (a-h,o-z)

      exb88 = 0.0d0
      eclyp = 0.0d0

      exup = 0.0d0
      exdn =0.0d0

!     spin-up exchange energy

      if (up.gt.1.0d-18) then
         call exbe88(up,grup,exup)
         exup = exup*up/(up+dn)
      endif

!     spin-down exchange energy

      if (dn.gt.1.0d-18) then
         call exbe88(dn,grdn,exdn)
         exdn = exdn*dn/(up+dn)
      endif

      exb88 = exup + exdn
      
      if ((up+dn).gt.1.0d-18) then
         call corlyp(up,dn,grup,grdn,gr,eclyp)
      endif

      return
      end
!-------------------------------------------------------------------
      subroutine exbe88(rho,gr,exbec)

      implicit real*8 (a-h,o-z)

      parameter ( PI = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( FOURTHRD = 4.0d0*thrd)
      parameter ( CX = 0.930525736D0 )
      parameter ( BETA = 0.0042D0 )

      xx = dabs(gr)/(rho**FOURTHRD)
      arsi = dlog( xx + dsqrt(xx*xx + 1.0d0) )
      denom = CX*(1.0d0 + 6.0d0*BETA*xx*arsi)
      totfac = 1.0d0 + BETA*xx*xx/denom
      exbec = -CX*totfac*(rho**THRD)
      
      return
      end

!-------------------------------------------------------------------

      subroutine corlyp(up,dn,grup,grdn,gr,eclyp)
      implicit real*8 (a-h,o-z)

      parameter ( PI = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( FOURTHRD = 4.0d0*thrd)
      parameter ( FAC = 3.09366772628013593097d0 )
!     FAC = (3*PI**2)**(1/3)
      parameter ( CF = 2.871234000d0 )
!     CF = 0.3*FAC*FAC

      parameter ( AA = 0.04918d0 )
      parameter ( BB = 0.132d0 )
      parameter ( CC = 0.2533d0 )
      parameter ( DD = 0.349d0 )

      rho = up + dn
      
      if (rho.gt.1.0d-18) then 
         rhmth = rho**(-THRD)
         denom = 1.0d0 + DD*rhmth
         om = dexp( - CC*rhmth)/denom
         om = om*(rhmth**11.0d0)
         del = CC*rhmth + DD*rhmth/denom

         ec1lyp = - AA*4.0d0*up*dn/(denom*rho)

         brac1 = up**(8.0d0*THRD) + dn**(8.0d0*THRD)
         brac1 = 2.0d0**(11.0d0*THRD)*CF*brac1 
         brac1 = brac1 + (47.0d0 - 7.0d0*del)*gr*gr/18.0d0
         ter1 = (2.5d0 - del/18.0d0)*(grup*grup + grdn*grdn)
         brac1 = brac1 - ter1 
         ter1 = (del - 11.0d0)/9.0d0
         ter1 = ter1*(up*grup*grup + dn*grdn*grdn)/rho
         brac1 = brac1 - ter1
         brac1 = up*dn*brac1

         brac2 = - 2.0d0*THRD*(rho*gr)**2.0d0
         brac2 = brac2 + (2.0d0*THRD*rho*rho - up*up)*grdn*grdn
         brac2 = brac2 + (2.0d0*THRD*rho*rho - dn*dn)*grup*grup

         ec2lyp = - AA*BB*om*(brac1 + brac2)

         eclyp = (ec1lyp + ec2lyp)/rho
      else
         eclyp = 0.0d0 
      endif

      return 
      end

!-------------------------------------------------------------------

