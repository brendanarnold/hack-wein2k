      subroutine B88(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap, &
                     exupb88,exdnb88,vxupb88,vxdnb88)

      implicit real*8(a-h,o-z)
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)

      rho2=2.d0*up
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrup/(2.d0*fk*rho2)
        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
        call exchb88(rho2,s,u,v,exupb88,vxupb88)
      else
        exupb88=0.d0
        vxupb88=0.d0
      endif
      
      rho2=2.d0*dn
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrdn/(2.d0*fk*rho2)
        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
        call exchb88(rho2,s,u,v,exdnb88,vxdnb88)
      else
        exdnb88=0.d0
        vxdnb88=0.d0
      endif

      return
      end

!-----------------------------------------------------------------------

      SUBROUTINE EXCHB88(rho,S,U,V,EX,VX)

      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      
! construct LDA exchange energy density
      exunif = AX*rho**THRD

! construct B88 enhancement factor
      ALPHAB88 = 0.00451357747125D0
      BETAB88 = 0.0252D0
      BB88 = 2D0*(6D0*PI**2D0)**(1D0/3D0)
      C1B88 = ALPHAB88*BB88**2D0
      C2B88 = BETAB88*BB88
      S2 = S*S
      P1B88 = C1B88*S2
      P2B88 = 1D0 + C2B88*S*DLOG(BB88*S + DSQRT(1D0 + BB88**2D0*S2))
      FxB88 = 1D0 + P1B88/P2B88
      EX = exunif*FxB88

!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxB88/ ds
!  Fss=d Fs/ds
      IF (S .GT. 1D-18) THEN

        DP1B88 = 2D0*C1B88*S
        D2P1B88 = 2D0*C1B88
        DP2B88 = C2B88*DLOG(BB88*S + DSQRT(1D0 + BB88**2D0*S2)) + &
                 C2B88*S*BB88/DSQRT(1D0 + BB88**2D0*S2)
        D2P2B88 = 2D0*C2B88*BB88/DSQRT(1D0 + BB88**2D0*S2) - &
                  C2B88*BB88**3D0*S2/(1D0 + BB88**2D0*S2)**1.5D0
        Fs = (1D0/S)*(DP1B88/P2B88 - P1B88*DP2B88/P2B88**2D0)
        Fss = -Fs/S + (1D0/S)*(D2P1B88/P2B88 - &
              2D0*DP1B88*DP2B88/P2B88**2D0 - P1B88*D2P2B88/P2B88**2D0+&
              2D0*P1B88*DP2B88**2D0/P2B88**3D0)

      ELSE

        Fs = 0D0
        Fss = 0D0

      ENDIF

      VX = exunif*(THRD4*FxB88-(U-THRD4*S2*s)*FSS-V*FS)

      RETURN
      END
