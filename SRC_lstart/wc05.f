      subroutine WC05(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap, &
                      exupwc05,exdnwc05,vxupwc05,vxdnwc05)

      implicit real*8(a-h,o-z)
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)

      rho2=2.d0*up
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrup/(2.d0*fk*rho2)
        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
        call exchwc05(rho2,s,u,v,exupwc05,vxupwc05)
      else
        exupwc05=0.d0
        vxupwc05=0.d0
      endif
      
      rho2=2.d0*dn
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrdn/(2.d0*fk*rho2)
        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
        call exchwc05(rho2,s,u,v,exdnwc05,vxdnwc05)
      else
        exdnwc05=0.d0
        vxdnwc05=0.d0
      endif

      return
      end

!-----------------------------------------------------------------------

      SUBROUTINE EXCHWC05(rho,S,U,V,EX,VX)

      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(CKAPPA=0.804D0, CMU=0.21951D0, CC=0.0079325D0)
      
! construct LDA exchange energy density
      exunif = AX*rho**THRD

! construct WC05 enhancement factor
      S2 = S*S

      XWC05 = 10D0/81D0*S2 + (CMU - 10D0/81D0)*S2*DEXP(-S2) + &
              DLOG(1D0 + CC*S2*S2)

      FxWC05 = 1D0 + CKAPPA - CKAPPA/(1D0 + XWC05/CKAPPA)

      EX = exunif*FxWC05

!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxWC05/ ds
!  Fss=d Fs/ds
      YWC05 = 20D0/81D0 + (2D0*CMU - 20D0/81D0)*(1D0 - S2)*DEXP(-S2) +&
              4D0*CC*S2/(1D0 + CC*S2*S2)

      FS = YWC05/(1D0 + XWC05/CKAPPA)**2D0

      DYWC05 = (4D0*CMU - 40D0/81D0)*(S2*S - 2D0*S)*DEXP(-S2) + &
               8D0*CC*S*(1D0 - CC*S2*S2)/(1D0 + CC*S2*S2)**2D0

      FSS = DYWC05/(1D0 + XWC05/CKAPPA)**2D0 - &
            2D0/CKAPPA*S*YWC05**2D0/(1D0 + XWC05/CKAPPA)**3D0

      VX = exunif*(THRD4*FxWC05-(U-THRD4*S2*S)*FSS-V*FS)

      RETURN
      END
