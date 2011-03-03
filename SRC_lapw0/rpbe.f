      SUBROUTINE RPBE(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap, &
                      exup,exdn,vxup,vxdn)

!B. Hammer, L. B. Hansen, and J. K. Norskov, Phys. Rev. B 59, 7413 (1999).

      implicit real*8(a-h,o-z)
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)

      rho2=2.d0*up
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrup/(2.d0*fk*rho2)
        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
        call exchrpbe(rho2,s,u,v,exup,vxup)
      else
        exup=0.d0
        vxup=0.d0
      endif
      
      rho2=2.d0*dn
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrdn/(2.d0*fk*rho2)
        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
        call exchrpbe(rho2,s,u,v,exdn,vxdn)
      else
        exdn=0.d0
        vxdn=0.d0
      endif

      return
      end

!-----------------------------------------------------------------------

      SUBROUTINE EXCHRPBE(rho,S,U,V,EX,VX)

      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(CKAPPA=0.804D0,CMU=0.21951D0)
      
! construct LDA exchange energy density
      exunif = AX*rho**THRD

! construct RPBE enhancement factor
      S2 = S*S

      Fx = 1D0 + CKAPPA*(1D0 - DEXP(-(CMU/CKAPPA)*S2))

      EX = exunif*Fx

!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d Fx/ ds
!  Fss=d Fs/ds
      FS = 2D0*CMU*DEXP(-(CMU/CKAPPA)*S2)

      FSS = -4D0*(CMU*CMU/CKAPPA)*S*DEXP(-(CMU/CKAPPA)*S2)

      VX = exunif*(THRD4*Fx-(U-THRD4*S2*S)*FSS-V*FS)

      RETURN
      END
