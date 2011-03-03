      SUBROUTINE COUPLX(l,theta,phi,ipr)
!.....THIS IS THE ROUTINE CALCULATES THE MATRIX ELEMENT OF
!     <Ylm' ms'|X|Ylm ms>. 
!     in this version X=l(dzeta) i.e. component of the angular
!     momentum in direction dzeta = [sin(theta)cos(phi),
!                                    sin(theta)sin(phi),cos(theta)]
!
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (lmax=3)
      complex*16 eiphi,dimagx,ai
      complex*16 cp(-lmax:lmax,-lmax:lmax)
      dimension cr(-lmax:lmax,-lmax:lmax),ci(-lmax:lmax,-lmax:lmax)
!**********************************************************************
!
      DO 11 MI=-L,L
      DO 11 MF=-L,L
        CP(MI,MF)=(0.D0,0.d0)
   11 CONTINUE
      dimagx=(0.0d0,1.0d0)
      eiphi=cos(phi)+dimagx*sin(phi)
!
      if(ipr.gt.1) &
      write(6,*)' real and imag part of Ldzeta operator, L=',L
      DO 21 MF=-L,L
!
! <mf||mf>
      a=dfloat(mf)*cos(theta)
      cp(mf,mf)=a
! <mf||mf-1>
      if(mf.gt.-l)then
      ai=sin(theta)*sqrt(dfloat((l+mf)*(l-mf+1)))/(2.d0*eiphi)
      cp(mf,mf-1)=ai
      endif
! <mf||mf+1>
      if(mf.lt.l)then
      ai=sin(theta)*eiphi*sqrt(dfloat((l-mf)*(l+mf+1)))/2.d0
      cp(mf,mf+1)=ai
      endif
! <mf-1||mf>
      if(mf.gt.-l)then
      ai=sin(theta)*eiphi*sqrt(dfloat((l+mf)*(l-mf+1)))/2.d0
      cp(mf-1,mf)=ai
      endif
! <mf+1||mf>
      if(mf.lt.l)then
      ai=sin(theta)*sqrt(dfloat((l-mf)*(l+mf+1)))/(2.d0*eiphi)
      cp(mf,mf+1)=ai
      endif
!
   21 CONTINUE
!    X(l,s) = L(dzeta)
      DO 12 MI=-L,L
      DO 12 MF=-L,L
      ci(mi,mf)=aimag(cp(mi,mf))
      cr(mi,mf)=cp(mi,mf)
12    continue
      if(ipr.gt.1)then
      write(6,*)' real part'
      do 14 mi=-L,L
      write(6,100)(cr(mi,mf),mf=-L,L)
 100  format(9f8.4)
 14   continue
      write(6,*)' imag part'
      do 15 mi=-L,L
      write(6,100)(ci(mi,mf),mf=-L,L)
 15   continue
      write(6,*)
      endif
      RETURN
      END
