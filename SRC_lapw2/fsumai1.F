SUBROUTINE fsumai1(ekink,nwave)

  USE param; USE defs
  USE struk
  USE reclat; USE com
  IMPLICIT NONE

!_REAL      REAL*8           factor
!_COMPLEX      COMPLEX*16           factor
  INTEGER                :: ia,nwave,jneq,j,jj,ik
  INTEGER                :: G(3),NST,STG(3,NSYM)
  real*8                 :: fj(0:lmax2),nkzz(3),rotkzz(3)
  REAL*8                 :: lenkzz,fstar(3),frmt,pha,gl
  complex*16             :: ekink(*),TAUP(NSYM)
  real*8,allocatable     :: fsur(:,:)

  allocate (fsur(0:3,nat))
  fsur=zero
  ia=1
  DO jneq=1,nat
     frmt=4.d0*pi*rmt(jneq)**2
     ! skip G=(0,0,0): no contribution to ekin
     DO j=2,nwave
        nkzz=MATMUL(br1,kzz(1:3,j))
        lenkzz=SQRT(SUM(nkzz(1:3)**2))
        CALL SPHBES(lmax2,RMT(jneq)*lenkzz,fj)
        g(1:3)=kzz(1:3,j)
        CALL STERN(g,nst,STG,taup)
        fstar=zero
        DO jj=1,nst
           nkzz(1:3)=stg(1:3,jj)
           CALL ROTATE(nkzz,rotij(1,1,ia),rotkzz)
           !     Calculate phase factor
           pha=2.0d0*pi*SUM(pos(1:3,ia)*nkzz(1:3))
           nkzz(1:3)=MATMUL(br1(1:3,1:3),rotkzz(1:3))
           CALL ROTATE (nkzz,ROTLOC(1,1,jneq),rotkzz)  
           !     Calculate |G|
           gl=SQRT(SUM(rotkzz(1:3)**2))
           !     Update surface force
           !cc   swaped sin and cos because there's an i in the equation
!_REAL          factor=ekink(j)*(-sin(pha))*taup(jj)*fj(1)/gl
!_COMPLEX          factor=ekink(j)*cmplx(-sin(pha),cos(pha))*taup(jj)*fj(1)/gl
           fstar(1:3)=fstar(1:3)+factor*rotkzz(1:3)
        ENDDO
        fsur(1:3,jneq)=fsur(1:3,jneq)+fstar(1:3)/nst*frmt
     ENDDO
     fsur(0,jneq)=SQRT(SUM(fsur(1:3,jneq)**2))
     WRITE(6,77) jneq,1,'+ SUR',(fsur(ik,jneq),ik=0,3)
     ia=ia+mult(jneq)
     WRITE(21,79) jneq,jneq,(fsur(ik,jneq)*1000,ik=0,3)
  ENDDO
  DEALLOCATE(fsur)
77 FORMAT(2i3,a7,4e15.7)
79 FORMAT(':FSU',i3.3,':',1x,i3,'.ATOM',4f15.3)
END SUBROUTINE fsumai1

