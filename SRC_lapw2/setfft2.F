SUBROUTINE setfft2(nkk,ifft1,ifft2,ifft3,a,tp,ii,s,t,z,keigen)
  !
  ! Sets the kinetic energy fft-fields from the eigenvectors a and k-vectors
  ! GM 6/6-00
  USE defs
  USE struk
  USE param
  IMPLICIT NONE

!_REAL      REAL*8           A(nmat)
!_COMPLEX      COMPLEX*16           A(nmat)
  INTEGER, intent(in) :: nkk,ifft1,ifft2,ifft3,ii
  INTEGER     KEIGEN(3,NMAT)
  REAL*8, intent(in) :: s,t,z
  COMPLEX*16, intent(out) :: tp(IFFT1,IFFT2,IFFT3)
  
  INTEGER             :: I, I1, I2, I3
  REAL*8              :: TPI
  
  tpi=two*pi
  tp=zeroc
  
  DO I=1,NKK
     i1=keigen(1,I)
     i2=keigen(2,I)
     i3=keigen(3,I)
     IF(I1.LT.0) I1=I1+IFFT1                                         
     IF(I2.LT.0) I2=I2+IFFT2                                         
     IF(I3.LT.0) I3=I3+IFFT3
     TP(I1+1,I2+1,I3+1)=((s+dble(keigen(1,i)))*BR1(ii,1) &
          +(t+dble(keigen(2,i)))*BR1(ii,2) &
          +(z+dble(keigen(3,i)))*BR1(ii,3))*a(i)
  ENDDO
END SUBROUTINE setfft2

