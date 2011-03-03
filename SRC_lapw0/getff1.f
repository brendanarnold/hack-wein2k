      SUBROUTINE GETFF1(nkk,ifft1,ifft2,ifft3,rhoksp,fft,kzz,inst)
!
!.....sets rhoksp from fft-field
!
      use struct
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16 FFT(IFFT1,IFFT2,IFFT3,2),rhoksp(nkk,2)          
      dimension kzz(3,NKK),inst(NKK)
      DO 22 J=1,NKK                                                     
!     transfer karth into internal coord.
      if(ortho) then
         ii1=nint((avec(1,1)*kzz(1,j)+avec(1,2)*kzz(2,j)+ &
         avec(1,3)*kzz(3,j))/aa)
         ii2=nint((avec(2,1)*kzz(1,j)+avec(2,2)*kzz(2,j)+ &
         avec(2,3)*kzz(3,j))/bb)
         ii3=nint((avec(3,1)*kzz(1,j)+avec(3,2)*kzz(2,j)+ &
         avec(3,3)*kzz(3,j))/cc)
      else
         ii1=kzz(1,j)
         ii2=kzz(2,j)
         ii3=kzz(3,j)
      end if
      IF(IABS(ii1).GE.IFFT1/2) GOTO 22                             
      IF(IABS(ii2).GE.IFFT2/2) GOTO 22                             
      IF(IABS(ii3).GE.IFFT3/2) GOTO 22                             
      I1=ii1                                                       
        IF(I1.LT.0) I1=I1+IFFT1                                         
      I2=ii2                                                       
        IF(I2.LT.0) I2=I2+IFFT2                                         
      I3=ii3                                                       
        IF(I3.LT.0) I3=I3+IFFT3                                         
      rhoksp(J,1)=FFT(I1+1,I2+1,I3+1,1)*INST(J)                  
      rhoksp(J,2)=FFT(I1+1,I2+1,I3+1,2)*INST(J)                       
 22   continue
      return
      end
