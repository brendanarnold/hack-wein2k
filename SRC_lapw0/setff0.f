      SUBROUTINE SETFF0(ifft1,ifft2,ifft3,UST,fft)
!
!.....sets values from ust into fft-field,just one spin!!
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      COMPLEX*16 FFT(IFFT1,IFFT2,IFFT3)         
      COMPLEX*16 UST(-IFFT1:IFFT1,-IFFT2:IFFT2,-IFFT3:IFFT3)
!
      DO 9 I1=1,IFFT1
            IKX=I1-1                                                    
            IF(I1.GT.IFFT1/2) IKX=IKX-IFFT1                             
      DO 9 I2=1,IFFT2
            IKY=I2-1                                                    
            IF(I2.GT.IFFT2/2) IKY=IKY-IFFT2                             
      DO 9 I3=1,IFFT3                                                   
            IKZ=I3-1                                                    
            IF(I3.GT.IFFT3/2) IKZ=IKZ-IFFT3  
  9   FFT(I1,I2,I3)=conjg(UST(IKX,IKY,IKZ))                                           
      return
      end
