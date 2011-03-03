      SUBROUTINE REAN3 (NKK,KZZ,CVALUE,IFFT1,IFFT2,IFFT3,FFT,U, &
       u0)       
!                                                                       
!     REANALYSE THE k=0 FOURIER COEFFICENT USING ARRAY U (GENERATED        
!     IN REAN1), determine v0                                                 
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16         U,FFT                                          
      COMPLEX*16         CVALUE(NKK),vmean                                    
!                                                                       
      DIMENSION       U(ifft1,ifft2,ifft3)
      DIMENSION       FFT(IFFT1,IFFT2,IFFT3)                            
      DIMENSION       KZZ(3,NKK)                                        
!                                                                       
!--------------------------------------------------------------------   
!                   
      vmean=fft(1,1,1)
      fft(1,1,1)=0.d0
      cvalue(1)=0.d0                                                     
      DO 5000 I3=1,IFFT3                                             
         IKZ=I3-1                                                    
         IF(I3.GT.IFFT3/2) IKZ=IKZ-IFFT3                             
         DO 5000 I2=1,IFFT2                                             
            IKY=I2-1                                                    
            IF(I2.GT.IFFT2/2) IKY=IKY-IFFT2                             
            DO 5000 I1=1,IFFT1                                             
               IKX=I1-1                                                    
               IF(I1.GT.IFFT1/2) IKX=IKX-IFFT1                             
!lh               K1=IKX-KZZ(1,1) ! kzz(*,1) has to be zero in order to obtain
!lh               K2=IKY-KZZ(2,1) ! the integral of the effective
!lh               K3=IKZ-KZZ(3,1) ! potential without MTs
!lh               CVALUE(1)=CVALUE(1) + U(K1,K2,K3)*FFT(I1,I2,I3)               
               CVALUE(1)=CVALUE(1) + U(i1,i2,i3)*FFT(I1,I2,I3)               
                                               
 5000 CONTINUE                     
      cvalue(1)=cvalue(1)/u(1,1,1)
      fft(1,1,1)=-cvalue(1)
      u0= u(1,1,1)                                                
      write(6,*) ' v-mean, v0,u0',cvalue(1), vmean+cvalue(1),u(1,1,1)
      cvalue(1)=vmean+cvalue(1)
      RETURN                                                            
      END                                                               
