      SUBROUTINE REAN3 (NKK,KZZ,CVALUE,IFFT1,IFFT2,IFFT3,FFT,U, &
       u0)       
!                                                                       
!     REANALYSE THE k=0 FOURIER COEFFICENT USING ARRAY U (GENERATED        
!     IN REAN0), determine v0                                                 
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16         U,FFT                                          
      COMPLEX*16         CVALUE(1),vmean                                    
!                                                                       
      DIMENSION U(-IFFT1/2:IFFT1/2,-IFFT2/2:IFFT2/2,-IFFT3/2:IFFT3/2)         
      DIMENSION FFT(IFFT1,IFFT2,IFFT3)                            
      DIMENSION KZZ(3,NKK)                                        
!                                                                       
!--------------------------------------------------------------------   
!                   
         vmean=fft(1,1,1)
         fft(1,1,1)=0.d0
         cvalue(1)=0.d0                                                     
         jre=1  
         DO 5000 I1=1,IFFT1                                             
            IKX=I1-1                                                    
            IF(I1.GT.IFFT1/2) IKX=IKX-IFFT1                             
         DO 5000 I2=1,IFFT2                                             
            IKY=I2-1                                                    
            IF(I2.GT.IFFT2/2) IKY=IKY-IFFT2                             
         DO 5000 I3=1,IFFT3                                             
            IKZ=I3-1                                                    
            IF(I3.GT.IFFT3/2) IKZ=IKZ-IFFT3                             
      CVALUE(JRE)=CVALUE(JRE) + U(ikx,iky,ikz)*FFT(I1,I2,I3)               
!      if(abs(ikx).lt.4.and.abs(iky).lt.4.and.abs(ikz).lt.4)
!     *write(6,2939) ikx,iky,ikz,u(ikx,iky,ikz),fft(i1,i2,i3)
!     *,U(ikx,iky,ikz)*FFT(I1,I2,I3)
! 2939 format(3i3,3(2f12.8))
                                               
 4000 CONTINUE                                                          
 5000    CONTINUE                     
      cvalue(1)=cvalue(1)/u(0,0,0)
      fft(1,1,1)=-cvalue(1)
      u0= u(0,0,0)                                                
!      write(6,*) ' v-mean, v0,u0',cvalue(1), vmean+cvalue(1),u(0,0,0)
      cvalue(1)=vmean+cvalue(1)
      RETURN                                                            
      END                                                               
