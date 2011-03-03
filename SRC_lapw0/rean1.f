      SUBROUTINE REAN1 (NKK,KZZ,CVALUE,IFFT1,IFFT2,IFFT3,FFT,ispin,U)         
!                                                                       
!     REANALYSE THE FOURIER COEFFICENTS USING ARRAY U (GENERATED        
!     IN REAN1)                                                         
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16         U,FFT,ctem                                          
      COMPLEX*16         CVALUE(NKK)                                    
!                                                                       
      DIMENSION       U(-IFFT1:IFFT1,-IFFT2:IFFT2,-IFFT3:IFFT3)         
      DIMENSION       FFT(IFFT1,IFFT2,IFFT3,2)                            
      DIMENSION       KZZ(3,NKK)                                        
!                                                                       
!--------------------------------------------------------------------   
!                                                                       
         DO 5000 I1=1,IFFT1                                             
            IKX=I1-1                                                    
            IF(I1.GT.IFFT1/2) IKX=IKX-IFFT1                             
         DO 5000 I2=1,IFFT2                                             
            IKY=I2-1                                                    
            IF(I2.GT.IFFT2/2) IKY=IKY-IFFT2                             
         DO 5000 I3=1,IFFT3                                             
            IKZ=I3-1                                                    
            IF(I3.GT.IFFT3/2) IKZ=IKZ-IFFT3  
            ctem=fft(i1,i2,i3,ispin)                           
               DO 4000 JRE=1,NKK                                        
               K1=IKX-KZZ(1,JRE)                                        
               K2=IKY-KZZ(2,JRE)                                        
               K3=IKZ-KZZ(3,JRE)
      CVALUE(JRE)=CVALUE(JRE) + U(K1,K2,K3)*ctem               
                                                                        
 4000 CONTINUE                                                          
 5000    CONTINUE                                                       
      RETURN                                                            
      END                                                               
