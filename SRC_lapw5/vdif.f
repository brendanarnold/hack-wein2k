      SUBROUTINE VDIF(V1,V2,V)                                          
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V1(3),V2(3),V(3)                                        
      DO 1 I=1,3                                                        
   1  V(I)=V2(I)-V1(I)                                                  
      RETURN                                                            
      END                                                               
