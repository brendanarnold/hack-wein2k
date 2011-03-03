      FUNCTION DFAC(L)                                                  
!CC CALCULATES DOUBLE FACTORIALS OF ARGUMENT 2*L+1                      
      IMPLICIT REAL*8 (A-H,O-Z)
      DFAC = 1.d0                                                         
      IF (L.EQ.0) RETURN                                                
      X = 1.d0                                                            
      DO 10 I=1, L                                                      
      X = X+2.d0                                                          
   10 DFAC = DFAC*X                                                     
      RETURN                                                            
      END                                                               
