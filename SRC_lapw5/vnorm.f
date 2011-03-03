      REAL*8 FUNCTION VNORM(V)                                                 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(3)                                                    
      VNORM=0.0                                                         
      DO 1 I=1,3                                                        
    1 VNORM=VNORM+V(I)**2                                               
      VNORM=SQRT(VNORM)                                                 
      RETURN                                                            
      END                                                               
