      FUNCTION AINTER(R,P,RS)                                           
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(4),P(4)                                               
      N=4                                                               
      PS=0.d0                                                             
      DO 1 J=1,N                                                        
      TERM=1.d0                                                           
      DENOM=1.d0                                                          
      DO 2 I=1,N                                                        
      IF(I.EQ.J) GOTO 2                                                 
      DENOM=DENOM*(R(J)-R(I))                                           
      TERM=TERM*(RS-R(I))                                               
  2   CONTINUE                                                          
  1   PS=PS+TERM*P(J)/DENOM                                             
       AINTER=PS                                                        
      RETURN                                                            
      END                                                               
