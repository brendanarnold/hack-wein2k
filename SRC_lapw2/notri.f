      FUNCTION NOTRI(K,L,M)                                             
      IMPLICIT REAL*8 (A-H,O-Z)
      NOTRI=-1                                                          
      IF(MOD((K+L+M),2).EQ.1) GOTO 10                                   
      IF((K+L-M).LT.0) GOTO 10                                          
      IF((K-L+M).LT.0) GOTO 10                                          
      IF((M+L-K).LT.0) GOTO 10                                          
      NOTRI=1                                                           
   10 RETURN                                                            
      END                                                               
