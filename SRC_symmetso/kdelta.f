      LOGICAL FUNCTION KDELTA(K)                                        
!                                                                       
!.... TEST, IF K IS IN STAR OF G  (GENERATED IN STERN)                  
!                                                                       
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
      INTEGER        STG                                              
      COMMON /FOUR/  NST,STG(3,noper)                               
      DIMENSION      K(3)                                               
!---------------------------------------------------------------------  
      DO 1 I=1,NST                                                      
      DO 2 J=1,3 
      IF(STG(J,I).NE.K(J)) GOTO 1                                       
   2  CONTINUE                                                          
      KDELTA=.TRUE.                                                     
      RETURN                                                            
  1   CONTINUE                                                          
      KDELTA=.FALSE.                                                    
      RETURN                                                            
      END                                                               
