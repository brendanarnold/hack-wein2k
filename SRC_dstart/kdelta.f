      LOGICAL FUNCTION KDELTA(K)                                        
!                                                                       
!.... TEST, IF K IS IN STAR OF G  (GENERATED IN STERN)                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
!p    PARAMETER (NSYM=48)                                               
      COMPLEX*16       TAUP                                     
      INTEGER        G,STG                                              
      COMMON /FOUR/  G(3),NST,STG(3,NSYM),TAUP(NSYM)                           
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
