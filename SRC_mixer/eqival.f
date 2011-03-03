      LOGICAL FUNCTION EQIVAL (KVEC1,KVEC2)                             
!                                                                       
!----------------------------------------------------------------       
!---  EQIVAL IS TRUE IF KVEC1=KVEC2, ELSE EQIVAL IS    FALSE --         
!----------------------------------------------------------------       
!                                                                       
      DIMENSION KVEC1(3),KVEC2(3)                                       
!                                                                       
      EQIVAL=.TRUE.                                                     
      DO 10 J=1,3                                                       
         IF( KVEC1(J).NE.KVEC2(J) ) EQIVAL=.FALSE.                      
 10   CONTINUE                                                          
      RETURN                                                            
      END                                                               
