      SUBROUTINE INTERP(R,P,N,RS,PS,DPS,DERIV)                          
!                                                                       
!     ******************************************************************
!                                                                       
!     INTERPOLATE VIA LAGRANGE                                          
!                                                                       
!     ******************************************************************
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL DERIV,NODRIV                                              
      DIMENSION R(N),P(N)                                               
      NODRIV=.NOT.DERIV                                                 
      PS=0.E0                                                           
      DPS=0.E0                                                          
      DO 1 J=1,N                                                        
      TERM=1.E0                                                         
      DENOM=1.E0                                                        
      DTERM=0.E0                                                        
      DO 2 I=1,N                                                        
      IF(I.EQ.J) GO TO 2                                                
      DENOM=DENOM*(R(J)-R(I))                                           
      TERM=TERM*(RS-R(I))                                               
      IF(NODRIV) GO TO 2                                                
      DTERM1=1.E0                                                       
      DO 3 K=1,N                                                        
      IF(K.EQ.J.OR.K.EQ.I) GO TO 3                                      
      DTERM1=DTERM1*(RS-R(K))                                           
    3 CONTINUE                                                          
      DTERM=DTERM+DTERM1                                                
    2 CONTINUE                                                          
      IF(NODRIV) GO TO 1                                                
      DPS=DPS+DTERM*P(J)/DENOM                                          
    1 PS=PS+TERM *P(J)/DENOM                                            
      RETURN                                                            
      END                                                               
