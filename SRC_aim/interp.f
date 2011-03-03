      SUBROUTINE INTERP(R,P,N,RS,PS,DPS,DERIV)                          
!     
!     ******************************************************************
!     
!     INTERPOLATE VIA LAGRANGE                                          
!     
!     ******************************************************************
!     
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL DERIV
      DIMENSION R(N),P(N)                                               
      PS=0.D0                                                           
      DPS=0.D0                                                          
      DO  J=1,N                                                        
        TERM=1.E0                                                         
        DENOM=1.E0                                                        
        DTERM=0.E0                                                        
        DO  I=1,N                                                        
          IF(I.ne.J) then                                                
            DENOM=DENOM*(R(J)-R(I))                                           
            TERM=TERM*(RS-R(I))                                               
            IF(deriv) then
              DTERM1=1.D0
              DO K=1,N                                                        
                IF(.not.(K.EQ.J.OR.K.EQ.I)) then
                  DTERM1=DTERM1*(RS-R(K))         
                endif
              enddo                                                          
              DTERM=DTERM+DTERM1
            endif
          endif
        enddo
        IF(deriv) then                                                
          DPS=DPS+DTERM*P(J)/DENOM                                          
        endif
      PS=PS+TERM *P(J)/DENOM
      enddo
      RETURN                                                            
      END
                                                
