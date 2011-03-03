!
      SUBROUTINE CHFAC(A,N)                                                     
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      DIMENSION A(N,N+1)                                                        
      DO 10 I=1,N-1                                                             
      DO 10 J=I+1,N                                                             
   10 A(J,I)=A(I,J)                                                             
      CALL AINV(A,N,1.D-14,IER,N)                                               
      RETURN                                                                    
      END                                                                       
