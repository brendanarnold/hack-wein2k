      SUBROUTINE CHSLV(A,N,B)                                                   
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      DIMENSION A(N,N+1),B(N)                                                   
      DO 1 I=1,N                                                                 
    1 A(I,N+1)=B(I)                                                             
      DO 2 I=1,N                                                                 
      B(I)=0.D0                                                                 
      DO 2 J=1,N                                                                 
    2 B(I)=B(I)+A(J,I)*A(J,N+1)                                                 
      RETURN                                                                    
      END                                                                       
