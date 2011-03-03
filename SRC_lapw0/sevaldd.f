!
       DOUBLE PRECISION FUNCTION SEVALDD(N,U,X,B,C,D)
       INTEGER N
       REAL*8 U, X(N), B(N), C(N), D(N)
!----------------------------------------------------------
!        CALCULATE SECOND DERIVATIVE OF SPLINE FOR POINT U
!        WHERE X(I) .LT. U .LT. X(I+1)
!----------------------------------------------------------
!
       INTEGER I,J,K
       REAL*8 DX
       DATA I/1/
       IF(I .GE. N) I=1
       IF(U .LT. X(I)) GO TO 10
       IF(U .LE. X(I+1)) GO TO 30
!----------------------------------------------------------
!        SEARCHING OF THE INTERVAL
!----------------------------------------------------------
!
 10    I=1
       J=N+1
 20    K=(I+J)/2
       IF(U .LT. X(K)) J=K
       IF(U .GE. X(K)) I=K
       IF(J .GT. I+1) GO TO 20
!----------------------------------------------------------
!        THE FINAL PART
!----------------------------------------------------------
!
 30    DX=U-X(I)
       SEVALDD=2.D0*C(I)+6.D0*DX*D(I)
!
       RETURN
       END    
