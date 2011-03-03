!
       DOUBLE PRECISION FUNCTION SEVAL(N,U,X,Y,B,C,D)
       INTEGER N
       REAL*8 U, X(N), Y(N), B(N), C(N), D(N)
!----------------------------------------------------------
!        FUNCTION CALCULATES VALUE OF THE CUBIC SPLINE:
!        SEVAL=Y(I)+B(I)*(U-X(I))+C(I)*(U-X(I))**2+
!              D(I)*(U-X(I)**3,
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
       SEVAL=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!
       RETURN
       END    
