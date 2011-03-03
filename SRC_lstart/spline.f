      SUBROUTINE SPLINE(N,X,Y,B,C,D)
      INTEGER N
      REAL*8 X(N), Y(N), B(N), C(N), D(N)
!---------------------------------------------------------
!       CALCULATE COEFFICIENTS B(I),C(I),D(I),I=1,2,...,N
!       FOR CUBIC INTERPOLATIONAL SPLINE:
!       S(X)=Y(I)+B(I)*(X-X(I))+C(I)*(X-X(I))**2+
!            D(I)*(X-X(I))**3
!       FOR X(I).LE.X.LE.X(I+1),
!       WHERE: N - NUMBER OF POINTS (N.GE.2),
!              X - ABSCISES (IN STRICTLY INCREASING ORDER),
!              Y - ORDINATES.
!---------------------------------------------------------
!
      INTEGER NM1, IB, I
      REAL*8 T
!
      NM1=N-1
      IF(N .LT. 2) RETURN
      IF(N .LT. 3) GO TO 50
!---------------------------------------------------------
!        CONSTRUCT 3-DIAGONAL SYSTEM WITH B - DIAGONAL,
!        D - OVERDIAGONAL, C - RIGHT SIDES
!---------------------------------------------------------
!
      D(1)=X(2)-X(1)
      C(2)=(Y(2)-Y(1))/D(1)
      DO 10 I=2,NM1
       D(I)=X(I+1)-X(I)
       B(I)=2.*(D(I-1)+D(I))
       C(I+1)=(Y(I+1)-Y(I))/D(I)
       C(I)=C(I+1)-C(I)
 10   CONTINUE
!----------------------------------------------------------
!        BOUNDARY CONDITIONS
!----------------------------------------------------------
!
       B(1)= -D(1)
       B(N)= -D(N-1)
       C(1)=0.D00
       C(N)=0.D00
       IF(N .EQ. 3) GO TO 15
       C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
       C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
       C(1)=C(1)*D(1)**2/(X(4)-X(1))
       C(N)= -C(N)*D(N-1)**2/(X(N)-X(N-3))
!----------------------------------------------------------
!       UPWARD STEP
!----------------------------------------------------------
!
 15    DO 20 I=2,N
        T=D(I-1)/B(I-1)
        B(I)=B(I)-T*D(I-1)
        C(I)=C(I)-T*C(I-1)
 20    CONTINUE
!----------------------------------------------------------
!       BACKWARD SUBSTITUTION
!----------------------------------------------------------
!
       C(N)=C(N)/B(N)
       DO 30 IB =1,NM1
        I=N-IB
        C(I)=(C(I)-D(I)*C(I+1))/B(I)
 30    CONTINUE
!----------------------------------------------------------
!       CALCULATE COEFFICIENTS OF THE POLINOMIALS
!----------------------------------------------------------
!
       B(N)=(Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D00*C(N))
       DO 40 I=1,NM1
        B(I)=(Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D00*C(I))
        D(I)=(C(I+1)-C(I))/D(I)
        C(I)=3.D00*C(I)
 40    CONTINUE
       C(N)=3.D00*C(N)
       D(N)=D(N-1)
       RETURN
!
 50    B(1)=(Y(2)-Y(1))/(X(2)-X(1))
       C(1)=0.D00
       D(1)=0.D00
       B(2)=B(1)
       C(2)=0.D00
       D(2)=0.D00
!
       RETURN
       END
