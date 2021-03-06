      DOUBLE PRECISION FUNCTION DD7TPR(P, X, Y)
      implicit real*8 (a-h,o-z)
!
!  ***  RETURN THE INNER PRODUCT OF THE P-VECTORS X AND Y.  ***
!
      INTEGER P
      DOUBLE PRECISION X(P), Y(P)
!
      INTEGER I
      DOUBLE PRECISION ONE, SQTETA, T, ZERO
      EXTERNAL DR7MDC
      DOUBLE PRECISION DR7MDC
!
!  ***  DR7MDC(2) RETURNS A MACHINE-DEPENDENT CONSTANT, SQTETA, WHICH
!  ***  IS SLIGHTLY LARGER THAN THE SMALLEST POSITIVE NUMBER THAT
!  ***  CAN BE SQUARED WITHOUT UNDERFLOWING.
!
!/6
!     DATA ONE/1.D+0/, SQTETA/0.D+0/, ZERO/0.D+0/
!/7
      PARAMETER (ONE=1.D+0, ZERO=0.D+0)
      DATA SQTETA/0.D+0/
!/
!
      DD7TPR = ZERO
      IF (P .LE. 0) GO TO 999
      IF (SQTETA .EQ. ZERO) SQTETA = DR7MDC(2)
      DO 20 I = 1, P
         T = DMAX1(DABS(X(I)), DABS(Y(I)))
         IF (T .GT. ONE) GO TO 10
         IF (T .LT. SQTETA) GO TO 20
         T = (X(I)/SQTETA)*Y(I)
         IF (DABS(T) .LT. SQTETA) GO TO 20
 10      DD7TPR = DD7TPR + X(I)*Y(I)
 20   CONTINUE
!
 999  RETURN
!  ***  LAST CARD OF DD7TPR FOLLOWS  ***
      END
