      DOUBLE PRECISION FUNCTION DRLDST(P, D, X, X0)
      implicit real*8 (a-h,o-z)
!
!  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  ***
!  ***  NL2SOL VERSION 2.2  ***
!
      INTEGER P
      DOUBLE PRECISION D(P), X(P), X0(P)
!
      INTEGER I
      DOUBLE PRECISION EMAX, T, XMAX, ZERO
!/6
!     DATA ZERO/0.D+0/
!/7
      PARAMETER (ZERO=0.D+0)
!/
!
!  ***  BODY  ***
!
      EMAX = ZERO
      XMAX = ZERO
      DO 10 I = 1, P
         T = DABS(D(I) * (X(I) - X0(I)))
         IF (EMAX .LT. T) EMAX = T
         T = D(I) * (DABS(X(I)) + DABS(X0(I)))
         IF (XMAX .LT. T) XMAX = T
 10      CONTINUE
      DRLDST = ZERO
      IF (XMAX .GT. ZERO) DRLDST = EMAX / XMAX
 999  RETURN
!  ***  LAST CARD OF DRLDST FOLLOWS  ***
      END
