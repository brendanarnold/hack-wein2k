      DOUBLE PRECISION FUNCTION DV2NRM(P, X)
      IMPLICIT REAL*8 (a-h,o-z)
!
!  ***  RETURN THE 2-NORM OF THE P-VECTOR X, TAKING  ***
!  ***  CARE TO AVOID THE MOST LIKELY UNDERFLOWS.    ***
!
      INTEGER P
      DOUBLE PRECISION X(P)
!
      INTEGER I, J
      DOUBLE PRECISION ONE, R, SCALE, SQTETA, T, XI, ZERO
!/+
      DOUBLE PRECISION DABS, DSQRT
!/
      EXTERNAL DR7MDC
      DOUBLE PRECISION DR7MDC
!
!/6
!     DATA ONE/1.D+0/, ZERO/0.D+0/
!/7
      PARAMETER (ONE=1.D+0, ZERO=0.D+0)
      SAVE SQTETA
!/
      DATA SQTETA/0.D+0/
!
      IF (P .GT. 0) GO TO 10
         DV2NRM = ZERO
         GO TO 999
 10   DO 20 I = 1, P
         IF (X(I) .NE. ZERO) GO TO 30
 20      CONTINUE
      DV2NRM = ZERO
      GO TO 999
!
 30   SCALE = DABS(X(I))
      IF (I .LT. P) GO TO 40
         DV2NRM = SCALE
         GO TO 999
 40   T = ONE
      IF (SQTETA .EQ. ZERO) SQTETA = DR7MDC(2)
!
!     ***  SQTETA IS (SLIGHTLY LARGER THAN) THE SQUARE ROOT OF THE
!     ***  SMALLEST POSITIVE FLOATING POINT NUMBER ON THE MACHINE.
!     ***  THE TESTS INVOLVING SQTETA ARE DONE TO PREVENT UNDERFLOWS.
!
      J = I + 1
      DO 60 I = J, P
         XI = DABS(X(I))
         IF (XI .GT. SCALE) GO TO 50
              R = XI / SCALE
              IF (R .GT. SQTETA) T = T + R*R
              GO TO 60
 50           R = SCALE / XI
              IF (R .LE. SQTETA) R = ZERO
              T = ONE  +  T * R*R
              SCALE = XI
 60      CONTINUE
!
      DV2NRM = SCALE * DSQRT(T)
 999  RETURN
!  ***  LAST CARD OF DV2NRM FOLLOWS  ***
      END
