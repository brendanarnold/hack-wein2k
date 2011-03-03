      SUBROUTINE DV7CPY(P, Y, X)
      IMPLICIT REAL*8 (a-h,o-z)
!
!  ***  SET Y = X, WHERE X AND Y ARE P-VECTORS  ***
!
      INTEGER P
      DOUBLE PRECISION X(P), Y(P)
!
      INTEGER I
!
      DO 10 I = 1, P
 10      Y(I) = X(I)
      RETURN
      END
