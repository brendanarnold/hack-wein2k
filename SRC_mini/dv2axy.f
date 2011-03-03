      SUBROUTINE DV2AXY(P, W, A, X, Y)
      IMPLICIT REAL*8 (a-h,o-z)
!
!  ***  SET W = A*X + Y  --  W, X, Y = P-VECTORS, A = SCALAR  ***
!
      INTEGER P
      DOUBLE PRECISION A, W(P), X(P), Y(P)
!
      INTEGER I
!
      DO 10 I = 1, P
 10      W(I) = A*X(I) + Y(I)
      RETURN
      END
