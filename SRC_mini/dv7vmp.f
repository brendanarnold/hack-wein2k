      SUBROUTINE DV7VMP(N, X, Y, Z, K)
      IMPLICIT REAL*8 (a-h,o-z)
!
! ***  SET X(I) = Y(I) * Z(I)**K, 1 .LE. I .LE. N (FOR K = 1 OR -1)  ***
!
      INTEGER N, K
      DOUBLE PRECISION X(N), Y(N), Z(N)
      INTEGER I
!
      IF (K .GE. 0) GO TO 20
      DO 10 I = 1, N
 10      X(I) = Y(I) / Z(I)
      GO TO 999
!
 20   DO 30 I = 1, N
 30      X(I) = Y(I) * Z(I)
 999  RETURN
!  ***  LAST CARD OF DV7VMP FOLLOWS  ***
      END
