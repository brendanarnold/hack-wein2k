      SUBROUTINE DV7SCP(P, Y, S)
      IMPLICIT REAL*8 (a-h,o-z)
!
!  ***  SET P-VECTOR Y TO SCALAR S  ***
!
      INTEGER P
      DOUBLE PRECISION S, Y(P)
!
      INTEGER I
!
      DO 10 I = 1, P
 10      Y(I) = S
      RETURN
      END
