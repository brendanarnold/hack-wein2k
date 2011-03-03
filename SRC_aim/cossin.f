      SUBROUTINE cossin(SINTH,COSTH,SINPH,COSPH,AB,V)
      implicit real*8 (a-h,o-z)
      DOUBLE PRECISION   V(3)
!
      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      DOUBLE PRECISION   A, B, C, AB, ABC, ABMAX, ABCMAX
      DOUBLE PRECISION   D4LL1C, D2L13, PI
      DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
!
      PI = (4.0D+0)*ATAN(1.0D+0)
!        calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
!        Theta, Phi ... polar angles of vector V
!
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. 0.0D+0) THEN
         A = V(1)/ABMAX
         B = V(2)/ABMAX
         AB = 1.D0/SQRT(A*A+B*B)
         COSPH = A*AB
         SINPH = B*AB
      ELSE
         COSPH = 1.0D+0
         SINPH = 0.0D+0
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. 0.0D+0) THEN
         A = V(1)/ABCMAX
         B = V(2)/ABCMAX
         C = V(3)/ABCMAX
         AB = A*A + B*B
         ABC = 1.D0/SQRT(AB + C*C)
         COSTH = C*ABC
         SINTH = SQRT(AB)*ABC
      ELSE
         COSTH = 1.0D+0
         SINTH = 0.0D+0
      ENDIF
      ab=sqrt(V(1)*V(1)+V(2)*V(2))
  999 RETURN
!
!        End of 'cossin
!
      END
