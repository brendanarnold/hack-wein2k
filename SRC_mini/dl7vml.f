      SUBROUTINE DL7VML(N, X, L, Y)
      implicit real*8 (a-h,o-z)
!
!  ***  COMPUTE  X = L*Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR
!  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
!  ***  STORAGE.  ***
!
      INTEGER N
      DOUBLE PRECISION X(N), L(1), Y(N)
!     DIMENSION L(N*(N+1)/2)
      INTEGER I, II, IJ, I0, J, NP1
      DOUBLE PRECISION T, ZERO
!/6
!     DATA ZERO/0.D+0/
!/7
      PARAMETER (ZERO=0.D+0)
!/
!
      NP1 = N + 1
      I0 = N*(N+1)/2
      DO 20 II = 1, N
         I = NP1 - II
         I0 = I0 - I
         T = ZERO
         DO 10 J = 1, I
              IJ = I0 + J
              T = T + L(IJ)*Y(J)
 10           CONTINUE
         X(I) = T
 20      CONTINUE
 999  RETURN
!  ***  LAST CARD OF DL7VML FOLLOWS  ***
      END
