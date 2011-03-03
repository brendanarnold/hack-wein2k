      SUBROUTINE DL7UPD(BETA1, GAMMA1, L, LAMBDA, LPLUS, N, W, Z)
      implicit real*8 (a-h,o-z)
!
!  ***  COMPUTE LPLUS = SECANT UPDATE OF L  ***
!
!  ***  PARAMETER DECLARATIONS  ***
!
      INTEGER N
      DOUBLE PRECISION BETA(N), GAMMA(N), L(*), LAMBDA(N), LPLUS(*), &
                   W(N), Z(N)
!     DIMENSION L(N*(N+1)/2), LPLUS(N*(N+1)/2)
!
!--------------------------  PARAMETER USAGE  --------------------------
!
!   BETA = SCRATCH VECTOR.
!  GAMMA = SCRATCH VECTOR.
!      L (INPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE.
! LAMBDA = SCRATCH VECTOR.
!  LPLUS (OUTPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE, WHICH MAY
!             OCCUPY THE SAME STORAGE AS  L.
!      N (INPUT) LENGTH OF VECTOR PARAMETERS AND ORDER OF MATRICES.
!      W (INPUT, DESTROYED ON OUTPUT) RIGHT SINGULAR VECTOR OF RANK 1
!             CORRECTION TO  L.
!      Z (INPUT, DESTROYED ON OUTPUT) LEFT SINGULAR VECTOR OF RANK 1
!             CORRECTION TO  L.
!
!-------------------------------  NOTES  -------------------------------
!
!  ***  APPLICATION AND USAGE RESTRICTIONS  ***
!
!        THIS ROUTINE UPDATES THE CHOLESKY FACTOR  L  OF A SYMMETRIC
!     POSITIVE DEFINITE MATRIX TO WHICH A SECANT UPDATE IS BEING
!     APPLIED -- IT COMPUTES A CHOLESKY FACTOR  LPLUS  OF
!     L * (I + Z*W**T) * (I + W*Z**T) * L**T.  IT IS ASSUMED THAT  W
!     AND  Z  HAVE BEEN CHOSEN SO THAT THE UPDATED MATRIX IS STRICTLY
!     POSITIVE DEFINITE.
!
!  ***  ALGORITHM NOTES  ***
!
!        THIS CODE USES RECURRENCE 3 OF REF. 1 (WITH D(J) = 1 FOR ALL J)
!     TO COMPUTE  LPLUS  OF THE FORM  L * (I + Z*W**T) * Q,  WHERE  Q
!     IS AN ORTHOGONAL MATRIX THAT MAKES THE RESULT LOWER TRIANGULAR.
!        LPLUS MAY HAVE SOME NEGATIVE DIAGONAL ELEMENTS.
!
!  ***  REFERENCES  ***
!
! 1.  GOLDFARB, D. (1976), FACTORIZED VARIABLE METRIC METHODS FOR UNCON-
!             STRAINED OPTIMIZATION, MATH. COMPUT. 30, PP. 796-811.
!
!  ***  GENERAL  ***
!
!     CODED BY DAVID M. GAY (FALL 1979).
!     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED
!     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND
!     MCS-7906671.
!
!------------------------  EXTERNAL QUANTITIES  ------------------------
!
!  ***  INTRINSIC FUNCTIONS  ***
!/+
      DOUBLE PRECISION DSQRT
!/
!--------------------------  LOCAL VARIABLES  --------------------------
!
      INTEGER I, IJ, J, JJ, JP1, K, NM1, NP1
      DOUBLE PRECISION A, B, BJ, ETA, GJ, LJ, LIJ, LJJ, NU, S, THETA, &
                       WJ, ZJ
      DOUBLE PRECISION ONE, ZERO
!
!  ***  DATA INITIALIZATIONS  ***
!
!/6
!     DATA ONE/1.D+0/, ZERO/0.D+0/
!/7
      PARAMETER (ONE=1.D+0, ZERO=0.D+0)
!/
!
!+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
!
      NU = ONE
      ETA = ZERO
      IF (N .LE. 1) GO TO 30
      NM1 = N - 1
!
!  ***  TEMPORARILY STORE S(J) = SUM OVER K = J+1 TO N OF W(K)**2 IN
!  ***  LAMBDA(J).
!
      S = ZERO
      DO 10 I = 1, NM1
         J = N - I
         S = S + W(J+1)**2
         LAMBDA(J) = S
 10      CONTINUE
!
!  ***  COMPUTE LAMBDA, GAMMA, AND BETA BY GOLDFARB*S RECURRENCE 3.
!
      DO 20 J = 1, NM1
         WJ = W(J)
         A = NU*Z(J) - ETA*WJ
         THETA = ONE + A*WJ
         S = A*LAMBDA(J)
         LJ = DSQRT(THETA**2 + A*S)
         IF (THETA .GT. ZERO) LJ = -LJ
         LAMBDA(J) = LJ
         B = THETA*WJ + S
         GAMMA(J) = B * NU / LJ
         BETA(J) = (A - B*ETA) / LJ
         NU = -NU / LJ
         ETA = -(ETA + (A**2)/(THETA - LJ)) / LJ
 20      CONTINUE
 30   LAMBDA(N) = ONE + (NU*Z(N) - ETA*W(N))*W(N)
!
!  ***  UPDATE L, GRADUALLY OVERWRITING  W  AND  Z  WITH  L*W  AND  L*Z.
!
      NP1 = N + 1
      JJ = N * (N + 1) / 2
      DO 60 K = 1, N
         J = NP1 - K
         LJ = LAMBDA(J)
         LJJ = L(JJ)
         LPLUS(JJ) = LJ * LJJ
         WJ = W(J)
         W(J) = LJJ * WJ
         ZJ = Z(J)
         Z(J) = LJJ * ZJ
         IF (K .EQ. 1) GO TO 50
         BJ = BETA(J)
         GJ = GAMMA(J)
         IJ = JJ + J
         JP1 = J + 1
         DO 40 I = JP1, N
              LIJ = L(IJ)
          if(i.lt.0) print*, 'lplus',lplus(ij)
              LPLUS(IJ) = LJ*LIJ + BJ*W(I) + GJ*Z(I)
              W(I) = W(I) + LIJ*WJ
              Z(I) = Z(I) + LIJ*ZJ
              IJ = IJ + I
 40           CONTINUE
 50      JJ = JJ - J
 60      CONTINUE
!
 999  RETURN
!  ***  LAST CARD OF DL7UPD FOLLOWS  ***
      END
