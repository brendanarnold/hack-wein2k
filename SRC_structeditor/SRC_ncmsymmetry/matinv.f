

      SUBROUTINE MATINV (A,N,NP,B,AUX)

!* INVERTS A SQUARE MATRIX USING ROUTINES OF "NUMERICAL RECIPES",
!*   W.PRESS ET AL, CAMBRIDGE U.P. 1987.
!* ARGUMENTS:
!*   A   : INPUT SQUARE MATRIX
!*   N   : INPUT DIMENSION OF THE MATRIX
!*   NP  : INPUT FIRST PHYSICAL DIMENSION OF MATRICES A AND B
!*   B   : OUTPUT INVERSE MATRIX. MAKE B=A FOR 'IN PLACE' INVERSION
!*   AUX : AUXILIARY ARRAY OF MINIMUN SIZE N*N. IF A AND B ARE DIFFERENT
!*         AND MATRIX A NEEDS NOT BE PRESERVED, YOU CAN MAKE AUX=A
!* WRITTEN BY JOSE SOLER (JSOLER AT EMDUAM11) 22/4/90.

!*     FOR SINGLE PRECISION SIMPLY COMMENT OUT NEXT LINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=100,ZERO=0.D0,ONE=1.D0)
      DIMENSION A(NP,NP),B(NP,NP),AUX(N,N),INDEX(NMAX)
      IF (N.GT.NMAX) THEN
        WRITE(6,*) ' MATINV: DIMENSION TOO SMALL. N,NMAX=',N,NMAX
        STOP
      ENDIF
      DO 10 I=1,N
      DO 10 J=1,N
        AUX(J,I)=A(J,I)
   10 CONTINUE
      CALL LUDCMP (AUX,N,N,INDEX,D)
      DO 30 I=1,N
        DO 20 J=1,N
          B(J,I)=ZERO
   20   CONTINUE
        B(I,I)=ONE
        CALL LUBKSB (AUX,N,N,INDEX,B(1,I))
   30 CONTINUE
      END