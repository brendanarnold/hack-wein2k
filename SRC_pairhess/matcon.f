      SUBROUTINE MATCNO (A,LDA,N,CONDNO,IFLAG)
!
!  --------------------------------------------------------------------
!   MATCNO   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
!            DIVISION, NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY,
!            GAITHERSBURG, MARYLAND  20899
!
!   FOR: COMPUTING THE CONDITION NUMBER OF THE N BY N MATRIX A.  THE
!        CONDITION NUMBER IS KAPPA = NORM(A)*NORM[INV(A)].  HERE THE
!        INFINITY-NORM IS USED.  THE CONDITION NUMBER MEASURES HOW
!        "ILL-CONDITIONED" THE MATRIX A IS, AND THUS HOW ACCURATELY THE
!        SYSTEM A*X = B MAY BE SOLVED.  IF A IS SINGULAR THEN KAPPA IS
!        INFINITE AND AN ERROR FLAG IS RETURNED.  SEE THE REFERENCE FOR
!        A GOOD DISCUSSION OF THE PROBLEM.
!
!   SUBPROGRAMS CALLED: MATINV (STSPAC)
!                       ANORMI (ATTACHED)
!
!   CURRENT VERSION COMPLETED JUNE 28, 1989
!
!   REFERENCE: GOLUB, G.H. AND VAN LOAN, C.F., "MATRIX COMPUTATIONS",
!              JOHNS HOPKINS UNIVERSITY PRESS, BALTIMORE, 1985, P. 72.
!  --------------------------------------------------------------------
!   DEFINITION OF PASSED PARAMETERS:
!
!   * A(LDA,*) = MATRIX (SIZE N BY N) WHOSE CONDITION NUMBER IS TO BE
!                DETERMINED.  THE SECOND DIMENSION OF A, DEFINED IN THE
!                CALLING PROGRAM, MUST BE >=N [REAL]
!
!        * LDA = THE LEADING DIMENSION OF A (>=N) [INTEGER]
!
!          * N = SIZE OF MATRIX A (N<=LDA) [INTEGER]
!
!       CONDNO = THE CONDITION NUMBER OF A (KAPPA) ON RETURN WHEN
!                IFLAG=0 [REAL]
!
!        IFLAG = ERROR INDICATOR ON OUTPUT.  THE OUTPUT VALUE OF IFLAG
!                FROM SUBROUTINE MATINV IS PASSED DIRECTLY THROUGH THIS
!                SUBROUTINE.  IF IFLAG=0, THEN ALL IS WELL.  OTHERWISE,
!                SEE DOCUMENTATION OF MATINV FOR INTERPRETATION OF THE
!                THE ERROR FLAG IFLAG [INTEGER]
!
!   * INDICATES PARAMETERS REQUIRING INPUT VALUES
!  --------------------------------------------------------------------
!
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION A(LDA,*)
      FAC1 = ANORMI(A,LDA,N)
      CALL MATINV (A,LDA,N,IFLAG)
      IF (IFLAG.NE.0) THEN
         CONDNO = 0.999999E29
      ELSE
         FAC2 = ANORMI(A,LDA,N)
         CONDNO = FAC1*FAC2
      ENDIF
      RETURN
!
      END
!
!=======================================================================
!
      DOUBLE Precision FUNCTION ANORMI (A,LDA,N)
!
!   COMPUTE THE INFINITY-NORM OF THE N BY N MATRIX A
!
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION A(LDA,*)
      ANORMI = 0.0
      DO 20 I = 1, N
         T = 0.0
         DO 10 J = 1, N
            T = T+ABS(A(I,J))
   10    CONTINUE
         ANORMI = MAX(ANORMI,T)
   20 CONTINUE
      RETURN
!
      END
!
      SUBROUTINE MATINV (A,LDA,N,IFLAG)
!
!  --------------------------------------------------------------------
!   MATINV   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
!            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
!            MARYLAND  20899
!
!   FOR: COMPUTING THE INVERSE OF A GENERAL N BY N MATRIX IN PLACE,
!        I.E., THE INVERSE OVERWRITES THE ORIGINAL MATRIX.  THE STEPS
!        OF THE ALGORITHM ARE DESCRIBED BELOW AS THEY OCCUR.  ROW
!        INTERCHANGES ARE DONE AS NEEDED IN ORDER TO INCREASE THE
!        ACCURACY OF THE INVERSE MATRIX.  WITHOUT INTERCHANGES THIS
!        ALGORITHM WILL FAIL WHEN ANY OF THE LEADING PRINCIPAL
!        SUBMATRICES ARE SINGULAR OR WHEN THE MATRIX ITSELF IS
!        SINGULAR.  WITH INTERCHANGES THIS ALGORITHM WILL FAIL ONLY
!        WHEN THE MATRIX ITSELF IS SINGULAR.  THE LEADING PRINCIPAL
!
!                                   [A B C]
!        SUBMATRICES OF THE MATRIX  [D E F]  ARE  [A]  AND  [A B] .
!                                   [G H I]                 [D E]
!
!   SUBPROGRAMS CALLED: -NONE-
!
!   CURRENT VERSION COMPLETED JANUARY 15, 1987
!
!   REFERENCE: STEWART, G.W., 'INTRODUCTION TO MATRIX COMPUTATIONS',
!              ACADEMIC PRESS, INC., 1973
!  --------------------------------------------------------------------
!   DEFINITION OF PASSED PARAMETERS
!
!     * A = MATRIX (SIZE NXN) TO BE INVERTED (REAL)
!
!   * LDA = LEADING DIMENSION OF MATRIX A [LDA>=N] (INTEGER)
!
!     * N = NUMBER OF ROWS AND COLUMNS OF MATRIX A (INTEGER)
!
!   IFLAG = ERROR INDICATOR ON OUTPUT (INTEGER)   INTERPRETATION:
!           -2 -> TOO MANY ROW INTERCHANGES NEEDED - INCREASE MX
!           -1 -> N>LDA
!            0 -> NO ERRORS DETECTED
!            K -> MATRIX A FOUND TO BE SINGULAR AT THE KTH STEP OF
!                 THE CROUT REDUCTION (1<=K<=N)
!
!   * INDICATES PARAMETERS REQUIRING INPUT VALUES
!  --------------------------------------------------------------------
!
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (MX=300)
      DIMENSION A(LDA,*),IEX(MX,2)
      IFLAG = 0
!
!   CHECK CONSISTENCY OF PASSED PARAMETERS
!
      IF (N.GT.LDA) THEN
         IFLAG = -1
         RETURN
      ENDIF
!
!   COMPUTE A = LU BY THE CROUT REDUCTION WHERE L IS LOWER TRIANGULAR
!   AND U IS UNIT UPPER TRIANGULAR (ALGORITHM 3.4, P. 138 OF THE
!   REFERENCE)
!
      NEX = 0
      DO 70 K = 1, N
         DO 20 I = K, N
            S = A(I,K)
            DO 10 L = 1, K-1
               S = S-A(I,L)*A(L,K)
   10       CONTINUE
            A(I,K) = S
   20    CONTINUE
!
!   INTERCHANGE ROWS IF NECESSARY
!
         Q = 0.0
         L = 0
         DO 30 I = K, N
            R = ABS(A(I,K))
            IF (R.GT.Q) THEN
               Q = R
               L = I
            ENDIF
   30    CONTINUE
         IF (L.EQ.0) THEN
            IFLAG = K
            RETURN
         ENDIF
         IF (L.NE.K) THEN
            NEX = NEX+1
            IF (NEX.GT.MX) THEN
               IFLAG = -2
               RETURN
            ENDIF
            IEX(NEX,1) = K
            IEX(NEX,2) = L
            DO 40 J = 1, N
               Q = A(K,J)
               A(K,J) = A(L,J)
               A(L,J) = Q
   40       CONTINUE
         ENDIF
!
!   END ROW INTERCHANGE SECTION
!
         DO 60 J = K+1, N
            S = A(K,J)
            DO 50 L = 1, K-1
               S = S-A(K,L)*A(L,J)
   50       CONTINUE
            A(K,J) = S/A(K,K)
   60    CONTINUE
   70 CONTINUE
!
!   INVERT THE LOWER TRIANGLE L IN PLACE (SIMILAR TO ALGORITHM 1.5,
!   P. 110 OF THE REFERENCE)
!
      DO 100 K = N, 1, -1
         A(K,K) = 1.0/A(K,K)
         DO 90 I = K-1, 1, -1
            S = 0.0
            DO 80 J = I+1, K
               S = S+A(J,I)*A(K,J)
   80       CONTINUE
            A(K,I) = -S/A(I,I)
   90    CONTINUE
  100 CONTINUE
!
!   INVERT THE UPPER TRIANGLE U IN PLACE (ALGORITHM 1.5, P. 110 OF
!   THE REFERENCE)
!
      DO 130 K = N, 1, -1
         DO 120 I = K-1, 1, -1
            S = A(I,K)
            DO 110 J = I+1, K-1
               S = S+A(I,J)*A(J,K)
  110       CONTINUE
            A(I,K) = -S
  120    CONTINUE
  130 CONTINUE
!
!   COMPUTE INV(A) = INV(U)*INV(L)
!
      DO 160 I = 1, N
         DO 150 J = 1, N
            IF (J.GT.I) THEN
               S = 0.0
               L = J
            ELSE
               S = A(I,J)
               L = I+1
            ENDIF
            DO 140 K = L, N
               S = S+A(I,K)*A(K,J)
  140       CONTINUE
            A(I,J) = S
  150    CONTINUE
  160 CONTINUE
!
!   INTERCHANGE COLUMNS OF INV(A) TO REVERSE EFFECT OF ROW
!   INTERCHANGES OF A
!
      DO 180 I = NEX, 1, -1
         K = IEX(I,1)
         L = IEX(I,2)
         DO 170 J = 1, N
            Q = A(J,K)
            A(J,K) = A(J,L)
            A(J,L) = Q
  170    CONTINUE
  180 CONTINUE
      RETURN
      END

