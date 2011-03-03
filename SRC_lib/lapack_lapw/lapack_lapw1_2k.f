      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
!
!     End of DLACPY
!
      END
      DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A.
!
!  Description
!  ===========
!
!  DLANSY returns the value
!
!     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANSY as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is to be referenced.
!          = 'U':  Upper triangular part of A is referenced
!          = 'L':  Lower triangular part of A is referenced
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
!          set to zero.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
               ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL DLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      DLANSY = VALUE
      RETURN
!
!     End of DLANSY
!
      END
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                         T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                         WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) DOUBLE PRECISION array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) &
         RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C1'
!
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2
!
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C2'
!
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1
!
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C1'
!
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2'
!
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                              WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2' * W'
!
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C2'
!
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1'
!
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1' * W'
!
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of DLARFB
!
      END
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) DOUBLE PRECISION array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   VII
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
!
!              general case
!
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!
!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
                  CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ), &
                              V( I, 1 ), LDV, V( I, I ), 1, ZERO, &
                              T( 1, I ), 1 )
               ELSE
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
                  CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ), &
                              V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                              T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                           LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                     CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ), &
                                 V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO, &
                                 T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                     CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                                 V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, &
                                 T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                              T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
!
!     End of DLARFT
!
      END
      SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
!     ..
!
!  Purpose
!  =======
!
!  DLATRD reduces NB rows and columns of a real symmetric matrix A to
!  symmetric tridiagonal form by an orthogonal similarity
!  transformation Q' * A * Q, and returns the matrices V and W which are
!  needed to apply the transformation to the unreduced part of A.
!
!  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
!  matrix, of which the upper triangle is supplied;
!  if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
!  matrix, of which the lower triangle is supplied.
!
!  This is an auxiliary routine called by DSYTRD.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U': Upper triangular
!          = 'L': Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  NB      (input) INTEGER
!          The number of rows and columns to be reduced.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit:
!          if UPLO = 'U', the last NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements above the diagonal
!            with the array TAU, represent the orthogonal matrix Q as a
!            product of elementary reflectors;
!          if UPLO = 'L', the first NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements below the diagonal
!            with the array TAU, represent the  orthogonal matrix Q as a
!            product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= (1,N).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!          elements of the last NB columns of the reduced matrix;
!          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!          the first NB columns of the reduced matrix.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors, stored in
!          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!          See Further Details.
!
!  W       (output) DOUBLE PRECISION array, dimension (LDW,NB)
!          The n-by-nb matrix W required to update the unreduced part
!          of A.
!
!  LDW     (input) INTEGER
!          The leading dimension of the array W. LDW >= max(1,N).
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n) H(n-1) . . . H(n-nb+1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!  and tau in TAU(i-1).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!  and tau in TAU(i).
!
!  The elements of the vectors v together form the n-by-nb matrix V
!  which is needed, with W, to apply the transformation to the unreduced
!  part of the matrix, using a symmetric rank-2k update of the form:
!  A := A - V*W' - W*V'.
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5 and nb = 2:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  a   a   a   v4  v5 )              (  d                  )
!    (      a   a   v4  v5 )              (  1   d              )
!    (          a   1   v5 )              (  v1  1   a          )
!    (              d   1  )              (  v1  v2  a   a      )
!    (                  d  )              (  v1  v2  a   a   a  )
!
!  where d denotes a diagonal element of the reduced matrix, a denotes
!  an element of the original matrix that is unchanged, and vi denotes
!  an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IW
      DOUBLE PRECISION   ALPHA
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Reduce last NB columns of upper triangle
!
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
!
!              Update A(1:i,i)
!
               CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ), &
                           LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), &
                           LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)
!
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
!
!              Compute W(1:i-1,i)
!
               CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1, &
                           ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ), &
                              LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE, &
                              A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, &
                              W( 1, IW ), 1 )
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ), &
                              LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE, &
                              W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, &
                              W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1, &
                       A( 1, I ), 1 )
               CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
!
   10    CONTINUE
      ELSE
!
!        Reduce first NB columns of lower triangle
!
         DO 20 I = 1, NB
!
!           Update A(i:n,i)
!
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ), &
                        LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ), &
                        LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)
!
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, &
                            TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
!
!              Compute W(i+1:n,i)
!
               CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA, &
                           A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, &
                           A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ), &
                           LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, &
                           A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ), &
                           LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1, &
                       A( I+1, I ), 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
!
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of DLATRD
!
      END
      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DORGQL generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the last N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, NB, NBMIN, &
                         NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.
!
      NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.
!
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
!
!        Set A(m-kk+1:m,1:n-kk) to zero.
!
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the first or only block.
!
      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB, &
                            A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
               CALL DLARFB( 'Left', 'No transpose', 'Backward', &
                            'Columnwise', M-K+I+IB-1, N-K+I-1, IB, &
                            A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, &
                            WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows 1:m-k+i+ib-1 of current block
!
            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA, &
                         TAU( I ), WORK, IINFO )
!
!           Set rows m-k+i+ib:m of current block to zero
!
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORGQL
!
      END
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, NB, &
                         NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.
!
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.N ) &
         CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                      TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                            A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL DLARFB( 'Left', 'No transpose', 'Forward', &
                            'Columnwise', M-I+1, N-I-IB+1, IB, &
                            A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                            LDA, WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORGQR
!
      END
      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DORGTR generates a real orthogonal matrix Q which is defined as the
!  product of n-1 elementary reflectors of order N, as returned by
!  DSYTRD:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from DSYTRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from DSYTRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by DSYTRD.
!          On exit, the N-by-N orthogonal matrix Q.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSYTRD.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N-1).
!          For optimum performance LWORK >= (N-1)*NB, where NB is
!          the optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORGQL, DORGQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGTR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to DSYTRD with UPLO = 'U'
!
!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix
!
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
!
      ELSE
!
!        Q was determined by a call to DSYTRD with UPLO = 'L'.
!
!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix
!
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
            CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, &
                         LWORK, IINFO )
         END IF
      END IF
      RETURN
!
!     End of DORGTR
!
      END
      SUBROUTINE DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORM2L overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQLF in the last k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, MI, NI, NQ
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2L', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) &
           THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(1:m-k+i,1:n)
!
            MI = M - K + I
         ELSE
!
!           H(i) is applied to C(1:m,1:n-k+i)
!
            NI = N - K + I
         END IF
!
!        Apply H(i)
!
         AII = A( NQ-K+I, I )
         A( NQ-K+I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( 1, I ), 1, TAU( I ), C, LDC, &
                     WORK )
         A( NQ-K+I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of DORM2L
!
      END
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORM2R overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) &
           THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i)
!
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), &
                     LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of DORM2R
!
      END
      SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), &
                         WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DORMQL overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQLF in the last k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IINFO, IWS, LDWORK, MI, NB, &
                         NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2L, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.  NB may be at most NBMAX, where NBMAX
!     is used to define the local array T.
!
      NB = MIN( NBMAX, ILAENV( 1, 'DORMQL', SIDE // TRANS, M, N, K, &
           -1 ) )
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQL', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL DLARFT( 'Backward', 'Columnwise', NQ-K+I+IB-1, IB, &
                         A( 1, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H' is applied to C(1:m-k+i+ib-1,1:n)
!
               MI = M - K + I + IB - 1
            ELSE
!
!              H or H' is applied to C(1:m,1:n-k+i+ib-1)
!
               NI = N - K + I + IB - 1
            END IF
!
!           Apply H or H'
!
            CALL DLARFB( SIDE, TRANS, 'Backward', 'Columnwise', MI, NI, &
                         IB, A( 1, I ), LDA, T, LDT, C, LDC, WORK, &
                         LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORMQL
!
      END
      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), &
                         WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DORMQR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
                         MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.  NB may be at most NBMAX, where NBMAX
!     is used to define the local array T.
!
      NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K, &
           -1 ) )
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
                         LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H' is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H' is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H'
!
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                         IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, &
                         WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORMQR
!
      END
      SUBROUTINE DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), &
                         WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DORMTR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by DSYTRD:
!
!  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from DSYTRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from DSYTRD.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                               (LDA,M) if SIDE = 'L'
!                               (LDA,N) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by DSYTRD.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!
!  TAU     (input) DOUBLE PRECISION array, dimension
!                               (M-1) if SIDE = 'L'
!                               (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSYTRD.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LEFT, UPPER
      INTEGER            I1, I2, IINFO, MI, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORMQL, DORMQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) ) &
                THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMTR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( LEFT ) THEN
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      END IF
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to DSYTRD with UPLO = 'U'
!
         CALL DORMQL( SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C, &
                      LDC, WORK, LWORK, IINFO )
      ELSE
!
!        Q was determined by a call to DSYTRD with UPLO = 'L'
!
         IF( LEFT ) THEN
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         END IF
         CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, &
                      C( I1, I2 ), LDC, WORK, LWORK, IINFO )
      END IF
      RETURN
!
!     End of DORMTR
!
      END
      SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
                         ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, &
                         IFAIL, INFO )
!
!  -- LAPACK driver routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYEVX computes selected eigenvalues and, optionally, eigenvectors
!  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
!  selected by specifying either a range of values or a range of indices
!  for the desired eigenvalues.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  RANGE   (input) CHARACTER*1
!          = 'A': all eigenvalues will be found.
!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!                 will be found.
!          = 'I': the IL-th through IU-th eigenvalues will be found.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, the lower triangle (if UPLO='L') or the upper
!          triangle (if UPLO='U') of A, including the diagonal, is
!          destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  VL      (input) DOUBLE PRECISION
!  VU      (input) DOUBLE PRECISION
!          If RANGE='V', the lower and upper bounds of the interval to
!          be searched for eigenvalues. VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!
!  IL      (input) INTEGER
!  IU      (input) INTEGER
!          If RANGE='I', the indices (in ascending order) of the
!          smallest and largest eigenvalues to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The absolute error tolerance for the eigenvalues.
!          An approximate eigenvalue is accepted as converged
!          when it is determined to lie in an interval [a,b]
!          of width less than or equal to
!
!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!
!          where EPS is the machine precision.  If ABSTOL is less than
!          or equal to zero, then  EPS*|T|  will be used in its place,
!          where |T| is the 1-norm of the tridiagonal matrix obtained
!          by reducing A to tridiagonal form.
!
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!          If this routine returns with INFO>0, indicating that some
!          eigenvectors did not converge, try setting ABSTOL to
!          2*DLAMCH('S').
!
!          See "Computing Small Singular Values of Bidiagonal Matrices
!          with Guaranteed High Relative Accuracy," by Demmel and
!          Kahan, LAPACK Working Note #3.
!
!  M       (output) INTEGER
!          The total number of eigenvalues found.  0 <= M <= N.
!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          On normal exit, the first M elements contain the selected
!          eigenvalues in ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!          contain the orthonormal eigenvectors of the matrix A
!          corresponding to the selected eigenvalues, with the i-th
!          column of Z holding the eigenvector associated with W(i).
!          If an eigenvector fails to converge, then that column of Z
!          contains the latest approximation to the eigenvector, and the
!          index of the eigenvector is returned in IFAIL.
!          If JOBZ = 'N', then Z is not referenced.
!          Note: the user must ensure that at least max(1,M) columns are
!          supplied in the array Z; if RANGE = 'V', the exact value of M
!          is not known in advance and an upper bound must be used.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,8*N).
!          For optimal efficiency, LWORK >= (NB+3)*N,
!          where NB is the blocksize for DSYTRD returned by ILAENV.
!
!  IWORK   (workspace) INTEGER array, dimension (5*N)
!
!  IFAIL   (output) INTEGER array, dimension (N)
!          If JOBZ = 'V', then if INFO = 0, the first M elements of
!          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!          indices of the eigenvectors that failed to converge.
!          If JOBZ = 'N', then IFAIL is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, then i eigenvectors failed to converge.
!                Their indices are stored in array IFAIL.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, &
                         INDISP, INDIWO, INDTAU, INDWKN, INDWRK, ISCALE, &
                         ITMP1, J, JJ, LLWORK, LLWRKN, LOPT, NSPLIT
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, &
                         SIGMA, SMLNUM, TMP1, VLL, VUU
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, DLAMCH, DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DORGTR, DORMTR, DSCAL, DSTEBZ, &
                         DSTEIN, DSTEQR, DSTERF, DSWAP, DSYTRD, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
!
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
         INFO = -8
      ELSE IF( INDEIG .AND. IL.LT.1 ) THEN
         INFO = -9
      ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
         INFO = -10
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LWORK.LT.MAX( 1, 8*N ) ) THEN
         INFO = -17
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEVX', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      M = 0
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( N.EQ.1 ) THEN
         WORK( 1 ) = 7
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = A( 1, 1 )
         ELSE
            IF( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) THEN
               M = 1
               W( 1 ) = A( 1, 1 )
            END IF
         END IF
         IF( WANTZ ) &
            Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     Get machine constants.
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
!
!     Scale matrix to allowable range, if necessary.
!
      ISCALE = 0
      ABSTLL = ABSTOL
      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      END IF
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N
               CALL DSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL DSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 ) &
            ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
!
!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
!
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDWRK = INDD + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), &
                   WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
      LOPT = 3*N + WORK( INDWRK )
!
!     If all eigenvalues are desired and ABSTOL is less than or equal to
!     zero, then call DSTERF or DORGTR and SSTEQR.  If this fails for
!     some eigenvalue, then try DSTEBZ.
!
      IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND. &
          ( ABSTOL.LE.ZERO ) ) THEN
         CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         IF( .NOT.WANTZ ) THEN
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL DLACPY( 'A', N, N, A, LDA, Z, LDZ )
            CALL DORGTR( UPLO, N, Z, LDZ, WORK( INDTAU ), &
                         WORK( INDWRK ), LLWORK, IINFO )
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ, &
                         WORK( INDWRK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 30 I = 1, N
                  IFAIL( I ) = 0
   30          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 40
         END IF
         INFO = 0
      END IF
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.
!
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, &
                   WORK( INDD ), WORK( INDE ), M, NSPLIT, W, &
                   IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ), &
                   IWORK( INDIWO ), INFO )
!
      IF( WANTZ ) THEN
         CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), M, W, &
                      IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, &
                      WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         CALL DORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, &
                      LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
   40 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
      IF( WANTZ ) THEN
         DO 60 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 50 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   50       CONTINUE
!
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   60    CONTINUE
      END IF
!
!     Set WORK(1) to optimal workspace size.
!
      WORK( 1 ) = MAX( 7*N, LOPT )
!
      RETURN
!
!     End of DSYEVX
!
      END
      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
!  form T by an orthogonal similarity transformation: Q' * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, &
                         HALF = 1.0D0 / 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   ALPHA, TAUI
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A
!
         DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               A( I, I+1 ) = ONE
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
               CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, &
                           TAU, 1 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, &
                           LDA )
!
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
!
!        Reduce the lower triangle of A
!
         DO 20 I = 1, N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, &
                         TAUI )
            E( I ) = A( I+1, I )
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               A( I+1, I ) = ONE
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
               CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, &
                           A( I+1, I ), 1, ZERO, TAU( I ), 1 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ), &
                       1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, &
                           A( I+1, I+1 ), LDA )
!
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
!
      RETURN
!
!     End of DSYTD2
!
      END
      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRD reduces a real symmetric matrix A to real symmetric
!  tridiagonal form T by an orthogonal similarity transformation:
!  Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= 1.
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.
!
      NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
         NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN ) &
                  NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, &
                         LDWORK )
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A - V*W' - W*V'
!
            CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ), &
                         LDA, WORK, LDWORK, ONE, A, LDA )
!
!           Copy superdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
!
!        Reduce the lower triangle of A
!
         DO 40 I = 1, N - NX, NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), &
                         TAU( I ), WORK, LDWORK )
!
!           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
!           an update of the form:  A := A - V*W' - W*V'
!
            CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE, &
                         A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, &
                         A( I+NB, I+NB ), LDA )
!
!           Copy subdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ), &
                      TAU( I ), IINFO )
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DSYTRD
!
      END
      SUBROUTINE ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
                         ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, &
                         IWORK, IFAIL, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  ZHEEVX computes selected eigenvalues and, optionally, eigenvectors
!  of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can
!  be selected by specifying either a range of values or a range of
!  indices for the desired eigenvalues.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  RANGE   (input) CHARACTER*1
!          = 'A': all eigenvalues will be found.
!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!                 will be found.
!          = 'I': the IL-th through IU-th eigenvalues will be found.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, the lower triangle (if UPLO='L') or the upper
!          triangle (if UPLO='U') of A, including the diagonal, is
!          destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  VL      (input) DOUBLE PRECISION
!  VU      (input) DOUBLE PRECISION
!          If RANGE='V', the lower and upper bounds of the interval to
!          be searched for eigenvalues. VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!
!  IL      (input) INTEGER
!  IU      (input) INTEGER
!          If RANGE='I', the indices (in ascending order) of the
!          smallest and largest eigenvalues to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The absolute error tolerance for the eigenvalues.
!          An approximate eigenvalue is accepted as converged
!          when it is determined to lie in an interval [a,b]
!          of width less than or equal to
!
!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!
!          where EPS is the machine precision.  If ABSTOL is less than
!          or equal to zero, then  EPS*|T|  will be used in its place,
!          where |T| is the 1-norm of the tridiagonal matrix obtained
!          by reducing A to tridiagonal form.
!
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!          If this routine returns with INFO>0, indicating that some
!          eigenvectors did not converge, try setting ABSTOL to
!          2*DLAMCH('S').
!
!          See "Computing Small Singular Values of Bidiagonal Matrices
!          with Guaranteed High Relative Accuracy," by Demmel and
!          Kahan, LAPACK Working Note #3.
!
!  M       (output) INTEGER
!          The total number of eigenvalues found.  0 <= M <= N.
!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          On normal exit, the first M elements contain the selected
!          eigenvalues in ascending order.
!
!  Z       (output) COMPLEX*16 array, dimension (LDZ, max(1,M))
!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!          contain the orthonormal eigenvectors of the matrix A
!          corresponding to the selected eigenvalues, with the i-th
!          column of Z holding the eigenvector associated with W(i).
!          If an eigenvector fails to converge, then that column of Z
!          contains the latest approximation to the eigenvector, and the
!          index of the eigenvector is returned in IFAIL.
!          If JOBZ = 'N', then Z is not referenced.
!          Note: the user must ensure that at least max(1,M) columns are
!          supplied in the array Z; if RANGE = 'V', the exact value of M
!          is not known in advance and an upper bound must be used.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,2*N-1).
!          For optimal efficiency, LWORK >= (NB+1)*N,
!          where NB is the max of the blocksize for ZHETRD and for
!          ZUNMTR as returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  RWORK   (workspace) DOUBLE PRECISION array, dimension (7*N)
!
!  IWORK   (workspace) INTEGER array, dimension (5*N)
!
!  IFAIL   (output) INTEGER array, dimension (N)
!          If JOBZ = 'V', then if INFO = 0, the first M elements of
!          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!          indices of the eigenvectors that failed to converge.
!          If JOBZ = 'N', then IFAIL is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, then i eigenvectors failed to converge.
!                Their indices are stored in array IFAIL.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, LQUERY, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, &
                         INDISP, INDIWK, INDRWK, INDTAU, INDWRK, ISCALE, &
                         ITMP1, J, JJ, LLWORK, LOPT, LWKOPT, NB, NSPLIT
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, &
                         SIGMA, SMLNUM, TMP1, VLL, VUU
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANHE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, &
                         ZHETRD, ZLACPY, ZSTEIN, ZSTEQR, ZSWAP, ZUNGTR, &
                         ZUNMTR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LQUERY = ( LWORK.EQ.-1 )
!
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL ) &
               INFO = -8
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -9
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -10
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
            INFO = -15
         ELSE IF( LWORK.LT.MAX( 1, 2*N-1 ) .AND. .NOT.LQUERY ) THEN
            INFO = -17
         END IF
      END IF
!
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         NB = MAX( NB, ILAENV( 1, 'ZUNMTR', UPLO, N, -1, -1, -1 ) )
         LWKOPT = ( NB+1 )*N
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      M = 0
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( N.EQ.1 ) THEN
         WORK( 1 ) = 1
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = A( 1, 1 )
         ELSE IF( VALEIG ) THEN
            IF( VL.LT.DBLE( A( 1, 1 ) ) .AND. VU.GE.DBLE( A( 1, 1 ) ) ) &
                 THEN
               M = 1
               W( 1 ) = A( 1, 1 )
            END IF
         END IF
         IF( WANTZ ) &
            Z( 1, 1 ) = CONE
         RETURN
      END IF
!
!     Get machine constants.
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
!
!     Scale matrix to allowable range, if necessary.
!
      ISCALE = 0
      ABSTLL = ABSTOL
      VLL = VL
      VUU = VU
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N
               CALL ZDSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL ZDSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 ) &
            ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
!
!     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
!
      INDD = 1
      INDE = INDD + N
      INDRWK = INDE + N
      INDTAU = 1
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL ZHETRD( UPLO, N, A, LDA, RWORK( INDD ), RWORK( INDE ), &
                   WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
      LOPT = N + WORK( INDWRK )
!
!     If all eigenvalues are desired and ABSTOL is less than or equal to
!     zero, then call DSTERF or ZUNGTR and ZSTEQR.  If this fails for
!     some eigenvalue, then try DSTEBZ.
!
      IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND. &
          ( ABSTOL.LE.ZERO ) ) THEN
         CALL DCOPY( N, RWORK( INDD ), 1, W, 1 )
         INDEE = INDRWK + 2*N
         IF( .NOT.WANTZ ) THEN
            CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL DSTERF( N, W, RWORK( INDEE ), INFO )
         ELSE
            CALL ZLACPY( 'A', N, N, A, LDA, Z, LDZ )
            CALL ZUNGTR( UPLO, N, Z, LDZ, WORK( INDTAU ), &
                         WORK( INDWRK ), LLWORK, IINFO )
            CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL ZSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ, &
                         RWORK( INDRWK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 30 I = 1, N
                  IFAIL( I ) = 0
   30          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 40
         END IF
         INFO = 0
      END IF
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
!
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWK = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, &
                   RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, &
                   IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), &
                   IWORK( INDIWK ), INFO )
!
      IF( WANTZ ) THEN
         CALL ZSTEIN( N, RWORK( INDD ), RWORK( INDE ), M, W, &
                      IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, &
                      RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO )
!
!        Apply unitary matrix used in reduction to tridiagonal
!        form to eigenvectors returned by ZSTEIN.
!
         CALL ZUNMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, &
                      LDZ, WORK( INDWRK ), LLWORK, IINFO )
      END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
   40 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
      IF( WANTZ ) THEN
         DO 60 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 50 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   50       CONTINUE
!
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   60    CONTINUE
      END IF
!
!     Set WORK(1) to optimal complex workspace size.
!
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of ZHEEVX
!
      END
      SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  ZHETD2 reduces a complex Hermitian matrix A to real symmetric
!  tridiagonal form T by a unitary similarity transformation:
!  Q' * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          Hermitian matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the unitary
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the unitary matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) COMPLEX*16 array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ), &
                         HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      COMPLEX*16         ALPHA, TAUI
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZHEMV, ZHER2, ZLARFG
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETD2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A
!
         A( N, N ) = DBLE( A( N, N ) )
         DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
            ALPHA = A( I, I+1 )
            CALL ZLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               A( I, I+1 ) = ONE
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
               CALL ZHEMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, &
                           TAU, 1 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL ZAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL ZHER2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, &
                           LDA )
!
            END IF
            A( I, I+1 ) = E( I )
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
!
!        Reduce the lower triangle of A
!
         A( 1, 1 ) = DBLE( A( 1, 1 ) )
         DO 20 I = 1, N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            ALPHA = A( I+1, I )
            CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               A( I+1, I ) = ONE
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
               CALL ZHEMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, &
                           A( I+1, I ), 1, ZERO, TAU( I ), 1 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, A( I+1, I ), &
                       1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL ZHER2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, &
                           A( I+1, I+1 ), LDA )
!
            END IF
            A( I+1, I ) = E( I )
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
!
      RETURN
!
!     End of ZHETD2
!
      END
      SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZHETRD reduces a complex Hermitian matrix A to real symmetric
!  tridiagonal form T by a unitary similarity transformation:
!  Q**H * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the unitary
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the unitary matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) COMPLEX*16 array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= 1.
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, &
                         NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHER2K, ZHETD2, ZLATRD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.
!
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
         NX = MAX( NB, ILAENV( 3, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'ZHETRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN ) &
                  NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL ZLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, &
                         LDWORK )
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A - V*W' - W*V'
!
            CALL ZHER2K( UPLO, 'No transpose', I-1, NB, -CONE, &
                         A( 1, I ), LDA, WORK, LDWORK, ONE, A, LDA )
!
!           Copy superdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL ZHETD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
!
!        Reduce the lower triangle of A
!
         DO 40 I = 1, N - NX, NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL ZLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), &
                         TAU( I ), WORK, LDWORK )
!
!           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
!           an update of the form:  A := A - V*W' - W*V'
!
            CALL ZHER2K( UPLO, 'No transpose', N-I-NB+1, NB, -CONE, &
                         A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, &
                         A( I+NB, I+NB ), LDA )
!
!           Copy subdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL ZHETD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ), &
                      TAU( I ), IINFO )
      END IF
!
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZHETRD
!
      END
      SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
!          is accessed; if UPLO = 'L', only the lower trapezium is
!          accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) COMPLEX*16 array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
!
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLACPY
!
      END
      DOUBLE PRECISION FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLANHE  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  complex hermitian matrix A.
!
!  Description
!  ===========
!
!  ZLANHE returns the value
!
!     ZLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in ZLANHE as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          hermitian matrix A is to be referenced.
!          = 'U':  Upper triangular part of A is referenced
!          = 'L':  Lower triangular part of A is referenced
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, ZLANHE is
!          set to zero.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The hermitian matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced. Note that the imaginary parts of the diagonal
!          elements need not be set and are assumed to be zero.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J - 1
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
               VALUE = MAX( VALUE, ABS( DBLE( A( J, J ) ) ) )
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( DBLE( A( J, J ) ) ) )
               DO 30 I = J + 1, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
               ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is hermitian).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( A( J, J ) ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( A( J, J ) ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         DO 130 I = 1, N
            IF( DBLE( A( I, I ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( A( I, I ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      ZLANHE = VALUE
      RETURN
!
!     End of ZLANHE
!
      END
      DOUBLE PRECISION FUNCTION ZLANSY( NORM, UPLO, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLANSY  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  complex symmetric matrix A.
!
!  Description
!  ===========
!
!  ZLANSY returns the value
!
!     ZLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in ZLANSY as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is to be referenced.
!          = 'U':  Upper triangular part of A is referenced
!          = 'L':  Lower triangular part of A is referenced
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, ZLANSY is
!          set to zero.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
               ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL ZLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      ZLANSY = VALUE
      RETURN
!
!     End of ZLANSY
!
      END
      SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                         T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                         WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLARFB applies a complex block reflector H or its transpose H' to a
!  complex M-by-N matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'C': apply H' (Conjugate transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) COMPLEX*16 array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) COMPLEX*16 array, dimension (LDT,K)
!          The triangular K-by-K matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMM, ZLACGV, ZTRMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) &
         RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C1'
!
               DO 10 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2
!
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, &
                              K, M-K, ONE, C( K+1, 1 ), LDC, &
                              V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W'
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', &
                              M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK, &
                              LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                           'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2'
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
                              N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ), &
                              LDV, ONE, C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                           'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C2'
!
               DO 70 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1
!
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, &
                              K, M-K, ONE, C, LDC, V, LDV, ONE, WORK, &
                              LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W'
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', &
                              M-K, N, K, -ONE, V, LDV, WORK, LDWORK, &
                              ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                           'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK, &
                           LDWORK )
!
!              C2 := C2 - W'
!
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - &
                                     DCONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1'
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
                              N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE, &
                              C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                           'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK, &
                           LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C1'
!
               DO 130 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1'
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                           'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2'
!
                  CALL ZGEMM( 'Conjugate transpose', &
                              'Conjugate transpose', N, K, M-K, ONE, &
                              C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                              WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2' * W'
!
                  CALL ZGEMM( 'Conjugate transpose', &
                              'Conjugate transpose', M-K, N, K, -ONE, &
                              V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1'
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
                           'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2'
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
                              K, N-K, ONE, C( 1, K+1 ), LDC, &
                              V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C2'
!
               DO 190 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2'
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                           'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK, &
                           LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1'
!
                  CALL ZGEMM( 'Conjugate transpose', &
                              'Conjugate transpose', N, K, M-K, ONE, C, &
                              LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1' * W'
!
                  CALL ZGEMM( 'Conjugate transpose', &
                              'Conjugate transpose', M-K, N, K, -ONE, V, &
                              LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - &
                                     DCONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2'
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
                           'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK, &
                           LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1'
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
                              K, N-K, ONE, C, LDC, V, LDV, ONE, WORK, &
                              LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of ZLARFB
!
      END
      SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLARFT forms the triangular factor T of a complex block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) COMPLEX*16 array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) COMPLEX*16 array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX*16         VII
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZLACGV, ZTRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
!
!              general case
!
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!
!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
                  CALL ZGEMV( 'Conjugate transpose', N-I+1, I-1, &
                              -TAU( I ), V( I, 1 ), LDV, V( I, I ), 1, &
                              ZERO, T( 1, I ), 1 )
               ELSE
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
                  IF( I.LT.N ) &
                     CALL ZLACGV( N-I, V( I, I+1 ), LDV )
                  CALL ZGEMV( 'No transpose', I-1, N-I+1, -TAU( I ), &
                              V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                              T( 1, I ), 1 )
                  IF( I.LT.N ) &
                     CALL ZLACGV( N-I, V( I, I+1 ), LDV )
               END IF
               V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL ZTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                           LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                     CALL ZGEMV( 'Conjugate transpose', N-K+I, K-I, &
                                 -TAU( I ), V( 1, I+1 ), LDV, V( 1, I ), &
                                 1, ZERO, T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                     CALL ZLACGV( N-K+I-1, V( I, 1 ), LDV )
                     CALL ZGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                                 V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, &
                                 T( I+1, I ), 1 )
                     CALL ZLACGV( N-K+I-1, V( I, 1 ), LDV )
                     V( I, N-K+I ) = VII
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL ZTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                              T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
!
!     End of ZLARFT
!
      END
      SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to
!  Hermitian tridiagonal form by a unitary similarity
!  transformation Q' * A * Q, and returns the matrices V and W which are
!  needed to apply the transformation to the unreduced part of A.
!
!  If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a
!  matrix, of which the upper triangle is supplied;
!  if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a
!  matrix, of which the lower triangle is supplied.
!
!  This is an auxiliary routine called by ZHETRD.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER
!          Specifies whether the upper or lower triangular part of the
!          Hermitian matrix A is stored:
!          = 'U': Upper triangular
!          = 'L': Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  NB      (input) INTEGER
!          The number of rows and columns to be reduced.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit:
!          if UPLO = 'U', the last NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements above the diagonal
!            with the array TAU, represent the unitary matrix Q as a
!            product of elementary reflectors;
!          if UPLO = 'L', the first NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements below the diagonal
!            with the array TAU, represent the  unitary matrix Q as a
!            product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!          elements of the last NB columns of the reduced matrix;
!          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!          the first NB columns of the reduced matrix.
!
!  TAU     (output) COMPLEX*16 array, dimension (N-1)
!          The scalar factors of the elementary reflectors, stored in
!          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!          See Further Details.
!
!  W       (output) COMPLEX*16 array, dimension (LDW,NB)
!          The n-by-nb matrix W required to update the unreduced part
!          of A.
!
!  LDW     (input) INTEGER
!          The leading dimension of the array W. LDW >= max(1,N).
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n) H(n-1) . . . H(n-nb+1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!  and tau in TAU(i-1).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!  and tau in TAU(i).
!
!  The elements of the vectors v together form the n-by-nb matrix V
!  which is needed, with W, to apply the transformation to the unreduced
!  part of the matrix, using a Hermitian rank-2k update of the form:
!  A := A - V*W' - W*V'.
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5 and nb = 2:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  a   a   a   v4  v5 )              (  d                  )
!    (      a   a   v4  v5 )              (  1   d              )
!    (          a   1   v5 )              (  v1  1   a          )
!    (              d   1  )              (  v1  v2  a   a      )
!    (                  d  )              (  v1  v2  a   a   a  )
!
!  where d denotes a diagonal element of the reduced matrix, a denotes
!  an element of the original matrix that is unchanged, and vi denotes
!  an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE, HALF
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ), &
                         HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IW
      COMPLEX*16         ALPHA
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEMV, ZHEMV, ZLACGV, ZLARFG, ZSCAL
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Reduce last NB columns of upper triangle
!
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
!
!              Update A(1:i,i)
!
               A( I, I ) = DBLE( A( I, I ) )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ), &
                           LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), &
                           LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)
!
               ALPHA = A( I-1, I )
               CALL ZLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = ALPHA
               A( I-1, I ) = ONE
!
!              Compute W(1:i-1,i)
!
               CALL ZHEMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1, &
                           ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE, &
                              W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO, &
                              W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, &
                              A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, &
                              W( 1, IW ), 1 )
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE, &
                              A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO, &
                              W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, &
                              W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, &
                              W( 1, IW ), 1 )
               END IF
               CALL ZSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*ZDOTC( I-1, W( 1, IW ), 1, &
                       A( 1, I ), 1 )
               CALL ZAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
!
   10    CONTINUE
      ELSE
!
!        Reduce first NB columns of lower triangle
!
         DO 20 I = 1, NB
!
!           Update A(i:n,i)
!
            A( I, I ) = DBLE( A( I, I ) )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ), &
                        LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ), &
                        LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            A( I, I ) = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)
!
               ALPHA = A( I+1, I )
               CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, &
                            TAU( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
!
!              Compute W(i+1:n,i)
!
               CALL ZHEMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA, &
                           A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, &
                           W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO, &
                           W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ), &
                           LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, &
                           A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, &
                           W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ), &
                           LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*ZDOTC( N-I, W( I+1, I ), 1, &
                       A( I+1, I ), 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
!
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLATRD
!
      END
      SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNGQL generates an M-by-N complex matrix Q with orthonormal columns,
!  which is defined as the last N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by ZGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by ZGEQLF in the last k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQLF.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, &
                         NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2L
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'ZUNGQL', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'ZUNGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.
!
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
!
!        Set A(m-kk+1:m,1:n-kk) to zero.
!
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the first or only block.
!
      CALL ZUNG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL ZLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB, &
                            A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
               CALL ZLARFB( 'Left', 'No transpose', 'Backward', &
                            'Columnwise', M-K+I+IB-1, N-K+I-1, IB, &
                            A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, &
                            WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows 1:m-k+i+ib-1 of current block
!
            CALL ZUNG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA, &
                         TAU( I ), WORK, IINFO )
!
!           Set rows m-k+i+ib:m of current block to zero
!
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of ZUNGQL
!
      END
      SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNGQR generates an M-by-N complex matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by ZGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by ZGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQRF.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
                         LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'ZUNGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'ZUNGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.N ) &
         CALL ZUNG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                      TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL ZLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                            A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL ZLARFB( 'Left', 'No transpose', 'Forward', &
                            'Columnwise', M-I+1, N-I-IB+1, IB, &
                            A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                            LDA, WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL ZUNG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of ZUNGQR
!
      END
      SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNGTR generates a complex unitary matrix Q which is defined as the
!  product of n-1 elementary reflectors of order N, as returned by
!  ZHETRD:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from ZHETRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from ZHETRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by ZHETRD.
!          On exit, the N-by-N unitary matrix Q.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= N.
!
!  TAU     (input) COMPLEX*16 array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZHETRD.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= N-1.
!          For optimum performance LWORK >= (N-1)*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNGQL, ZUNGQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'ZUNGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'ZUNGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to ZHETRD with UPLO = 'U'
!
!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix
!
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL ZUNGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
!
      ELSE
!
!        Q was determined by a call to ZHETRD with UPLO = 'L'.
!
!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix
!
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
            CALL ZUNGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, &
                         LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNGTR
!
      END
      SUBROUTINE ZUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNM2L overwrites the general complex m-by-n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'C', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'C',
!
!  where Q is a complex unitary matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by ZGEQLF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'C': apply Q' (Conjugate transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          ZGEQLF in the last k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQLF.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the m-by-n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) COMPLEX*16 array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, MI, NI, NQ
      COMPLEX*16         AII, TAUI
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNM2L', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      IF( ( LEFT .AND. NOTRAN .OR. .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) or H(i)' is applied to C(1:m-k+i,1:n)
!
            MI = M - K + I
         ELSE
!
!           H(i) or H(i)' is applied to C(1:m,1:n-k+i)
!
            NI = N - K + I
         END IF
!
!        Apply H(i) or H(i)'
!
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = DCONJG( TAU( I ) )
         END IF
         AII = A( NQ-K+I, I )
         A( NQ-K+I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( 1, I ), 1, TAUI, C, LDC, WORK )
         A( NQ-K+I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of ZUNM2L
!
      END
      SUBROUTINE ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNM2R overwrites the general complex m-by-n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'C', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'C',
!
!  where Q is a complex unitary matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by ZGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'C': apply Q' (Conjugate transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          ZGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQRF.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the m-by-n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) COMPLEX*16 array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) or H(i)' is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) or H(i)' is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i) or H(i)'
!
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = DCONJG( TAU( I ) )
         END IF
         AII = A( I, I )
         A( I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC, &
                     WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of ZUNM2R
!
      END
      SUBROUTINE ZUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNMQL overwrites the general complex M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'C':      Q**H * C       C * Q**H
!
!  where Q is a complex unitary matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by ZGEQLF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**H from the Left;
!          = 'R': apply Q or Q**H from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'C':  Transpose, apply Q**H.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          ZGEQLF in the last k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQLF.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IINFO, IWS, LDWORK, LWKOPT, &
                         MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      COMPLEX*16         T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNM2L
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         NB = MIN( NBMAX, ILAENV( 1, 'ZUNMQL', SIDE // TRANS, M, N, K, &
              -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMQL', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL ZUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL ZLARFT( 'Backward', 'Columnwise', NQ-K+I+IB-1, IB, &
                         A( 1, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H' is applied to C(1:m-k+i+ib-1,1:n)
!
               MI = M - K + I + IB - 1
            ELSE
!
!              H or H' is applied to C(1:m,1:n-k+i+ib-1)
!
               NI = N - K + I + IB - 1
            END IF
!
!           Apply H or H'
!
            CALL ZLARFB( SIDE, TRANS, 'Backward', 'Columnwise', MI, NI, &
                         IB, A( 1, I ), LDA, T, LDT, C, LDC, WORK, &
                         LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMQL
!
      END
      SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNMQR overwrites the general complex M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'C':      Q**H * C       C * Q**H
!
!  where Q is a complex unitary matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**H from the Left;
!          = 'R': apply Q or Q**H from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'C':  Conjugate transpose, apply Q**H.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          ZGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQRF.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
                         LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      COMPLEX*16         T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNM2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         NB = MIN( NBMAX, ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N, K, &
              -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMQR', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL ZLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
                         LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H' is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H' is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H'
!
            CALL ZLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                         IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, &
                         WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMQR
!
      END
      SUBROUTINE ZUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNMTR overwrites the general complex M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'C':      Q**H * C       C * Q**H
!
!  where Q is a complex unitary matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by ZHETRD:
!
!  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**H from the Left;
!          = 'R': apply Q or Q**H from the Right.
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from ZHETRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from ZHETRD.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'C':  Conjugate transpose, apply Q**H.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  A       (input) COMPLEX*16 array, dimension
!                               (LDA,M) if SIDE = 'L'
!                               (LDA,N) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by ZHETRD.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!
!  TAU     (input) COMPLEX*16 array, dimension
!                               (M-1) if SIDE = 'L'
!                               (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZHETRD.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >=M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, UPPER
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNMQL, ZUNMQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) &
                THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'ZUNMQL', SIDE // TRANS, M-1, N, M-1, &
                    -1 )
            ELSE
               NB = ILAENV( 1, 'ZUNMQL', SIDE // TRANS, M, N-1, N-1, &
                    -1 )
            END IF
         ELSE
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M-1, N, M-1, &
                    -1 )
            ELSE
               NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N-1, N-1, &
                    -1 )
            END IF
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( LEFT ) THEN
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      END IF
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to ZHETRD with UPLO = 'U'
!
         CALL ZUNMQL( SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C, &
                      LDC, WORK, LWORK, IINFO )
      ELSE
!
!        Q was determined by a call to ZHETRD with UPLO = 'L'
!
         IF( LEFT ) THEN
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         END IF
         CALL ZUNMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, &
                      C( I1, I2 ), LDC, WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMTR
!
      END
