      SUBROUTINE DSCGST( N, AB, LDAB, ABDIAG, PANELB, INFO )

        IMPLICIT NONE
!
!     .. Scalar Arguments ..
      INTEGER            INFO, ITYPE, LDAB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, N )
      DOUBLE PRECISION   ABDIAG( N )
      DOUBLE PRECISION   PANELB( N, * )
!     ..
!
!  Purpose
!  =======
!
!  DSCGST reduces a real symmetric-definite generalized eigenproblem
!  to standard form.
!
!  The problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U**T)*A*inv(U)
!
!  B must have been previously factorized as U**T*U by DPOTRF.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          Lower Half (see also ABDIAG): 
!          On entry, the symmetric matrix A.  
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!          Strict Upper half:
!          The triangular factor from the Cholesky factorization of B,
!          as returned by DPOTRF.
!          The diagonal of B is stored in ABDIAG.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB. LDAB >= max(1,N).
!
!  ABDIAG  (input/output) DOUBLE PRECISION array, dimension (N)
!          Diagonal elements of B.
!
!  PANELB  Workspace
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            K, KB, NB
      INTEGER            I, J, LDPANELB

!     ..
!     .. External Subroutines ..
      EXTERNAL           DSYGS2, DSYMM, DSYR2K, DTRMM, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDAB.LT.MAX( 1, N ) ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSCGST', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DSYGST', 'L', N, -1, -1, -1 )
!
!              Compute inv(U')*A*inv(U)      with A lower stored
!
      LDPANELB = N
      DO 50 K = 1, N, NB
         KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(k:n,k:n)
!
         DO I = 1, MIN(KB, N-K)
!         DO I = 1, KB
	   CALL DCOPY( N-K-I+1, AB( K+I-1, K+I ), LDAB,  &
                           PANELB( I+1, I ), 1 )
         END DO
	 CALL DCOPY( KB, ABDIAG( K ), 1, PANELB, LDPANELB+1 )
         CALL DSYGS2( 1, 'L', KB, AB( K, K ), LDAB, &
                               PANELB, LDPANELB, INFO )
         IF( K+KB.LE.N ) THEN
           CALL DTRSM( 'Right', 'L', 'Transpose', 'Non-unit', &
                       N-K-KB+1, KB, ONE, PANELB, LDPANELB, &
                       AB( K+KB, K ), LDAB )
           CALL DSYMM( 'Right', 'L', N-K-KB+1, KB, -HALF, &
                       AB( K, K ), LDAB,  &
                       PANELB( KB+1, 1 ), LDPANELB, ONE, &
                       AB( K+KB, K ), LDAB )
           CALL DSYR2K( 'L', 'No transpose', N-K-KB+1, KB, &
                       -ONE, AB( K+KB, K ), LDAB,  &
                       PANELB( KB+1, 1 ), LDPANELB,  &
                       ONE, AB( K+KB, K+KB ), LDAB )
           CALL DSYMM( 'Right', 'L', N-K-KB+1, KB, -HALF, &
                       AB( K, K ), LDAB,  &
                       PANELB( KB+1, 1 ), LDPANELB, ONE, &
                       AB( K+KB, K ), LDAB )
	   CALL DSWAP( N-K-KB+1, ABDIAG(K+KB), 1,  &
                                 AB(K+KB, K+KB), LDAB+1)
           CALL DTRSM( 'Left', 'U', 'Transpose', &
                       'Non-unit', N-K-KB+1, KB, ONE, &
                       AB( K+KB, K+KB ), LDAB, AB( K+KB, K ), &
                       LDAB )
	   CALL DSWAP( N-K-KB+1, ABDIAG(K+KB), 1,  &
                                 AB(K+KB, K+KB), LDAB+1)
         END IF
   50  CONTINUE
	
      RETURN
!
!     End of DSYGST
!
      END
