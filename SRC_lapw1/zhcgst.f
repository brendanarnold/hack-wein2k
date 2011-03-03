      SUBROUTINE ZHCGST( N, AB, LDAB, ABDIAG, PANELB, INFO )
      IMPLICIT NONE

!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDAB, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16   AB( LDAB, N )
      COMPLEX*16   ABDIAG( N )
      COMPLEX*16   PANELB( N, * )
!     ..
!
!  Purpose
!  =======
!
!  ZHCGST reduces a complex Hermitian-definite generalized eigenproblem
!  to standard form.
!
!  The problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U**H)*A*inv(U)
!
!  B must have been previously factorized as U**H*U by ZPOTRF.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          Lower Half (see also ABDIAG): 
!          On entry, the symmetric matrix A.  
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!          Strict Upper half:
!          The triangular factor from the Cholesky factorization of B**T,
!          as returned by ZPOTRF.
!          The diagonal of B is stored in ABDIAG.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB. LDAB >= max(1,N).
!
!  ABDIAG  (input/output) COMPLEX*16 array, dimension (N)
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
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      COMPLEX*16         CONE, HALF
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ), &
                         HALF = ( 0.5D+0, 0.0D+0 ) )

!     ..
!     .. Local Scalars ..
      INTEGER            K, KB, NB
      INTEGER            I, LDPANELB

!     ..
!     .. External Subroutines ..
      EXTERNAL           ZHEGS2, ZHEMM, ZHER2K, ZTRSM, XERBLA
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
         CALL XERBLA( 'ZHCGST', -INFO )
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
      NB = ILAENV( 1, 'ZHEGST', 'L', N, -1, -1, -1 )
!
!      conjugate overlap matrix
!
!_COMPLEX       DO J = 1, NVAA+NLO
!_COMPLEX         DO I = 1, J
!_COMPLEX          HS(I,J) = DCONJG(HS(I,J))
!_COMPLEX        END DO
!_COMPLEX       END DO
!
!              Compute inv(L)*A*inv(L')      with A lower stored
!
      LDPANELB = N
      DO 50 K = 1, N, NB
         KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(k:n,k:n)
!
!        Copy matrix panel of B into PANELB
         DO I = 1, MIN(KB, N-K)
!         DO I = 1, KB
	   CALL ZCOPY( N-K-I+1, AB( K+I-1, K+I ), LDAB,  &
                           PANELB( I+1, I ), 1 )
         END DO
	 CALL ZCOPY( KB, ABDIAG( K ), 1, PANELB, LDPANELB+1 )

!        all accesses to the array AB are accesses to matrix A
!        accesses to matrix B are via array PANELB
         CALL ZHEGS2( 1, 'L', KB, AB( K, K ), LDAB, &
                               PANELB, LDPANELB, INFO )
         IF( K+KB.LE.N ) THEN
           CALL ZTRSM( 'Right', 'L', 'Conjugate Transpose', 'Non-unit', &
                       N-K-KB+1, KB, CONE, PANELB, LDPANELB, &
                       AB( K+KB, K ), LDAB )
           CALL ZHEMM( 'Right', 'L', N-K-KB+1, KB, -HALF, &
                       AB( K, K ), LDAB,  &
                       PANELB( KB+1, 1 ), LDPANELB, CONE, &
                       AB( K+KB, K ), LDAB )
           CALL ZHER2K( 'L', 'No transpose', N-K-KB+1, KB, &
                       -CONE, AB( K+KB, K ), LDAB,  &
                       PANELB( KB+1, 1 ), LDPANELB,  &
                       ONE, AB( K+KB, K+KB ), LDAB )
           CALL ZHEMM( 'Right', 'L', N-K-KB+1, KB, -HALF, &
                       AB( K, K ), LDAB,  &
                       PANELB( KB+1, 1 ), LDPANELB, CONE, &
                       AB( K+KB, K ), LDAB )

!          Access to matrix B => Swap diagonal elements into array AB
	   CALL ZSWAP( N-K-KB+1, ABDIAG(K+KB), 1,  &
                                 AB(K+KB, K+KB), LDAB+1)
!                                   Not Conjugate Transpose
!                      matrix B                matrix A
           CALL ZTRSM( 'Left', 'U', 'Transpose', &
                       'Non-unit', N-K-KB+1, KB, CONE, &
                       AB( K+KB, K+KB ), LDAB, AB( K+KB, K ), &
                       LDAB )
!          Swap diagonal elements of matrix B out of array AB
	   CALL ZSWAP( N-K-KB+1, ABDIAG(K+KB), 1,  &
                                 AB(K+KB, K+KB), LDAB+1)
         END IF
   50  CONTINUE
!
!      conjugate overlap matrix
!
!_COMPLEX       DO J = 1, NVAA+NLO
!_COMPLEX         DO I = 1, J
!_COMPLEX          HS(I,J) = DCONJG(HS(I,J))
!_COMPLEX        END DO
!_COMPLEX       END DO

	
      RETURN
!
!     End of ZHEGST
!
      END
