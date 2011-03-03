      SUBROUTINE DSBEIN1( N, A, LDA, ABAND, M, W, Z, LDZ, RDB, &
                          WORK, LWORK, IFAIL, INFO, IBLOCK, ISPLIT )
!
!
!                                                                    1997
!  Based on LAPACK routines                              Dieter Kvasnicka
!                                                   Theoretical Chemistry
!                                         University of Technology Vienna
!
      implicit none
!     .. Scalar Arguments ..
      INTEGER            N, LDA, M, LDZ, RDB, LWORK, INFO
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( M )
      INTEGER            IBLOCK( N ), ISPLIT( N )
!       ABAND: Leading Dimension LDZ: used by BANDV
!              second Dimension: 2*(rdb+1)
      DOUBLE PRECISION   A( LDA, N ), ABAND( LDZ, 2*(RDB+1) ), W( N ),  &
                         WORK( LWORK ), Z( LDZ, M )
!     ..
!
!  Purpose
!  =======
!
!  DSBEIN1 computes the eigenvectors of a real symmetric band
!  matrix corresponding to specified eigenvalues, using inverse
!  iteration.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA, N)
!          The band matrix is stored in the diagonal and 
!          subdiagonals of the matrix A
!
!  LDA     (input) INTEGER
!          Leading dimension of A
!
!  ABAND   (working array) DOUBLE PRECISION array, dimension (LDZ, 2*(RDB+1))
!          Local copy of the band matrix
!
!  M       (input) INTEGER
!          The number of eigenvectors to be found.  0 <= M <= N.
!
!  W       (input) DOUBLE PRECISION array, dimension (N)
!          The first M elements of W contain the eigenvalues for
!          which eigenvectors are to be computed.  
!          (The eigenvalues                                          )
!          (should be grouped by split-off block and ordered from    )
!          (smallest to largest within the block.  ( The output array)
!          (W from DSTEBZ with ORDER = 'B' is expected here. )       )
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)
!          The computed eigenvectors.  The eigenvector associated
!          with the eigenvalue W(i) is stored in the i-th column of
!          Z.  Any vector which fails to converge is set to its current
!          iterate after a number of iterations.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
!
!  LWORK   (input) INTEGER
!          Size of the working array WORK
!
!  IFAIL   (output) INTEGER array, dimension (M)
!          On normal exit, all elements of IFAIL are zero.
!          If one or more eigenvectors fail to converge after
!          a number of iterations, then their indices are stored in
!          array IFAIL.
!          Used only if the number of subdiagonal of the band matrix
!          is one.
!
!  INFO    (output) INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, then i eigenvectors failed to converge
!               in a number of iterations.  Their indices are stored in
!               array IFAIL.
!
!  IBLOCK  (input) INTEGER array, dimension (N)
!          The submatrix indices associated with the corresponding
!          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
!          the first submatrix from the top, =2 if W(i) belongs to
!          the second submatrix, etc.  ( The output array IBLOCK
!          from DSTEBZ is expected here. )
!          Used only if the number of subdiagonal of the band matrix
!          is one.
!
!  ISPLIT  (input) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to
!          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!          through ISPLIT( 2 ), etc.
!          ( The output array ISPLIT from DSTEBZ is expected here. )
!          Used only if the number of subdiagonal of the band matrix
!          is one.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          (ONE = 1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, IERR
      DOUBLE PRECISION   tarray(3)
!     ..
!     .. Local Arrays ..
      INTEGER		 IWORK( 5000 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE

      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSBEIN1', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     Initialize band matrix
!
!      DO 15 I = 1, LDZ
!	DO 13 J = 1, (RDB+1)*2
!	  ABAND(I,J) = 0.0
!  13    CONTINUE
!  15  CONTINUE
!
!     Copy band matrix 
!
      call dtime4(tarray)
      DO 30 J = 1, RDB+1
	 DO 20 I = RDB-J+2, N
	    ABAND( I, J ) = A( I, I+J-RDB-1 )
!	    ABAND( I, J ) = A( I+J-RDB-1, I )
  20     CONTINUE
  30  CONTINUE

      call printtime4(' Copy matrix ', tarray)
!
!     Compute Eigenvectors using inverse iteration by calling 
!     the EISPACK routine BANDV
!
      IF(rdb.gt.1) THEN
!	write(6,*) ' Eigenvalues:'
!	write(6,*) (W(I), I=1,M)
        CALL BANDV1( LDZ, N, RDB+1, ABAND, 0.0d0, M, W, Z, IERR, &
                    LWORK - N, WORK( N+1 ), WORK )
        IF( IERR.NE.0 ) THEN
	  write(6,*) 'IERR after BANDV: ', IERR
	  STOP 'ERROR IN INVERSE ITERATION ROUTINE BANDV'
	END IF
        call printtime4(' BANDV ', tarray)

!
!     Use LAPACK routines for a tridiagonal matrix
!
      ELSE 
!	write(6,*) 'D: ', (ABAND(I,2), I=1,N)
!	write(6,*) 'E: ', (ABAND(I+1,1), I=1,N)
	CALL DSTEIN( N, ABAND(1,2), ABAND(2,1), M, W, &
                     IBLOCK, ISPLIT, Z, LDZ, &
                     WORK, IWORK, IFAIL, IINFO)
	IF( IINFO.NE.0 ) write(6,*) 'INFO after DSTEIN: ', IINFO
        call printtime4(' DSTEIN ', tarray)

      END IF

!
!


      RETURN
!
!     End of DSBEIN1
!
      END
