      SUBROUTINE DSYRDT4( UPLO, N, A, LDA, D, E,  &
                         TAU, ABAND, LDABAND, WORK, LWORK, HB,  &
                         RDB, DRPTOL, INFO)
        IMPLICIT NONE
!
!
!
!  Block packed algorithms                                           1997
!  Based on LAPACK routines                              Dieter Kvasnicka
!                                                   Theoretical Chemistry
!                                         University of Technology Vienna
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N, LDA, LDABAND, LWORK, HB, RDB, INFO
      DOUBLE PRECISION   DRPTOL
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, N ), D( N ), E( N ), TAU( N ), &
                         ABAND( LDABAND, N ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DSYRDT4 reduces a real symmetric matrix A to real symmetric
!  tridiagonal form T by an orthogonal similarity transformation:
!  Q**T * A * Q = T.
!
!  Note: this routine has three parameters which do not exist in the 
!        call to DSPTRD:
!        WORK, LWORK, HB
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'L':  Lower triangle of A is stored;
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors.
!
!  LDA     (input) INTEGER
!          Leading dimension of A
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors 
!
!  ABAND   (workspace) DOUBLE PRECISION array, dimension (LDABAND, N)
!
!  LDABAND (input) INTEGER
!          Leading dimension of ABAND
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
!
!  LWORK   (input) INTEGER
!          LWORK >= N*max(HB,B)*4.
!
!  HB      (input) INTEGER
!          The blocksize used in the outer blocking.
!
!  RDB     (input) INTEGER
!          The blocksize used in the inner blocking.
!
!  DRPTOL  (input) DOUBLE PRECISION
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
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            NB, IINFO, I, J
      DOUBLE PRECISION   tarray(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
!      real               dtime2
!      integer            IJ2K
      EXTERNAL           LSAME, ILAENV
!      EXTERNAL           dtime2, IJ2K
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( UPPER ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LWORK.LT.(N*MAX(HB,RDB)*4) ) THEN
         INFO = -8
      ELSE IF( HB.LT.1 ) THEN
	 INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYRDT4', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 3
         RETURN
      END IF
!
!     Determine the block size.
!
      NB=MIN(HB,N)

      IF( RDB .EQ. 1 ) THEN
        call dtime4(tarray)
	CALL DSYTRD2( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
	call printtime4(' DSYTRD ', tarray)
!c	write(6,*) 'D: ', (D(I), I=1,N)
!	write(6,*) 'E: ', (E(I), I=1,N)
       ELSE
!
!     Reduce bandwidth of matrix to RDB
!
!             Works with lower stored matrices only 
!
          call dtime4(tarray)
          CALL DSYRB4( UPLO, N, RDB, A, LDA,  &
                       NB, TAU, WORK, LWORK, IINFO )
	  IF( IINFO.LT.0 ) write(6,*) 'INFO after DSYRB4: ', IINFO
	  call printtime4(' DSYRB4 ', tarray)
!
!     Reduce band matrix to tridiagonal form 
!       (works with lower stored matrices only)
!
!     Copy Band Matrix
!
        IF( UPLO .EQ. 'U' ) THEN
	  DO 220 J=1,N
	    DO 210 I=MAX(J-RDB,1),J
	      ABAND(RDB+1+I-J,J)=A(I,J)
 210        CONTINUE
 220      CONTINUE
	ELSE
	  DO 270 J=1,N
	    DO 260 I=J,MIN( N, J+RDB )
	      ABAND( 1+I-J , J ) = A( I, J )
 260        CONTINUE
 270      CONTINUE
	END IF

	CALL DSBTRD( 'N', UPLO, N, RDB, ABAND, LDABAND, D, E, &
                       0, N, WORK, IINFO )  
        call printtime4(' DSBTRD ', tarray)

      END IF

      RETURN
!
!     End of DSYRDT4
!
      END
