        SUBROUTINE DSYWYV( UPLO, N, K, KI, A, LDA,  &
                          WI, LDWI, YI, LDYI,  &
                          YO, LDYO, VO, LDVO,  &
                          VI, LDVI,  &
                          WORK, LWORK )
!
! ----------------------------------------------------------------------
!
! Purpose:
!
!   Generate a matrix VI to find a rank-2k update representation
!   for the two-sided WY block transform
!   i.e., replace the transform ( I + WI * YI' )' * A * ( I + WI * YI' )
!   by a representation A + VI YI' + YI VI'.
!   If (K.GE.0) an additional update has to be performed using VO and YO
!   before computing V.
!
!   The reason for introducing the VO and YO terms is that the matrix A
!   used in computing VI
!   must be replaced by (A + VO YO' + YO VO'), because the rank-2k updates using
!   VO and YO have been delayed.
!   So the computation which is really performed is 
!     ( I + WI YI' )' * (A + VO YO' + YO VO') * ( I + WI YI' )

!   The method used is
!   VI = A WI
!   VI = VI + YO VO' WI + VO YO' WI
!   VI = VI + 1/2 YI VI' WI
! 
!   i.e. VI = A WI + YO VO' WI + VO YO' WI +
!             + 1/2 YI (A WI + YO VO' WI + VO YO' WI)' WI 
!   thus
!   A + VI YI' + YI VI' =
!   = A + (A WI + YO VO' WI + VO YO' WI + 1/2 YI (A WI + YO VO' WI + VO YO' WI)' WI) YI' +
!       + YI (A WI + YO VO' WI + VO YO' WI + 1/2 YI (A WI + YO VO' WI + VO YO' WI)' WI)' =
!   = A + A WI YI ' + YO VO' WI YI' + VO YO' WI YI' + 
!         + 1/2 YI WI' A' WI YI' + 1/2 YI WI' VO YO' WI YI' + 1/2 YI WI' YO VO' WI YI' + 
!       + YI WI' A' + YI WI' VO YO' + YI WI' YO VO' + 
!         + 1/2 YI WI' A WI YI' + 1/2 YI WI' YO VO' WI YI' + 1/2 YI WI' VO YO' WI YI' =
!   = (since A = A')
!     A + A WI YI' + YI WI' A + YI WI' A WI YI' + 
!       + YO VO' WI YI' + VO YO' WI YI' + YI WI' VO YO' + YI WI' YO VO' +
!       + YI WI' VO YO' WI YI' + YI WI' YO VO' WI YI'
!   This is equal to ( I + WI YI' )' * (A + VO YO' + YO VO') * ( I + WI YI' )
!
!   neglecting the VO and YO terms, this yields
!     A + A WI YI' + YI WI' A + YI WI' A WI YI' , which is equal to 
!     ( I + WI YI' )' * A * ( I + WI YI' )
!
!
! Note: >>> This routine performs no sanity checks on its arguments. <<<
!
! Parameters:
!
        implicit none
        CHARACTER*1           UPLO
        INTEGER               N, K, KI, LDA, LDWI, LDYI, LDVI, LDYO,  &
                              LDVO, LWORK
        DOUBLE PRECISION      A( LDA, N ), WI( LDWI, KI ),  &
                              YI( LDYI, KI ), &
                              YO( LDYO, K ), VO( LDVO, K ), &
                              VI( LDVI, KI ), &
                              WORK( LWORK )
!
!   uplo    (input)     Is the upper or lower triangle of A stored ?
!                       uplo = 'U' : Upper triangle.
!                            = 'L' : Lower triangle.
!
!   n       (input)     Size of the block to transform.
!                       n >= 0.
!
!   k       (input)     Number of columns of YO and VO.
!                       k >= 0.
!
!   ki      (input)     Number of nontrivial Householder transforms in
!                       the block WY transform.
!                       ki >= 0.
!
!   a       (input)     The ( n X n ) symmetric matrix A.
!                       Size : ( lda, n ).
!
!   lda     (input)     Leading dimension of the array a.
!                       lda >= max( n, 1 ).
!
!   wi      (input)     Contains the block reflector WI.
!                       Size : ( ldwi, ki ).
!
!   ldwi    (input)     Leading dimension of the array WI.
!                       ldwi >= max( n, 1 ).
!
!   yi      (input)     Contains the block reflector YI.
!                       Size : ( ldyi, ki ).
!
!   ldyi    (input)     Leading dimension of the array YI.
!                       ldyi >= max( n, 1 ).
!
!   vi      (output)    Contains the block reflector VI.
!                       Size : ( ldvi, ki ).
!
!   ldvi    (input)     Leading dimension of the array VI.
!                       ldvi >= max( n, 1 ).
!
!   yo      (input)     Contains the block reflector YO.
!                       Size : ( ldyo, k ).
!
!   ldyo    (input)     Leading dimension of the array YO.
!                       ldyo >= max( n, 1 ).
!
!   vo      (output)    Contains the block reflector VO.
!                       Size : ( ldvo, k ).
!
!   ldvo    (input)     Leading dimension of the array VO.
!                       ldvo >= max( n, 1 ).
!
!   work    (workspace) Size : lwork >= max(K,KI)*KI.
!
!   lwork   (input)     lwork >= max(K,KI)*KI.
!
! ----------------------------------------------------------------------
!
! Routines called:
!
        EXTERNAL              DSYMM, DGEMM
!
!   dsymm   BLAS
!   dgemm   BLAS
!
! Local variables:
!
        DOUBLE PRECISION      ZERO, ONE, HALF
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
	PARAMETER             ( HALF = 0.5D0 )
!
! ----------------------------------------------------------------------
!
        IF ( ( N .GT. 0 ) .AND. ( KI .GT. 0 ) ) THEN
!         --- VI = A * WI 
          CALL DSYMM( 'Left', UPLO, N, KI, &
                      ONE, A, LDA, WI, LDWI, ZERO, VI, LDVI )
!         --- VI = VI + YO * VO' * WI + VO * YO' * WI
          IF( (K .GT. 0) ) THEN
!           --- WORK = VO' * WI
            CALL DGEMM( 'T', 'N', K, KI, N,  &
                        ONE, VO, LDVO, WI, LDWI, &
                        ZERO, WORK, K )
!           --- VI = VI + YO * WORK
            CALL DGEMM( 'N', 'N', N, KI, K, &
                        ONE, YO, LDYO, WORK, K, &
                        ONE, VI, LDVI )
!           --- WORK = YO' * WI
            CALL DGEMM( 'T', 'N', K, KI, N, &
                        ONE, YO, LDYO, WI, LDWI, &
                        ZERO, WORK, K )
!           --- VI = VI + VO * WORK
            CALL DGEMM( 'N', 'N', N, KI, K,  &
                        ONE, VO, LDVO, WORK, K, &
                        ONE, VI, LDVI )
          END IF
!         --- WORK = 1/2 * VI' * WI 
          CALL DGEMM( 'Transpose', 'NoTranspose', KI, KI, N, &
                      HALF, VI, LDVI, WI, LDWI, &
                      ZERO, WORK, KI )
!         --- VI = VI + YI * WORK 
          CALL DGEMM( 'NoTranspose', 'NoTranspose', N, KI, KI, &
                      ONE, YI, LDYI, WORK, KI, &
                      ONE, VI, LDVI )
!
        ENDIF
!
        RETURN
        END
