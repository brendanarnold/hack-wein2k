! **********************************************************************
!
        SUBROUTINE DSYRB4( UPLO, N, RDB, A, LDA,  &
                           NB, TAU, WORK, LWORK, INFO )
!
! ----------------------------------------------------------------------
!
! Purpose:
!
!   Given a symmetric n-times-n matrix A with either its lower or upper
!   triangle explicitly stored, DSYRB4 reduces A to a symmetric banded
!   matrix B with RDB sub(super)diagonals by using a sequence of
!   orthogonal similarity transformations. 
!
!   A must be given in full storage (with either its lower or upper
!   triangle explicitly stored) and will be overwritten with B and
!   the "mantissae" of the Householder vectors used for the reduction.
!
!   Works only for lower stored matrices.
!
!
! Parameters:
!
        implicit none
        CHARACTER*1           UPLO 
        INTEGER               N, RDB, LDA, NB, LWORK, INFO
        DOUBLE PRECISION      A( LDA, N ), TAU( N ), WORK( LWORK )
!
!   uplo    (input)     Reduce the upper or lower triangle of A ?
!                       uplo = 'U' : Upper triangle.
!                            = 'L' : Lower triangle.
!
!   n       (input)     The order of the matrix A.
!                       n >= 0.
!
!   RDB       (input)     The number of sub(super-)diagonals of the
!                       reduced matrix B.
!                       1 <= RDB < n, if n > 1.
!                       RDB = 0     , if n <= 1.
!
!   a       (in/out)    On entry, this array contains the symmetric
!                       matrix A with either its lower (uplo = 'L') or
!                       upper (uplo = 'U') triangle stored.
!                       On exit, the main diagonal and the first RDB
!                       sub(super-)diagonals contain the banded
!                       matrix B, and the Householder vectors that were
!                       used in the reduction are stored in the
!                       remaining diagonals of the lower (upper)
!                       triangle.
!                       Size : ( lda, n ).
!
!   lda     (input)     Leading dimension of the array a.
!                       lda >= max( n, 1 ).
!
!   nb      (input)     Blocking factor, as suggested by the user, for
!                       the reduction.
!                       nb = RDB : Do no two-level blocking.
!                       nb > RDB : Try a rank-2*nb update
!                                (if the workspace is not sufficient
!                                then a smaller blocking factor will be
!                                used).
!
!   tau     (output)    Scaling factors for the Householder transforms
!                       used in the reduction.
!                       Size : n - RDB - 1.
!
!   work    (workspace) Size : lwork.
!
!   lwork   (input)     Size of the workspace buffer (must be provided
!                       on entry).
!                       lwork >= max(nb,RDB)*4*n.
!
!   info    (output)    On exit, info indicates the consistence of the
!                       arguments:
!
!                       info = - 1 : uplo is neither 'U' nor 'L' (upper
!                                    or lower case).
!                            = - 2 : n is out of range (negative).
!                            = - 3 : RDB is out of range ( < 0, or = 0
!                                    while n > 1, or >= n ).
!                            = - 5 : lda is out of range ( < n or < 1).
!                            = - 6 : nb is out of range (negative).
!                            = - 9 : lwork is too small (see above).
!
!                       info >=  1 : All arguments are OK.
!                                    In this case, info returns the
!                                    blocking factor that was used in
!                                    the reduction and in U's update.
!
! ----------------------------------------------------------------------
!
! Routines called:
!
        EXTERNAL              LSAME, DGEQRL, DGEWYG, DGEWY
!        EXTERNAL              NBDFLT
!       external              DSYWY
        LOGICAL               LSAME
!        INTEGER               NBDFLT
!
!   lsame   BLAS
!   nbdflt  determine default blocking factor
!   dgeqrl  QR or QL factorization
!   dgewyg  generate WY factors
!   dgewy   apply WY transform from one side
!   dsywy   apply two-sided WY transform to a symmetric matrix
!
        INTRINSIC             MAX, MIN
!        INTEGER               MAX, MIN
!
! Local variables:
!
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0,  &
                                ONE = 1.0D0 )
!
        LOGICAL               LOWER 
        INTEGER               NB0, IW, Neffective, BO, BI, NO,  &
                              K, SO, RO, LENO, WIDTHO,  &
                              IWO, IYO, IVO, IWORK
!
!
! ----------------------------------------------------------------------
!
! The workspace is used as follows:
!
! - The first nb * n <= nb0 * ( n - RDB ) elements hold the block W of
!   the current block WY transform.
!   Here, nb <= nb0 denotes the width of the current block column and
!   n <= n - RDB denotes its height; nb0 is the final blocking factor.
! - The next nb * n elements hold the block Y.
! - The next nb * n elements hold the block V.
! - Another ( n + nb ) * nb <= n * nb0 elements are used as
!   workspace in the routines dgeqrl, dgewyg, dsywyv, dgewy, and dsywy.
!
! ----------------------------------------------------------------------
!
!            --- check for errors in input ---
!
        LOWER = LSAME( UPLO, 'L' )
!
        IF ( .NOT. LOWER ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( ( N .LE. 1 ) .AND. ( RDB .NE. 0 ) ) .OR. &
                 ( ( N .GT. 1 ) .AND. &
                   ( ( RDB .LT. 1 ) .OR. ( RDB .GE. N ) ) ) ) &
        THEN
          INFO = - 3
        ELSEIF ( LDA .LT. MAX( N, 1 ) ) THEN
          INFO = - 5
        ELSEIF ( NB .LT. 0 ) THEN
          INFO = -6
        ELSEIF ( LWORK .LT. ( 4*MAX(NB,RDB)*N ) ) THEN
          INFO = -9
        ELSE
          INFO = 1
        ENDIF
!
        IF ( INFO .NE. 1 )     GOTO 999
!
!            --- check for quick return ---
!
        IF ( ( N .LE. 1 ) .OR. ( RDB .EQ. ( N - 1 ) ) )      GOTO 999
!
!            --- compute the maximum blocking factor for the
!                given workspace                             ---
!
        NB0 = LWORK / ( 4*N )
!
!            --- adjust to the user-supplied or default value ---
!
        IF ( NB .GT. 0 ) &
        THEN
          NB0 = MIN( NB, NB0 )
        ENDIF


!
! ......................................................................
!
        IW = 1
!
        IF ( LOWER ) THEN
!
!              --- reduce the lower triangle ---
!

!         N is the order of the matrix A

!         Neffective is the number of columns of matrix A that have to be reduced
!	  Neffective = N - RDB
	  Neffective = N - RDB - 1

!         BO is the blocking of the outer loop
          BO = NB0

!         BI is the blocking of the inner loop
          BI = RDB 

!         NO is the number of outer blocks ( =ceiling(Neffective/BO) ) 
	  NO = (Neffective-1)/BO+1

!         IWO is the outer block's W
          IWO = 1
!         IYO is the outer block's Y
          IYO = IWO + BO * N
!         IVO is the outer block's V
          IVO = IYO + BO * N
!         IWORK is a work array for subprograms
          IWORK = IVO + BO * N

          DO K = 1, NO

!           SO is the first column of the current outer block
	    SO = (K-1) * BO + 1
!           RO is the last column of the current outer block
            RO = MIN( K * BO, Neffective )

!           Length of the current outer block
	    LENO = N - SO - BI + 1
!           Width of the current outer block
            WIDTHO = RO - SO + 1

            CALL DSYRB5L( BI, LENO, WIDTHO, WORK( IWO ), &
                   WORK( IYO ), WORK( IVO ),  &
                   WORK( IWORK), LWORK - IWORK + 1, &
                   A( SO, SO ), LDA, N - SO + 1, TAU( SO ) )

          END DO
!
        ENDIF
  999   RETURN
        END
