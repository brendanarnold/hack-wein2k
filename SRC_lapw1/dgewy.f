        SUBROUTINE DGEWY( SIDE, M, N, K, A, LDA, &
                          W, LDW, Y, LDY, &
                          WORK )
        IMPLICIT NONE
!
! ----------------------------------------------------------------------
!
! Purpose:
!
!   Apply a WY block transform to the matrix A (from the left or from
!   the right), i.e., update A as
!
!      A  :=  ( I + W * Y' )' * A   , or
!
!      A  :=  A * ( I + W * Y' )    , resp.
!
! Date:   95/10/13
! Author: Bruno Lang
!         University of Wuppertal
!         na.blang@na-net.ornl.gov
!
! Note: >>> This routine performs no sanity checks on its arguments. <<<
!
! Parameters:
!
        CHARACTER*1           SIDE
        INTEGER               M, N, K, LDA, LDW, LDY
        DOUBLE PRECISION      A( LDA, N ), W( LDW, K ), Y( LDY, K )
        DOUBLE PRECISION      WORK( * )
!
!   side    (input)     Apply the transforms from the left or from
!                       the right ?
!                       side = 'L' : From the left.
!                            = 'R' : From the right.
!
!   m       (input)     Number of rows of the matrix A.
!                       m >= 0.
!
!   n       (input)     Number of columns of the matrix A.
!                       n >= 0.
!
!   k       (input)     Number of columns of the matrices W and Y.
!                       k >= 0.
!
!   a       (in/out)    On entry, the matrix A.
!                       On exit, ( I + W * Y' )' * A, if side = 'L', or
!                                A * ( I + W * Y' ),  if side = 'R'.
!                       Size : ( lda, n ).
!
!   lda     (input)     Leading dimension of the array a.
!                       lda >= max( m, 1 ).
!
!   w       (input)     The matrix W.
!                       Size : ( ldw, k ).
!
!   ldw     (input)     Leading dimension of the array w.
!                       ldw >= max( m, 1 ), if side = 'L'.
!                           >= max( n, 1 ), if side = 'R'.
!                       (W has row dimension m or n, resp.)
!
!   y       (input)     The matrix Y.
!                       Size : ( ldy, k ).
!
!   ldy     (input)     Leading dimension of the array y.
!                       ldy >= max( m, 1 ), if side = 'L'.
!                           >= max( n, 1 ), if side = 'R'.
!                       (Y has row dimension m or n, resp.)
!
!   work    (workspace) Size : n * k, if side = 'L'.
!                              m * k, if side = 'R'.
!
! ----------------------------------------------------------------------
!
! Routines called:
!
        EXTERNAL              LSAME, DGEMM
        LOGICAL               LSAME
!
!   lsame   BLAS
!   dgemm   BLAS
!
! Local variables:
!
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
!
! ----------------------------------------------------------------------
!
        IF ( ( K .GT. 0 ) .AND. ( M .GT. 0 ) .AND. ( N .GT. 0 ) ) &
        THEN
!
          IF ( LSAME( SIDE, 'L' ) ) &
          THEN
!
!                --- work = W' * A ---
!
            CALL DGEMM( 'Transpose', 'NoTranspose', K, N, M, &
                        ONE, W, LDW, A, LDA, ZERO, WORK, K )
!
!                --- A = Y * work + A ---
!
            CALL DGEMM( 'NoTranspose', 'NoTranspose', M, N, K, &
                        ONE, Y, LDY, WORK, K, ONE, A, LDA )
!
          ELSE
!
!                --- work = A * W ---
!
            CALL DGEMM( 'NoTranspose', 'NoTranspose', M, K, N, &
                        ONE, A, LDA, W, LDW, ZERO, WORK, M )
!
!                --- A = A + work * Y' ---
!
            CALL DGEMM( 'NoTranspose', 'Transpose', M, N, K, &
                        ONE, WORK, M, Y, LDY, ONE, A, LDA )
!
          ENDIF
!
        ENDIF
!
        RETURN
        END
