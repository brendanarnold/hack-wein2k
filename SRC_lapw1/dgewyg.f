        SUBROUTINE DGEWYG( UPLO, M, N, A, LDA, TAU, &
                           W, LDW, Y, LDY, K, &
                           WORK )
        IMPLICIT NONE
!
! ----------------------------------------------------------------------
!
! Purpose:
!
!   Given a set of Householder vectors as returned by dgeqrl, this
!   routine returns W and Y such that
!
!      Q  =  I  +  W * Y'
!
!   where Q is the product of the Householder matrices.
!
!   Every Householder transformation is represented by a vector vj and
!   a scalar tau( j ) such that
!
!      H( j )  =  I  -  tau( j ) * ( 1, vj' )' * ( 1, vj' )  .
!
!   The matrix Q is accumulated as
!
!      Q( j + 1 )  =  I  +  W( j + 1 ) * Y( j + 1 )'  =  Q( j ) * H( j )
!
!   where
!
!      Y( j + 1 )  =  [ Y( j ), y ]   with  y' = [ 0, vj' ]   and
!
!      W( j + 1 )  =  [ W( j ), w ]   with  w = - tau( j ) * Q( j ) * y.
!
! Date:   96/02/16
! Author: Bruno Lang
!         University of Wuppertal
!         na.blang@na-net.ornl.gov
!
! Note: >>> This routine performs no sanity checks on its arguments. <<<
!
! Parameters:
!
        CHARACTER*1           UPLO
        INTEGER               M, N, LDA, LDW, LDY, K
        DOUBLE PRECISION      A( LDA, * ), TAU( * ), &
                              W( LDW, * ), Y( LDY, * ), WORK( * )
!
!   uplo    (input)     Are the Housholder vectors stored in the upper
!                       or lower triangle of A ?
!                       uplo = 'U' : Upper triangle.
!                            = 'L' : Lower triangle.
!
!   m       (input)     Row dimension of the matrix A.
!                       m >= 0.
!
!   n       (input)     Column dimension of A.
!                       0 <= n <= m.
!
!   a       (input)     The matrix containing the Householder vectors.
!                       Size : ( lda, n ).
!
!   lda     (input)     Leading dimension of the array a.
!                       lda >= max( m, 1 ).
!
!   tau     (input)     Array containing the scaling factors for the
!                       Householder transformations.
!                       Size : n.
!
!   w       (output)    On exit, the first k columns of w contain the
!                       block reflector W.
!                       Size : ( ldw, n ).
!
!   ldw     (input)     Leading dimension of the array w.
!                       ldw >= max( m, 1 ).
!
!   y       (output)    On exit, the first k columns of y contain the
!                       block reflector Y.
!                       Size : ( ldy, n ).
!
!   ldy     (input)     Leading dimension of the array y.
!                       ldy >= max( m, 1 ).
!
!   k       (output)    On exit, k is the "active blocksize"
!                       = the number of valid columns in W and Y
!                       = the number of nontrivial Householder
!                         transforms contained in W and Y
!                       = the number of nonzero tau's in the input.
!                       0 <= k <= n.
!
!   work    (workspace) Size : n - 1.
!
! ----------------------------------------------------------------------
!
! Routines called:
!
        EXTERNAL              LSAME, DCOPY, DGEMV
        LOGICAL               LSAME
!
!   lsame   BLAS
!   dcopy   BLAS
!   dgemv   BLAS
!
        INTRINSIC             MIN
!        INTEGER               MIN
!
! Local variables:
!
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
!
        INTEGER               J, I
        DOUBLE PRECISION      TMP
!
! ----------------------------------------------------------------------
!
        IF ( ( M .GT. 0 ) .AND. ( N .GT. 0 ) ) &
        THEN
!
          K = 0
!
          IF ( LSAME( UPLO, 'Lower' ) ) &
          THEN
!
!                --- vectors come from a QR decomposition ---
!
            DO 110 J = 1, N
!
              TMP = - TAU( J )
!
              IF ( TMP .NE. ZERO ) &
              THEN
                K = K + 1
!
!                    --- Y( :, k ) := vj ---
!
                CALL DCOPY( J - 1, ZERO, 0, Y( 1, K ), 1 )
                Y( J, K ) = ONE
                CALL DCOPY( M - J, A( J + 1, J ), 1, Y( J + 1, K ), 1 )
!
!                    --- W( :, k ) := - tau( j ) * vj ---
!
                CALL DCOPY( J - 1, ZERO, 0, W( 1, K ), 1 )
                DO 100 I = J, M
                  W( I, K ) = TMP * Y( I, K )
  100           CONTINUE
!
!                     --- apply previous transformations ---
!
                IF ( K .GT. 1 ) &
                THEN
                  CALL DGEMV( 'Transpose', M, K - 1, &
                              ONE, Y, LDY, W( 1, K ), 1, &
                              ZERO, WORK, 1 )
                  CALL DGEMV( 'NoTranspose', M, K - 1, &
                              ONE, W, LDW, WORK, 1, &
                              ONE, W( 1, K ), 1 )
                ENDIF
!
              ENDIF
!
  110       CONTINUE
!
          ELSE
!
!                --- vectors come from a QL decomposition ---
!
            DO 210 J = N, 1, -1
!
              TMP = - TAU( J )
!
              IF ( TMP .NE. ZERO ) &
              THEN
                K = K + 1
!
!                    --- Y( :, k ) := vj ---
!
                CALL DCOPY( N - J, ZERO, 0, Y( M - N + J + 1, K ), 1 )
                Y( M - N + J, K ) = ONE
                CALL DCOPY( M - N + J - 1, A( 1, J ), 1, Y( 1, K ), 1 )
!
!                    --- W( :, k ) := - tau( j ) * vj ---
!
                CALL DCOPY( N - J, ZERO, 0, W( M - N + J + 1, K ), 1 )
                DO 200 I = 1, M - N + J
                  W( I, K ) = TMP * Y( I, K )
  200           CONTINUE
!
!                     --- apply previous transformations ---
!
                IF ( K .GT. 1 ) &
                THEN
                  CALL DGEMV( 'Transpose', M, K - 1, &
                              ONE, Y, LDY, W( 1, K ), 1, &
                              ZERO, WORK, 1 )
                  CALL DGEMV( 'NoTranspose', M, K - 1, &
                              ONE, W, LDW, WORK, 1, &
                              ONE, W( 1, K ), 1 )
                ENDIF
!
              ENDIF
!
  210       CONTINUE
!
          ENDIF
!
        ENDIF
!
        RETURN
        END
