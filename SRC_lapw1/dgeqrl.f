        SUBROUTINE DGEQRL( DECOMP, M, N, A, LDA, DRPTOL, &
                           TAU, WORK )
        IMPLICIT NONE
!
! ----------------------------------------------------------------------
!
! Purpose:
!
!   Compute the QR or QL factorization of the (m X n) matrix A,
!   where m >= n.
!
!   The matrix Q is represented as a product of elementary reflectors
!
!      Q  =  H( 1 ) * ... * H( n )        (QR decomposition)
!
!   or
!
!      Q  =  H( n ) * ... * H( 1 )        (QL decomposition).
!
!   Each H( j ) has the form
!
!      H( j )  =  I  -  tau * v * v'
!
!   where tau is a scalar, and v is a vector of length m.
!
!   In the QR case,
!     v( 1:j-1 ) = 0   and   v( j ) = 1   and, on exit, v( j+1:m ) is
!     stored in a( j+1:m, j ) and tau in tau( j ).
!
!   In the QL case,
!     v( m-n+j+1:m ) = 0   and   v( m-n+j ) = 1   and, on exit,
!     v( 1:m-n+j-1 ) is stored in a( 1:m-n+j-1, j ) and tau in tau( j ).
!
! Date:   96/02/20
! Author: Bruno Lang
!         University of Wuppertal
!         na.blang@na-net.ornl.gov
!
! Note: >>> This routine performs no sanity checks on its arguments. <<<
!
! Parameters:
!
        CHARACTER*1           DECOMP
        INTEGER               M, N, LDA
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), TAU( * ), WORK( * )
!
!   decomp  (input)     Which decomposition is required ?
!                       decomp = 'L' : QL decomposition.
!                              = 'R' : QR decomposition.
!
!   m       (input)     Number of rows of the matrix A.
!                       m >= 0.
!
!   n       (input)     Number of columns of the matrix A.
!                       0 <= n <= m.
!
!   a       (in/out)    On entry, the ( m X n ) matrix A.
!                       If decomp = 'R' then, on exit, the elements on
!                       and above the diagonal contain the ( n X n )
!                       upper triangular matrix R; the elements below
!                       the diagonal, together with the array tau,
!                       represent the orthogonal matrix Q as a product
!                       of Householder transformations.
!                       If decomp = 'L' then, on exit, the (n X n)
!                       matrix L is stored in the lower triangle at
!                       the "bottom" of A, and the Householder vectors
!                       are stored in the remaining trapezoid.
!                       Size : ( lda, n ).
!
!   lda     (input)     Leading dimension of the array a.
!                       lda >= min( m, 1 ).
!
!   tau     (output)    The scaling factors of the Householder trans-
!                       formations.
!                       Size : n.
!
!   drptol  (input)     Threshold for dropping Householder transforms.
!                       If the norm of the vector to be zeroed is
!                       already less than drptol then the transform is
!                       skipped.
!                       drptol >= 0.0.
!
!   work    (workspace) Size : n.
!
! ----------------------------------------------------------------------
!
! Routines called:
!
        EXTERNAL              LSAME, DLBRFG, DLARF
        LOGICAL               LSAME
!
!   dlbrfg  generate Householder vector
!   dlarf   LAPACK
!
        INTRINSIC             MIN
!        INTEGER               MIN
!
! Local variables:
!
        DOUBLE PRECISION      ONE
        PARAMETER             ( ONE = 1.0D0 )
!
        INTEGER               I, J
        DOUBLE PRECISION      TMP
!
! ----------------------------------------------------------------------
!
        IF ( N .GT. 0 ) &
        THEN
!
          IF ( LSAME( DECOMP, 'L' ) ) &
          THEN
!
!                --- QL factorization ---
!
            DO 100 J = N, 1, -1
!
!                  --- generate Householder vector to zero out
!                      A( 1:m-n+j, j)                    ---
!
              I = M - N + J
              CALL DLBRFG( I, A( I, J ), A( 1, J ), 1, TAU( J ), &
                           DRPTOL )
!
!                  --- apply transformation from the left ---
!
              TMP = A( I, J )
              A( I, J ) = ONE
              CALL DLARF( 'Left', I, J - 1, A( 1, J ), 1, TAU( J ), &
                          A, LDA, WORK )
              A( I, J ) = TMP
!
  100       CONTINUE
!
          ELSE
!
!                --- QR factorization ---
!
            DO 200 J = 1, N
!
!                  --- generate Householder vector to zero out
!                      A( j+1:m, j )                           ---
!
!              write(6,*) 'In dgeqrl, vor dlbrfg, tau(j) = ', tau(j)
!	      write(6,*) (a(i,j), i=j+1,min(j+4,m))
              CALL DLBRFG( M - J + 1, A( J, J ), &
                           A( MIN( J + 1, M ), J ), 1, TAU( J ), &
                           DRPTOL )
!
!                  --- apply transformation from the left ---
!
              TMP = A( J, J )
              A( J, J ) = ONE
!              write(6,*) 'In dgeqrl, vor dlarf, tau(j) = ', tau(j)
              CALL DLARF( 'Left', M - J + 1, N - J, A( J, J ), 1, &
                          TAU( J ), A( J, J + 1 ), LDA, WORK )
!              write(6,*) 'In dgeqrl, nach dlarf, tau(j) = ', tau(j)
              A( J, J ) = TMP
!
  200       CONTINUE
!
          ENDIF
!
        ENDIF
!
        RETURN
        END
