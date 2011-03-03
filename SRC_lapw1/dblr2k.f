        SUBROUTINE DBLR2K( LENI, WIDTHI, KO,  &
                           V, LDV, Y, LDY, A, LDA )
	IMPLICIT NONE
!
!
!     .. Scalar Arguments ..
	INTEGER            LENI, WIDTHI, KO,  &
                           LDV, LDY, LDA
!     
!     .. Array Arguments ..
	DOUBLE PRECISION   V( LDV, KO ), Y( LDY, KO ),  &
                           A( LDA, WIDTHI )
!
!
!  Purpose
!  =======
!
!  DBLR2K computes a rank-2k update of the trapezoidal matrix A
!
!  A = A + V Y' + Y V'
!
!  This is done by a symmetric rank-2k update of the upper 
!  (symmetric) part of A, followed by two matrix multiplications
!  to compute the unsymmetric part of A.
!
!
!  Arguments
!  =========
!
!  LENI      (input) INTEGER
!            First dimension of the matrix stored in the array A
!
!  WIDTHI    (input) INTEGER
!            Second dimenstion of the matrix stored in the array A
!
!  KO        (input) INTEGER
!            Inner dimension of the rank-2k update (number of 
!            columns of the arrays V and Y)
!
!  V         (input) DOUBLE PRECISION array, dimension (LDV, KO)
!            Matrix V
!
!  LDV       (input) INTEGER
!            Leading dimension of V
!
!  Y         (input) DOUBLE PRECISION array, dimension (LDY, KO)
!            Matrix Y
!
!  LDY       (input) INTEGER
!            Leading dimension of Y
!
!  A         (input/output) DOUBLE PRECISION array, 
!                           dimension (LDA, WIDTHI)
!
!  LDA       (input) INTEGER
!            Leading dimension of A
!
! =====================================================================
!
!     .. Parameters ..
	DOUBLE PRECISION   ONE
	PARAMETER          ( ONE = 1.0D0 )

!
!     .. Executable Statements ..
!
	IF( KO .GT. 0 ) THEN
!         -- A = A + V Y' + Y V', symmetric part
          CALL DSYR2K( 'L', 'N', WIDTHI, KO, &
                     ONE, V, LDV,  &
                          Y, LDY, &
                     ONE, A, LDA )

!         -- A = A + V Y' + Y V', non-symmetric part
!           -- A = A + V Y'
          IF( ( LENI - WIDTHI ) .GT. 0 ) THEN
	    CALL DGEMM( 'N', 'T', LENI - WIDTHI, WIDTHI, KO,  &
                       ONE, V( 1 + WIDTHI, 1 ), LDV, Y, LDY, &
                       ONE, A( 1 + WIDTHI, 1 ), LDA)
!           -- A = A + Y V'
	    CALL DGEMM( 'N', 'T', LENI - WIDTHI, WIDTHI, KO,  &
                       ONE, Y( 1 + WIDTHI, 1 ), LDY, V, LDV, &
                       ONE, A( 1 + WIDTHI, 1 ), LDA)
          END IF
	END IF

	RETURN
	END
