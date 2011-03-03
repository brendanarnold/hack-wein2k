        SUBROUTINE DLBRFG( N, ALPHA, X, INCX, TAU, DRPTOL )
        IMPLICIT NONE
!
! ----------------------------------------------------------------------
!
! Purpose:
!
!   dlbrfg generates an elementary reflector H of order n, such that
!
!      H * ( alpha ) = ( beta ) ,   H' * H = I ,
!          (   x   )   (   0  )
!
!   where alpha and beta are scalars, and x is an ( n - 1 )-element
!   vector. H is represented in the form
!
!      H = I - tau * ( 1 ) * ( 1 v' ) ,
!                    ( v )
!
!   where tau is a scalar and v is an ( n - 1 )-element vector.
!
!   If the elements of x are all zero, or the 2-norm of x doesn't exceed
!   drptol, then x is simply zeroed out, tau is set to 0 indicating
!   that H is the identity matrix. In this case, the transformation is
!   effectively skipped.
!
!   Otherwise  1 <= tau <= 2.
!
!   This routine is a modified version of dlbrfg by Bischof/Sun which
!   in turn is a modified version of the LAPACK routine dlarf.
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
        INTEGER               N, INCX
        DOUBLE PRECISION      ALPHA, TAU, DRPTOL
        DOUBLE PRECISION      X( * )
!
!   n       (input)     The order of the elementary reflector.
!                       n >= 0.
!
!   alpha   (in/out)    On entry, the value alpha.
!                       On exit, it is overwritten with the value beta.
!
!   x       (in/out)    On entry, the vector x.
!                       On exit, the "x elements" are overwritten with
!                       the vector v.
!                       Size : 1 + ( n - 2 ) * incx.
!
!   incx    (input)     The increment between elements of x.
!                       incx > 0.
!
!   tau     (output)    The scaling factor tau.
!
!   drptol  (input)     Threshold for dropping Householder transforms.
!                       If the norm of the vector to eliminate is
!                       already <= drptol then the transform is skipped.
!                       drptol >= 0.0.
!                       If you don't know, use 0.0.
!
! ----------------------------------------------------------------------
!
! Routines called:
!
        EXTERNAL              DNRM2, DCOPY, DLAPY2, DLAMCH, DSCAL
        DOUBLE PRECISION      DNRM2, DLAPY2, DLAMCH
!
!   dlapy2        LAPACK
!   dlamch        LAPACK
!   dnrm2         BLAS
!   dcopy         BLAS
!   dscal         BLAS
!
        INTRINSIC             SIGN, ABS
!        DOUBLE PRECISION      SIGN, ABS
!
! Local variables:
!
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
!
        DOUBLE PRECISION      XNORM, BETA, SAFMIN, RSAFMN
        INTEGER               CNT, J
!
!   xnorm   2-norm of x
!   beta    the computed value beta
!   safmin  safe minimum (such that relative accuracy can be guaranteed)
!   rsafmn  1 / safmin
!   cnt     counts the rescaling steps
!
! ----------------------------------------------------------------------
!
!            --- check for quick exit ---
!
        IF ( N .LE. 1 ) &
        THEN
          TAU = ZERO
        ELSE
!
          XNORM = DNRM2( N - 1, X, INCX )
!
          IF ( XNORM .LE. DRPTOL ) &
          THEN
!
!                --- H = I ---
!
            TAU = ZERO
            CALL DCOPY( N - 1, ZERO, 0, X, INCX )
!
          ELSE
!
!                --- general case ---
!
            BETA = - SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
!
            IF ( ABS( BETA ) .LT. SAFMIN ) &
            THEN
!
!                    --- xnorm, beta may be inaccurate; scale x
!                        and recompute them                     ---
!
              RSAFMN = ONE / SAFMIN
              CNT = 0
!
   10         CONTINUE
                CNT = CNT + 1
                ALPHA = ALPHA * RSAFMN
                CALL DSCAL( N - 1, RSAFMN, X, INCX )
                BETA = BETA * RSAFMN
              IF ( ABS( BETA ) .LT. SAFMIN )     GOTO 10
!
!                  --- now abs( beta ) is at least safmin ---
!
              XNORM = DNRM2( N - 1, X, INCX )
              BETA = - SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
              TAU = ( BETA - ALPHA ) / BETA
              CALL DSCAL( N - 1, ONE / ( ALPHA - BETA ), X, INCX )
!
!                  --- if alpha is subnormal, it may lose
!                      relative accuracy                  ---
!
              ALPHA = BETA
!
              DO 20 J = 1, CNT
                ALPHA = ALPHA * SAFMIN
   20         CONTINUE
!
            ELSE
!
              TAU = ( BETA - ALPHA ) / BETA
              CALL DSCAL( N - 1, ONE / ( ALPHA - BETA ), X, INCX )
              ALPHA = BETA
!
            ENDIF
!
          ENDIF
!
        ENDIF
!
        RETURN
        END
