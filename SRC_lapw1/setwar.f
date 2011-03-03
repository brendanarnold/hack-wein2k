      SUBROUTINE SETWAR
!
      use out, only      : IORD, NKK, KZZ1, KKK, IMAT, TAU, TAUP, init_out, &
                           WARP, WARP1, WARP2, WARP3
      use totpot, only   : JA, JB, JC, POTK
      use reallocate, only: doreallocate
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!     ..................................................................
!
!        'SETWAR' sets up the complex array WARP:
!        The array-indices are reciprocal lattice vector coordinates and
!        the values are the according plane wave total potential
!        coefficients.
!        WARP is used in 'WARPIN' (which is called by 'HAMILT') to
!        calculate the warpin terms in the Hamilton-matrix.
!
!     ..................................................................
!
!        Local Scalars
!
      INTEGER            IND, JK, JS, JX, JY, JZ
      INTEGER            K1, K2, K3
      INTEGER            K1old, K2old, K3old
!
!        External Subroutines
!
      EXTERNAL           STERN
!
!        Intrinsic Functions
!
      INTRINSIC          ABS
!
      K1=KMAX1START
      K2=KMAX2START
      K3=KMAX3START
      call init_out(IORD, K1, K2, K3)
 5    CONTINUE
      DO 60 JK = 1, NKK
         KZZ1(1) = JA(JK)
         KZZ1(2) = JB(JK)
         KZZ1(3) = JC(JK)
         CALL STERN(IND,IORD,IMAT,KZZ1,TAU,KKK,TAUP)
         if(jk.lt.10) WRITE (6,6000) KZZ1(1), KZZ1(2), KZZ1(3), IND
         DO 50 JS = 1, IND
            IF (ABS(KKK(1,JS)) .GT. K1) THEN
	       K1old = K1
	       K1 = K1 + K1 / 2
	       call doreallocate(WARP, K1old, K1, K2, K2, K3, K3)
               WARP1=K1
               GOTO 5
            ENDIF
            IF (ABS(KKK(2,JS)) .GT. K2) THEN
	       K2old = K2
	       K2 = K2 + K2 / 2
	       call doreallocate(WARP, K1, K1, K2old, K2, K3, K3)
               WARP2=K2
               GOTO 5
            ENDIF
            IF (ABS(KKK(3,JS)) .GT. K3) THEN
	       K3old = K3
	       K3 = K3 + K3 / 2
	       call doreallocate(WARP, K1, K1, K2, K2, K3old, K3)
               WARP3=K3
               GOTO 5
            ENDIF
            WARP(KKK(1,JS),KKK(2,JS),KKK(3,JS)) = POTK(JK)*conjg(TAUP(JS))
            IF (JK .LT. 9) THEN
               WRITE(6,6040) JS,KKK(1,JS),KKK(2,JS),KKK(3,JS),TAUP(JS), &
                             WARP(KKK(1,JS),KKK(2,JS),KKK(3,JS))
            ENDIF
 50      CONTINUE
 60   CONTINUE
!
      RETURN
!
 6000 FORMAT(3X,'K=',3I5,'  IND=',I2)
 6040 FORMAT(20X,I2,'. WAVE=',3I5,'    TAUP=',2F10.5,/,48X,'WARPING=', &
             2F10.5)
!
!        End of 'SETWAR'
!
      END
