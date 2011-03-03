      SUBROUTINE STERN(NST,IORD,IMAT,KZZ1,TAU,KKK,TAUP)
!
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            IORD, NST
      INTEGER            KZZ1(3), KKK(3,IORD), IMAT(3,3,IORD)
      DOUBLE PRECISION   TAU(3,IORD)
      COMPLEX*16         TAUP(IORD)
!
!     ..................................................................
!
!        Stern generates the star of rec lattice vector KZZ1(I).
!        The star vectors are stored in KKK, the star-size in NST,
!        IMAT contains the symmetry-matrices.
!
!     ..................................................................
!
!        Parameters
!
      COMPLEX*16         IMAG
      PARAMETER          (IMAG = (0.0D+0,1.0D+0))
!
!        Locals
!
      INTEGER            I, J, K, L, M
      INTEGER            G(3), INDEX(NSYM)
      DOUBLE PRECISION   TPI, TK
!
!        Intrinsic Functions
!
      INTRINSIC          ATAN, DCMPLX, EXP
!
      TPI = 8.0D+0*ATAN(1.0D+0)
      G(1) = KZZ1(1)
      G(2) = KZZ1(2)
      G(3) = KZZ1(3)
      NST = 0
!
!         start loop over all symmetry operations
!
      DO 60 I = 1, IORD
         TK = 0.0D+0
         DO 20 J = 1, 3
            TK = TK + TAU(J,I)*G(J)*TPI
            K = 0
            DO 10 L = 1, 3
               K = IMAT(J,L,I)*G(L) + K
   10       CONTINUE
!cc               tk = tk + tau(j,i)*k*tpi
            KKK(J,I) = K
   20    CONTINUE
         IF (NST .NE. 0) THEN
!
!        proof, if the vector KKK(J,I) is a new starmember or not
!
            DO 40 M = 1, NST
               DO 30 J = 1, 3
                  IF (KKK(J,M) .NE. KKK(J,I)) GOTO 40
   30          CONTINUE
!
!        KKK(J,I) is not a new starmember, it is equiv to KKK(J,M).
!        But the tauphase of KKK(J,I) can be new.  Therefore the
!        already defined phase TAUP(M) is averaged with the phase
!        of KKK(J,M).
!
!               TAUP(M) = TAUP(M) + EXP(DCMPLX(TK)*IMAG)
!_REAL         TAUP(M) = TAUP(M) + COS(TK)
!_COMPLEX         TAUP(M) = TAUP(M) + DCMPLX(COS(TK),SIN(TK))
               INDEX(M) = INDEX(M) + 1
               GOTO 60
   40       CONTINUE
         ENDIF
!
!        new vector found
!
         NST = NST + 1
         DO 50 J = 1, 3
            KKK(J,NST) = KKK(J,I)
   50    CONTINUE
!         TAUP(NST) = EXP(DCMPLX(TK)*IMAG)
!_REAL         TAUP(NST) = COS(TK)
!_COMPLEX         TAUP(NST) = DCMPLX(COS(TK),SIN(TK))

         INDEX(NST) = 1
   60 CONTINUE
      DO 70 I = 1, NST
         TAUP(I) = TAUP(I)/INDEX(I)
   70 CONTINUE
!
      RETURN
!
!        End of 'STERN'
!
      END
