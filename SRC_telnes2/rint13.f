!BOP
! !ROUTINE: Rint13
! !INTERFACE:
      SUBROUTINE RINT13(C1,C2,A,B,X,Y,S,JATOM)
! !USES:
      use struct, only : ro,jrj,dh
! !INPUT/OUTPUT PARAMETERS:
!   c1,c2   :   derivatives of feinstruktur constants
!   jatom   :   atom whose wave function is to be integrated
!   a,b     :   large and small components of final state radial basis function
!   x,y     :   large and small components of initial state wave function
!   S       :   solution of the integral
! !DESCRIPTION:
!     Perform radial integrals of the form :
!     S = Integral (r=0..muffin tin radius)  [  (c1*A*X+c2*B*Y)  ]  dr
!     Where j is the spherical Bessel function of the first kind,
!     and jatom determines the muffin tin radius, the radial mesh, and Q.
!     The metric $r^2$ is already in the wave functions.
! !REVISION HISTORY:
!     Some guy D.D. Koelling appears to have been involved?                         
!     Updated November 2004 (Kevin Jorissen)
!EOP
	   implicit none
!  IN/OUTPUT
      INTEGER,intent(in) ::            JATOM
	  real*8,intent(in) :: C1, C2, A(*), B(*), X(*), Y(*)
	  real*8,intent(out) :: S
!  LOCALS
      INTEGER            J, J1, JRI
      DOUBLE PRECISION   D, DX, GTEM1, GTEM2, P1, P2, R, R1, RNOT, Z2
      DOUBLE PRECISION   Z4
!        DH(i)   - step size of the logarithmic radial mesh for atom i
!        JRJ(i)  - number of radial (logarithmic) mesh points for atom i
!        RO(i)   - first radial mesh point for atom i

      RNOT = RO(JATOM)
      DX = DH(JATOM)
      JRI = JRJ(JATOM)
      D = EXP(DX)
!
      J = 3 - MOD(JRI,2)
      J1 = J - 1
      R = RNOT*(D**(J-1))
      R1 = R/D
      Z4 = dble(0)
      Z2 = dble(0)

 10   CONTINUE
         GTEM1 = C1*A(J)*X(J)
         GTEM2 = C2*B(J)*Y(J)
         Z4 = Z4 + R*(GTEM1+GTEM2)
         IF (Z4 .EQ. 0.0D+0) Z4 = 2.0D-55
         R = R*D
         J = J + 1
         IF (J .GE. JRI) GOTO 20
         Z2 = Z2 + R*(C1*A(J)*X(J) + C2*B(J)*Y(J))
         R = R*D
         J = J + 1

      GOTO 10
   20 CONTINUE
      P1 = RNOT*(C1*A(1)*X(1) + C2*B(1)*Y(1))
      P2 = R1*(C1*A(J1)*X(J1) + C2*B(J1)*Y(J1))
      S = 2*Z2 + 4*Z4 + R*(C1*A(J)*X(J) + C2*B(J)*Y(J)) + P2
      S = (DX*S+P1)/3.0D+0
      IF (J1 .GT. 1) S = S + 0.5D+0*DX*(P1+P2)

      RETURN
      END
