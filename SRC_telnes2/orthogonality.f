!BOP
! !ROUTINE: Orthogonality
! !INTERFACE:
      SUBROUTINE Orthogonality
! !USES:
      use initialstate1
      use initialstate2
  	  use dimension_constants
	  use struct, only : jrj
	  use input
	  use leftovers
	  use energygrid
	  use radial_functions
	  use program_control,only : headers
! !DESCRIPTION:
!     Orthogonality calculates the product of the initial state wave function
!     and final state radial functions for each orbital quantum number l of the
!     final state and for each eigenenergy of the final state.
!     The results are written to file 46.
!     The integration is done by calling rint13.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none
!     local variables
      INTEGER IRadius, elPr,LPrime,i,j
      DOUBLE PRECISION ScalProd(0:lmax), ScalProd2(0:lmax)
      DOUBLE PRECISION TrueDP(NRAD+41), TrueDQ(NRAD+41)
      DOUBLE PRECISION TrueDPdp(NRAD+41), TrueDQdq(NRAD+41)
      

      write(6,'(/,a)') 'Output from subroutine CalculOrtho :'
!     We put in TrueDP and TrueDQ the value of r*u (DP = r^2 * u, ...)
      DO I=1, JRJ(NATOM)
         TrueDP(I) = DP(I) / DR(I)
         TrueDQ(I) = DQ(I) / DR(I)
      ENDDO
      IF (LC.GT.0) THEN
         DO I=1, JRJ(NATOM)
            TrueDPdp(I) = DPdp(I) / DR(I)
            TrueDQdq(I) = DQdq(I) / DR(I)
         ENDDO
      ENDIF


!     We write in case.ortho the orthogonality integrals

      IF (LC.EQ.0) THEN
	     if (headers) then
         WRITE(46,*) '## Orthogonality integrals, u initial at l=0 j=l+1/2'
         WRITE(46,*) '## Radial function of the end state: w_l.'
         WRITE(46,*) '## Energy, int u.w_0, int u.w_1, int u.w_2, int u.w_3'
		 endif
         DO J=JEMIN, JEMAX, 5
            DO LPrime=0, lmax
               CALL RINT13(CFEIN1, CFEIN2, &
                    TrueDP, TrueDQ, &
                    Uell(1, J, LPrime), UPell(1, J, LPrime), &
                    ScalProd(LPrime), NATOM)
            ENDDO
            WRITE(46,1) ENE(J), (ScalProd(LPrime), LPrime=0, 3)
         ENDDO
      ELSE 
         WRITE(46,*) '## Orthogonality integrals'
         WRITE(46,*) '## Radial function of the initial state at l and' &
                     ,' j=l+1/2 (u), or at j=l-1/2 (v)'
         WRITE(46,*) '## Radial function of the end state: w_l.', &
              'Energy, int u.w_0, int u.w_1, int u.w_2, int u.w_3, ', &
              'int v.w_0, int v.w_1, int v.w_2, int v.w_3'
         DO J=JEMIN, JEMAX, 5
            DO LPrime=0, lmax
               CALL RINT13(CFEIN1, CFEIN2, &
                    TrueDP, TrueDQ, &
                    Uell(1, J, LPrime), UPell(1, J, LPrime), &
                    ScalProd(LPrime), NATOM)
               CALL RINT13(CFEIN1, CFEIN2, &
                    TrueDPdp, TrueDQdq, &
                    Uell(1, J, LPrime), UPell(1, J, LPrime), &
                    ScalProd2(LPrime), NATOM)
            ENDDO
            WRITE(46,1) ENE(J), (ScalProd(LPrime), LPrime=0, 3), &
                 (ScalProd2(LPrime), LPrime=0, 3)
         ENDDO
      ENDIF

 1    FORMAT(F10.5,8E14.5)
 4713 FORMAT(F10.5,7e14.5)

      RETURN
      END


