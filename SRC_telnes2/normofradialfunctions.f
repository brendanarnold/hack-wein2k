!BOP
! !ROUTINE: NormOfRadialFunctions
! !INTERFACE:
      SUBROUTINE NormOfRadialFunctions
! !USES:
      use energygrid
	  use leftovers
	  use radial_functions
	  use input,only : natom, UseThisL
	  use program_control,only : verbosity
	  use dimension_constants,only : lmax
! !DESCRIPTION:
!     Routine Radialfunctions has calculated the APW radial basis functions
!     for all l-values and all energies.
!     This routine calculates the norm of each function, simply by multiplying
!     each by itself (this is done by rint13).
!     A little output (every thirtieth energy) is written to the main output file 6.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none
      
!     Local variables:
      INTEGER L1, L2, J
      DOUBLE PRECISION Sol
      
	  write(6,'(/,a)') 'Output from routine CalculNormeInt:'
	  write(6,'(a)') 'Calculating products of radial functions for the final states.'

      DO J=JEMIN, JEMAX
         IF (MOD(J, 30).EQ.1.and.verbosity.ge.2) write(6,'(a,f9.3,a)') 'Energy =', ENE(J), ':'
         DO L1 = 0,lmax
            DO L2 = 0, L1
			if(UseThisL(L1).and.UseThisL(L2)) then  ! only calculate what we need

               CALL RINT13(CFEIN1, CFEIN2,  &
                    Uell(1, J, L1), UPell(1, J, L1),  &
                    Uell(1, J, L2), UPell(1, J, L2), Sol, NATOM)  
               IF (ABS(Sol).LT.1.D-12) write(6,1) J, L1, L2, abs(Sol)
               IF (MOD(J, 30).EQ.1.and.verbosity.ge.2) write(6,2) L1, L2, Sol
               UellUell(J, L1, L2) = Sol
               UellUell(J, L2, L1) = Sol   !     because of the symmetry in L1, L2

			endif
            ENDDO
         ENDDO
      ENDDO
      
 1    FORMAT ('JENE=',i3,' : WARNING: | U', I1, ' . U', I1, &
      	 ' | = ', F10.8,'< 10^-12')
 2    FORMAT (' U', I1, ' . U', I1, ' = ', F10.4)
      
      RETURN
      END


 
