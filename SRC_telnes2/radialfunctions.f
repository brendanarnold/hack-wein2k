!BOP
! !ROUTINE: RadialFunctions
! !INTERFACE:
      SUBROUTINE RadialFunctions
! !USES:
      use energygrid
	  use radial_functions
	  use dimension_constants
	  use work_atpar
	  use input,only : natom, UseThisL
	  use program_control,only : verbosity, headers
	  use struct
	  use initialstate1
! !DESCRIPTION:
!  Calculation of the radial APW basis functions.
!  Calls atpar for doing the actual work.
!  Large work arrays are allocated and destroyed before returning.
!  Some output to file 10 for high verbosity.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none
!     local variables
      INTEGER I,l,JHalf,J

      write(6,'(/,a)') 'Output from subroutine CalculRadial:'
!     Initialise arrays
      call make_abv(nrad,iemax)
	  call make_uel(nrad,iemax,lmax)

      DO l=0, lmax
	  if (UseThisL(l)) then  ! only calculate what we need
!     Calculation of the radial functions for l:
!     Atpar writes in A1(NRAD, IEMAX0) and B1(..) the radial functions Ul'(energy)
         CALL ATPAR(NATOM,l)
         DO I = 1, NRAD
            DO J = JEMIN, JEMAX
               Uell(I, J, l) = A1(I, J)
               UPell(I, J, l) = B1(I, J)
            ENDDO
         ENDDO
	  endif
      ENDDO
      call destroy_abv

      if(verbosity.ge.2) then
!     we write s- p- d- f- radial functions at 3 energies (0, middle, end)
!     The energies are written in the header:
      JHalf = (JEMAX-JEMIN) / 2
      if (headers) WRITE(10,2) ENE(JEMIN), ENE(JEMIN+JHalf), ENE(JEMIN+2*JHalf)
      DO I = 1, JRJ(NATOM)
         WRITE(10,1) DR(I),(Uell(I,Jemin,l),l=0,lmax),(Uell(I,Jhalf+jemin,l),l=0,lmax),(Uell(I,Jemax,l),l=0,lmax)
      ENDDO
	  endif   ! verbosity

      
 1    FORMAT(F10.5, 12E14.5)
 2    FORMAT('# first energy :', F10.5, 'second energy :', F10.5,  &
           'third energy :', F10.5)

      RETURN
      END


