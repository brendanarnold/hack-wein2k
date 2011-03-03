!BOP
! !ROUTINE: ReadDOS
! !INTERFACE:
      SUBROUTINE ReadDOS
! !USES:
      use dimension_constants,only : lmax
      use densityofstates
	  use energygrid
! !DESCRIPTION:
!     Read the DOS in the file case.dos.
!     Define the IndexLM array.
!	  Called for averaged spectra.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP
      implicit none
	  integer i,j
      
!     Read partial density of states from *.dos :
      READ (57,4712) IEMAX
	  call make_dos(iemax,lmax)
	  call make_ene(iemax)

	  do i=1,iemax
        read(57,*,end=913,err=913) ene(i),(dos(i,j),j=0,lmax)
      enddo ! i=1,iemax
      
      RETURN

 913  WRITE(6,*) 'ERROR WHILE READING DOS FILE.'
      stop

 4712 FORMAT(/,38x,I5,/)
      END


