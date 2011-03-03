!BOP
! !ROUTINE: CalculateEnergySpectrum
! !INTERFACE:
      SUBROUTINE CalculateEnergySpectrum
! !USES:
      use constants
	  use dimension_constants,only : lmax
	  use input, only : LC
      use program_control,only : averaged
	  use energygrid, only : jemax
	  use spectra_normal,only : makespec_normal
	  use spectra_hyperfine,only : makespec_hyperfine
      use crash
! !DESCRIPTION:
!    Governs the calculation of energy differential spectra.
!    Depending on calculation type - orientation sensitive or orientation
!    averaged - different routines are used.
!    First the necessary arrays are allocated.
!    Then the scattering cross section is calculated.
!    Finally, output is written.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP

	  implicit none

      
      IF (averaged) THEN
	     call makespec_normal(jemax,lmax,min(3,3*LC+1))    ! min(3,3*LC+1) equals 1 if LC=0, and equals 3 if LC > 0. - A little ugly, I admit it ...
         call errflg(errfn,'Error in AveragedSpectrum')
         call AveragedEnergySpectrum
         call errflg(errfn,'Error while writing output in WriteAveragedSpectrum')
         call WriteAveragedEnergySpectrum
      ELSE
	     call makespec_hyperfine(jemax,lmax,min(3,3*LC+1))
         CALL errflg(errfn,'Error in OrientedSpectrum')
         CALL OrientedEnergySpectrum
		 call errflg(errfn,'Error while writing output in WriteOrientedSpectrum')
		 call WriteOrientedEnergySpectrum
      ENDIF
!     I allocated, I calculated, I wrote.

      RETURN
      END

