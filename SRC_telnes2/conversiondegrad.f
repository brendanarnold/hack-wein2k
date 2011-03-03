!BOP
! !ROUTINE: ConversionDegRad
! !INTERFACE:
      real*8 FUNCTION ConversionDegRad(Angle)
! !INPUT/OUTPUT PARAMETERS:
!   Angle   :  angle in (decimal) degrees
! !DESCRIPTION:
!   Converts an angle from decimal degrees to radians.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  INPUT
	  real*8,intent(in) :: Angle

      ConversionDegRad = Angle / 57.29577951308232D0
      RETURN
      END

