!BOP
! !ROUTINE: BessJ0
! !INTERFACE:
      real*8 function BessJ0(X)
! !INPUT/OUTPUT PARAMETERS:
!   x    :  argument of the Bessel function
! !DESCRIPTION:
!   Calculates the spherical Bessel function of the first kind of order 0.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  IN/OUTPUT
      real*8,intent(in) :: x
	  real*8,parameter  :: eps=1.d-6

      if (X.LT.EPS) then
         bessj0=dble(1)-x**2/dble(6)
      else
         bessj0=DSIN(X)/X
      end if
      return
      end
