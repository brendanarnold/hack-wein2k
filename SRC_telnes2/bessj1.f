!BOP
! !ROUTINE: BessJ1
! !INTERFACE:
      real*8 function BessJ1(X)
! !INPUT/OUTPUT PARAMETERS:
!   x    :  argument of the Bessel function
! !DESCRIPTION:
!   Calculates the spherical Bessel function of the first kind of order 1.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  IN/OUTPUT
      real*8,intent(in) :: x
      real*8,parameter  :: eps=1.d-6

      if (X.LT.EPS) then
         bessj1=x/dble(3)*(dble(1)-x**2/dble(10))
      else
         bessj1=DSIN(X)/X**2 - DCOS(X)/X
      end if
      return
      end
      
