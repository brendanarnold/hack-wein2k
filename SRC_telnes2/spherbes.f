!BOP
! !ROUTINE: SpherBes
! !INTERFACE:
    real*8 function SpherBes(x,lambda)
! !INPUT/OUTPUT PARAMETERS:
!   x      :   argument of the spherical Bessel function
!   lambda :   order of the spherical Bessel function
! !DESCRIPTION:
!   Calculates the spherical Bessel function of the first kind ;
!   I think the algorithm goes back to some num rec.
!   Except for small arguments, other routines are invoked to do
!   the real work.
! !REVISION HISTORY:
!   Created by M. Nelhiebel for his SQQMat program.
!   Adapted by PH Louf 20 11 97
!   Updated November 2004 (Kevin Jorissen)
!EOP

      implicit none
!   INPUT
      integer,intent(in) :: lambda
	  real*8,intent(in) :: x
!   LOCALS
      real*8  bessj1,bessj0,bessj,dblfactrl
      real*8, parameter ::  eps=1.d-6

      if (x.lt.eps) then
         spherbes=x**lambda / dblfactrl(2*lambda+1)
      else
         if (lambda.eq.0) then
	    spherbes=bessj0(x)
         else if (lambda.eq.1) then
            spherbes=bessj1(x)
         else
            spherbes=bessj(lambda,x)
         end if
      end if

      return
      end
