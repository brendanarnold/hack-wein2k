!BOP
! !ROUTINE: Rotatecrist
! !INTERFACE:
    subroutine RotateCrist (V, Alpha, Beta, Gamma)
! !INPUT/OUTPUT PARAMETERS:
!   v    :   vector before and after rotation
!   alpha,beta,gamma : Euler angles defining the rotation
! !DESCRIPTION:
!   Rotation of a vector, where the rotation is defined by Euler angles.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  IN/OUTPUT
      real*8,intent(in) :: Alpha, Beta, Gamma
	  real*8 v(3)
!  LOCAL
      real*8 vout(3)
      
      call Rotation (V, Alpha, Beta, Gamma, Vout)
      v=vout

      RETURN
      END
      
