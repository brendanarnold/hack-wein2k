!BOP
! !ROUTINE: UseBr1
! !INTERFACE:
    subroutine UseBr1 (krot, k)
! !USES:
    use rotation_matrices,only : br1
! !INPUT/OUTPUT PARAMETERS:
!   k    :   vector before rotation
!   krot :   vector after rotation
! !DESCRIPTION:
!   Transformation of a vector from laboratory basis to carthesian crystal basis
!   by multiplication with Bravais matrix of the crystal.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP

	  implicit none      
      real*8,intent(in)  :: krot(3)
      real*8,intent(out) :: k(3)

      k(1)=krot(1)*br1(1,1) + krot(2)*br1(1,2) + krot(3)*br1(1,3)   
      k(2)=krot(1)*br1(2,1) + krot(2)*br1(2,2) + krot(3)*br1(2,3)   
      k(3)=krot(1)*br1(3,1) + krot(2)*br1(3,2) + krot(3)*br1(3,3)   

      RETURN
      END
