!BOP
! !ROUTINE: Transpose
! !INTERFACE:
    subroutine Transpose(M, Mout)
! !INPUT/OUTPUT PARAMETERS:
!   M    :   3*3 real*8 matrix to be transposed
!   Mout :   transposed matrix
! !DESCRIPTION:
!   Calculate the transpose of a matrix M -> Mout .
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      
      implicit none      
      real*8,intent(in)  :: M(3,3)
      real*8,intent(out) :: Mout(3,3)
      
      Mout(1:3,1) = M(1,1:3)
      Mout(1:3,2) = M(2,1:3)
      Mout(1:3,3) = M(3,1:3)
      
      RETURN
      END
      
      
