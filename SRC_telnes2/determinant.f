!BOP
! !ROUTINE: Determinant 
! !INTERFACE:
      SUBROUTINE Determinant(M, Det)
! !INPUT/OUTPUT PARAMETERS:
!     M  :   3\*3 real input matrix
!     Det :  determinant of M
! !DESCRIPTION:
!     calculate the determinant of a 3\*3 real matrix.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP      

      implicit none
	  real*8,intent(in) :: M(3,3)
	  real*8,intent(out) :: Det
            
      Det = M(1, 1) * M(2, 2) * M(3, 3) &
           + M(1, 2) * M(2, 3) * M(3, 1) &
           + M(2, 1) * M(3, 2) * M(1, 3) &
           - M(3, 1) * M(2, 2) * M(1, 3) &
           - M(2, 1) * M(1, 2) * M(3, 3) &
           - M(1, 1) * M(3, 2) * M(2, 3)
      
      RETURN
      END

      
