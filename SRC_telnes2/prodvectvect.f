!BOP
! !ROUTINE: ProdVectVect
! !INTERFACE:
      SUBROUTINE ProdVectVect(V1, V2, Vout)
! !INPUT/OUTPUT PARAMETERS:
!     V1 :  input vector
!     V2 :  input vector
!     Vout : output vector, Vout \= V1 \x V2
! !DESCRIPTION:
!     Calculates the vector product V1 x V2.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
	  real*8,intent(in)  :: v1(3),v2(3)
	  real*8,intent(out) :: vout(3)
      
      Vout(1) = V1(2) * V2(3) - V1(3) * V2(2)
      Vout(2) = V1(3) * V2(1) - V1(1) * V2(3)
      Vout(3) = V1(1) * V2(2) - V1(2) * V2(1)
      
      RETURN
      END
      
