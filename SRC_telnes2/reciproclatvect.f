!BOP
! !ROUTINE: ReciprocLatVect
! !INTERFACE:
      SUBROUTINE ReciprocLatVect (V1, V2, Volume, Vout)
! !INPUT/OUTPUT PARAMETERS:
!   v1  :  input vector
!   v2  :  input vector
!   volume : volume
!   vout : output vector
! !DESCRIPTION:
!     calculate 2 Pi * V1 x V2 / Volume = Vout
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
	  real*8,intent(in) ::  V1(3), V2(3), Volume
	  real*8,intent(out) :: Vout(3)
      
      CALL ProdVectVect (V1, V2, Vout)
      Vout(1) = Vout(1) * 6.28318530718D0 / Volume
      Vout(2) = Vout(2) * 6.28318530718D0 / Volume
      Vout(3) = Vout(3) * 6.28318530718D0 / Volume
      
      RETURN
      END


