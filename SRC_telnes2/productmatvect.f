!BOP
! !ROUTINE: ProductMatVect
! !INTERFACE:
      SUBROUTINE ProductMatVect (M, Vin, Vout)
! !INPUT/OUTPUT PARAMETERS:
!     M  :   3*3 input matrix
!     Vin :  input vector
!     Vout : output vector, Vout = M Vin
! !DESCRIPTION:
!     Calculates the product of a matrix and a vector.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
	  real*8,intent(in) ::  M(3,3), Vin(3)
	  real*8,intent(out) :: Vout(3)
      INTEGER I, K
      
      DO I=1, 3
         Vout(I) = 0.D0
         DO K=1, 3
            Vout(I) = Vout(I) + M(I,K) * Vin(K)
         ENDDO
      ENDDO
      RETURN
      END
      
