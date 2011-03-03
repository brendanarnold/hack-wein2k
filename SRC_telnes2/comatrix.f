!BOP
! !ROUTINE: CoMatrix
! !INTERFACE:
      SUBROUTINE CoMatrix(M, Mout)
! !INPUT/OUTPUT PARAMETERS:
!   M   :   3*3 real input matrix
!   Mout:   3*3 real output matrix, comatrix of M
! !DESCRIPTION:
!   Calculates comatrices of 3*3 real matrices.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  IN/OUTPUT
      real*8,intent(in)  ::   M(3,3)
	  real*8,intent(out) ::   Mout(3,3)
!  LOCALS
      INTEGER I, J, I1, I2, J1, J2
      
      DO I=1, 3
         I1 = MOD(I, 3)+1
         I2 = MOD(I+1, 3)+1
         DO J=1, 3
            J1 = MOD(J, 3)+1
            J2 = MOD(J+1, 3)+1
!     we don't have to multiply with (-1)^(I+J) with these definitions of 
!     I1, I2, J1 and J2.
            Mout(I, J) = M(I1, J1) * M(I2, J2) - M(I1, J2) * M(I2, J1)
         ENDDO
      ENDDO
      
      RETURN
      END
      
