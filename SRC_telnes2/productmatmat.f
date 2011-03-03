!BOP
! !ROUTINE: ProductMatMat
! !INTERFACE:
      SUBROUTINE ProductMatMat (M1, M2, M)
! !INPUT/OUTPUT PARAMETERS:
!     M1,M2  :  3*3 input matrices
!     M      :  3*3 output matrix, M = M1 * M2
! !DESCRIPTION:
!     Calculates the product of two 3*3 matrices.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
      real*8,intent(in)  :: M1(3,3),M2(3,3)
	  real*8,intent(out) :: M(3,3)
      INTEGER I, J, K
      
      DO I=1, 3
         DO J=1, 3
            M(I,J) = 0.D0
            DO K=1, 3
               M(I,J) = M(I,J)+M1(I,K) * M2(K,J)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END

