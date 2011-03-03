!BOP
! !ROUTINE: Identity
! !INTERFACE:
      SUBROUTINE IDENTITY(M)
! !INPUT/OUTPUT PARAMETERS:
!   M  :  3*3 matrix which is set to unity
! !DESCRIPTION:
!   Initialize a 3*3 matrix to unity.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  IN/OUTPUT
      DOUBLE PRECISION M(3,3)
!  LOCALS
      integer i,j
      
      DO I=1,3
	  DO J=1,3
	    IF (I.EQ.J) THEN
	      M(I,J)=dble(1)
	    ELSE
	      M(I,J)=dble(0)
	    END IF
	  ENDDO
      ENDDO

      RETURN
      END

