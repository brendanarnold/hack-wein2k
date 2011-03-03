!BOP
! !ROUTINE: Fact
! !INTERFACE:
      REAL*8 FUNCTION FACT(N)
! !INPUT/OUTPUT PARAMETERS:
!   n   :  argument
! !DESCRIPTION:
!   Calculates a factorial.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  INPUT
	  integer,intent(in) :: n
!  LOCALS
      INTEGER N1
      REAL*8 FA
      
      FA = dble(1)
      IF (N.GT.0) THEN
         DO N1 = 1, N
            FA = FA * DBLE(N1)
         ENDDO
      ENDIF
      FACT = FA
      RETURN
      END

