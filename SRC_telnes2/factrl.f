!BOP
! !ROUTINE: Factrl
! !INTERFACE:
      real*8 FUNCTION FACTRL(N)
! !INPUT/OUTPUT PARAMETERS:
!   n   :  argument
! !DESCRIPTION:
!   Calculates a factorial.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP

      real*8 AFact(33),GAMMLN
      SAVE NTOP,AFact
      DATA NTOP,AFact(1)/0,1.d0/
      IF (N.LT.0) THEN
        PAUSE 'negative factorial'
      ELSE IF (N.LE.NTOP) THEN
        FACTRL=AFact(N+1)
      ELSE IF (N.LE.32) THEN
        do J=NTOP+1,N
          AFact(J+1)=DBLE(J)*AFact(J)
        enddo
        NTOP=N
        FACTRL=AFact(N+1)
      ELSE
        FACTRL=DEXP(GAMMLN(DBLE(N)+1.d0))
      ENDIF
      RETURN
      END

