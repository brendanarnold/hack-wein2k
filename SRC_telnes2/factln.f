!BOP
! !ROUTINE: Factln
! !INTERFACE:
      real*8 FUNCTION FACTLN(N)
! !INPUT/OUTPUT PARAMETERS:
!   n   :  argument
! !DESCRIPTION:
!   Calculates a factorial.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      real*8 AFact(100),GAMMLN
      SAVE AFact
      DATA AFact/100*-1.d0/


      IF (N.LT.0) PAUSE 'negative factorial'
      IF (N.LE.99) THEN
        IF (AFact(N+1).LT.0.d0) AFact(N+1)=GAMMLN(DBLE(N)+1.d0)
        FACTLN=AFact(N+1)
      ELSE
        FACTLN=GAMMLN(DBLE(N)+1.d0)
      ENDIF
      RETURN
      END
