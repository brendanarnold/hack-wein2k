!BOP
! !ROUTINE: DblFactrl
! !INTERFACE:
      real*8 FUNCTION DblFactrl(N)
! !INPUT/OUTPUT PARAMETERS:
!     n   :   argument; n!! is calculated.
! !DESCRIPTION:
!    Calculates double factorials
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  INPUT
      integer,intent(in) :: n
!  LOCALS
      integer m,j,mtop
      real*8 AFact(33),FACTLN
      SAVE MTOP,AFact
      DATA MTOP,AFact(1)/0,1.d0/
!---------------N=2M+1... determine M
      M=(N-1)/2
      IF (M.LT.0) THEN
        PAUSE 'negative factorial'
      else if (MOD(N,2).EQ.0) then
	  DBLFACTRL=dble(0)
      ELSE IF (M.LE.MTOP) THEN
        DBLFACTRL=AFact(M+1)
      ELSE IF (M.LE.32) THEN
        do J=MTOP+1,M
          AFact(J+1)=DBLE(2*J+1)*AFact(J)
        enddo
        MTOP=M
        DBLFACTRL=AFact(M+1)
      ELSE
	  DBLFACTRL= DEXP(FACTLN(2*M+1)-FACTLN(M))/2.d0**M
      ENDIF
      RETURN
      END

