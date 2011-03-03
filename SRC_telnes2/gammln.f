!BOP
! !ROUTINE: GammLn
! !INTERFACE:
      real*8 FUNCTION GAMMLN(XX)
! !INPUT/OUTPUT PARAMETERS:
!   xx       :   argument
! !DESCRIPTION:
!   I don't know what this routine does.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP

      implicit none
!   INPUT
      real*8,intent(in) :: xx
!   LOCALS
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
	  integer j
      SAVE COF,STP
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0, &
          -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      do J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
      enddo
      GAMMLN=TMP+DLOG(STP*SER)
      RETURN
      END
