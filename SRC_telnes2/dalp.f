!BOP
! !ROUTINE: Dalp
! !INTERFACE:
      real*8 FUNCTION DALP (D1,D2,D3,D4)
! !INPUT/OUTPUT PARAMETERS:
!   D1  : INITIALE (N-1)
!   D2  : FINALE (N-1)
!   D3  : INITIALE (N)
!   D4  : FINALE (N) 
! !DESCRIPTION:
! PROCEDE DE PRATT POUR ACCELERER LA CONVERGENCE
! !REVISION HISTORY:
!   Originally taken from SRC_lcore.
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
! INPUT/OUTPUT
      real*8 d1,d2,d3,d4
! LOCAL
      real*8 d

      IF ((D1+D4).EQ.(D2+D3)) then
	    d=dble(0.5)
      else
        D=(D4-D2)/((D1+D4)-(D2+D3))                                       
        IF(D.LT.0.) d=dble(0)                                               
      endif
      DALP=D                                                            

      RETURN                                                            
      END                                                               
