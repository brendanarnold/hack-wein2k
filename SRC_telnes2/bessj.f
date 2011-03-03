!BOP
! !ROUTINE: BessJ
! !INTERFACE:
      real*8 FUNCTION BessJ(N,X)
! !INPUT/OUTPUT PARAMETERS:
!   n    :  order of the Bessel function
!   x    :  argument of the Bessel function
! !DESCRIPTION:
!   Calculates the spherical Bessel function of the first kind.
!   Uses Leeb-recursion.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  INPUT
      integer,intent(in) :: n
	  real*8,intent(in) :: x
!  LOCALS
      integer,parameter :: iacc=100
	  real*8 bj,bjm,a,bjp
	  real*8,external :: bessj0,bessj1
	  integer l,m,j


      IF(N.LT.2)PAUSE 'bad argument N in BESSJ'
      IF(X.GT.DBLE(N))THEN
        BJM=BESSJ0(X)
        BJ=BESSJ1(X)
        do J=1,N-1
          BJP=DBLE(2*J+1)/X*BJ-BJM
          BJM=BJ
          BJ=BJP
        enddo
        BESSJ=BJ
      ELSE
        M=2*((N+INT(DSQRT(DBLE(IACC*N))))/2)
        BESSJ=0.d0
!-----------use Leeb-rekursion: Numerische Verfahren der Kernphysik
!-----------rekuriere hinunter mit jN=1, jN+1=1, speichere dabei die
!-----------pseudo jL und normiere am Schluss die pseudo-j0 auf die 
!-----------exacte j0. jL mal der Normkonstante liefert
!-----------dann die exakte Loesung
	bj=dble(1)
	bjp=dble(1)
	do l=M,1,-1
	  bjm=DBLE(2*l+1)/x*bj-bjp
	  if (l.EQ.N) bessj=bj	! keep the not normalized value for bessj
	  bjp=bj
	  bj=bjm
  	enddo
	A=bjm/BESSJ0(x)
	bessj=bessj/A
      END IF
      RETURN
      END
