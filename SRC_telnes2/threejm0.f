!BOP
! !ROUTINE: threejm0
! !INTERFACE:
	subroutine threejm0(j1,j2,j3,frac2,jj1)
! !INPUT/OUTPUT PARAMETERS:
!   J1,J2,J3   :   integers specifying the 3j-symbol (J1   J2   J3)
!                                                      0    0    0
!   frac2,jj1  :   specifying the result as (-1)**jj1 * sqrt(frac2)
! !DESCRIPTION:
!   Calculates Wigner 3j symbols for m1=m2=m3=0 using eq. 3.7.9 and 3.7.16 (recursion) from Edmonds.
! !REVISION HISTORY:
!   Created M. Nelhiebel, 970811
!   Updated November 2004 (Kevin Jorissen)
!EOP

    implicit none
!  IN/OUTPUT
	real*8,intent(out)  ::  frac2
	integer,intent(in)  ::  j1,j2,j3
	integer,intent(out) ::  jj1
!  LOCAL VARIABLES
	integer                 jj2,jj3,j,jj,num1,num2,idenom1,idenom2

	frac2=dble(0)
	jj1=j1
	jj2=j2
	jj3=j3
	
	j=j1+j2+j3

	if (mod(j,2).ne.0) return				  !zero if j=odd
	if (.not.((abs(j1-j2).le.j3).and.((j1+j2).ge.j3))) return !triangle rule : 3j-symbol = 0

	frac2=dble(1)
	if (j.eq.0) return

!----------sort j1 j2 j3 in descending order
	
 10	if (jj1.lt.jj2) then
	   jj=jj1
	   jj1=jj2
	   jj2=jj
	end if
	if (jj1.lt.jj3)	then
	   jj=jj1
	   jj1=jj3
	   jj3=jj
	end if
	if (jj2.lt.jj3) then
	  jj=jj2
	  jj2=jj3
	  jj3=jj
	end if

!-----------formula 3.7.9 if jj3=0
	if (jj3.eq.0) then
	  frac2=frac2/dble(2*jj1+1)
	  return
	else
!-----------rekursion 3.7.16
	  num1=j-2*jj2-1
	  num2=j-2*jj3+2
	  frac2=frac2*dble(num1*num2)
	  idenom1=j-2*jj2
	  idenom2=j-2*jj3+1
	  frac2=frac2/dble(idenom1*idenom2)
	  jj3=jj3-1
	  jj2=jj2+1
	end if
	goto 10
	
	end
	