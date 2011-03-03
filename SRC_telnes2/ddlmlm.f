!BOP
! !ROUTINE: DDlmLM
! !INTERFACE:
      COMPLEX*16 FUNCTION DDlmLM (IEnergy, l1, m1, s1, l2, m2, s2)
! !USES:
	  use cross_dos,only : xdos
	  use dimension_constants,only : lmax
! !INPUT/OUTPUT PARAMETERS:
!     IEnergy           : integer specifying a particular module
!     l1,m1,s1,l2,m2,s2 : cross-dos component to be calculated
! !DESCRIPTION:
!    The cross-dos is stored in a compact xdos-array.  This function sort of
!    interfaces between array and user; it picks the elements of the xdos array
!    corresponding to a requested energy and quantum numbers l,m,l',m'.
!    That is, the function returns
!       D(l1,m1,s1;E) D(l2,m2,s2;E)*
!    where E is the IEnergy-th value in the energy grid,
!    and the final state wave function psi is expressed as
!    psi(r,theta,phi) = SUM (l=0,infinity; m=-l,+l)   [ u(l;r;E) D(l,m;E) Y(l,m;theta,phi) ]
!    where u are the radial eigenfunctions calculated in RadialFunctions, and Y are
!    spherical harmonics.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!   Added spin indices February 2005 (Kevin Jorissen)
!EOP
	  implicit none
!   INPUT      
      integer,intent(in) :: l1, m1, s1, l2, m2, s2, IEnergy
!   LOCALS
      integer I1, I2, IDD, b

      if(s1.eq.1.and.s2.eq.1) then
	    b = 0
	  elseif((s1.eq.1.and.s2.eq.2).or.(s1.eq.2.and.s2.eq.1)) then
	    b = ((lmax+1)**2*((lmax+1)**2+1)/2)   ! b = 136 for lmax=3
	  elseif(s1.eq.2.and.s2.eq.2) then
	    b = ((lmax+1)**2*((lmax+1)**2+1)/2) + (lmax+1)**4  ! b = 136 + 256 = 392 for lmax=3
	  else
	    stop 'unanticipated s1 or s2 in function ddlmlm'
	  endif
	  
      
      I1 = (l1+1)*(l1+1) - l1 + m1
      I2 = (l2+1)*(l2+1) - l2 + m2
      
      IF (I1.GE.I2) THEN
         IDD = I1*(I1-1)/2 + I2 + b
         DDlmLM =  xdos(ienergy,idd)
      ELSE
         IDD = I2*(I2-1)/2 + I1 + b
         DDlmLM = dconjg(xdos(ienergy,idd))
      ENDIF
      
      RETURN
      END


