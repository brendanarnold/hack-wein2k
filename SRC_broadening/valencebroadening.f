      subroutine ValenceBroadening(X,Y,yend,w,absorb,istep,wshift,E0,E1,E2,EF,delta,nimax)
!     VALENCE BROADENING : the array y is broadened by convolution with a Lorentz-function.
!     The result is in array yend.  Three different broadening schemes are available :
!     - w=0 : the width of the Lorentz does not depend on energy
!     - w=1 : the width of the Lorentz varies linearly with energy
!     - w=2 : the width of the Lorentz varies quadratically with energy
!     - w=3 : the width of the Lorentz is given by the scheme of Boucher et al.

!     Here, energy is assumed to be the content of the array x, offset by EF ('Fermi level').

      implicit none
!  INPUT
       integer,intent(in) :: nimax               ! size of energy mesh
       real*8,intent(in)  :: x(nimax),y(nimax)   ! unbroadened spectrum y on energy mesh x
	   integer,intent(in) :: w                   ! broadening method, w = 0..2
	   logical,intent(in) :: absorb              ! absorption or emission spectrum?
	   integer,intent(in) :: istep               ! offset in pixels
	   real*8,intent(in)  :: wshift              ! offset in eV
	   real*8,intent(in)  :: E0,E1,E2            ! parameters for energy dependent broadening
	   real*8,intent(in)  :: delta               ! energy step of the mesh
	   real*8,intent(in)  :: EF                  ! Fermi level ?
	   real*8,intent(out) :: yend(nimax)         ! broadened spectrum
!  LOCALS
       integer i1,i2
	   real*8 gamma,gamma0
	   real*8,parameter :: pi=3.14159265359



      yend=dble(0)
      if (W.gt.0.and.W.le.3) then
         if(ABSORB) then
            if(W.EQ.1) then
               write(6,*)'Valence broadening with W=E/10'
               write(6,*) 'offset: ',istep, ' channels (',wshift,' eV)'
            elseif(W.EQ.2) then
               write(6,*)'Valence broadening with W=E^2 (Muller et al., Phys.Rev.B 57 8181 (1998))'
			elseif(w.eq.3) then
			   write(6,*) 'Valence broadening according to Moureau et al., Phys.Rev.B 73 195111 (2006)'
			elseif(w.eq.0) then
			   write(6,*) 'Valence broadening set to a tiny tiny constant.'
            endif
         endif
         do i1=1,nimax

!           just not zero...    :-p
            gamma0=dble(0.0000001d0)
            gamma=gamma0
            if(ABSORB) then
!     ABSORPTION PART:
               if (w.eq.1) then
                  if(X(i1).gt.wshift) then
                     gamma=X(i1-istep)/dble(10)	
                  else
                     gamma=gamma0
                  endif
                  write(6,*) X(i1),gamma
               elseif (w.eq.2) then
			      if(X(i1-istep).gt.0) &
                  gamma=pi**2*dsqrt(dble(3))/dble(128)*E1*(X(i1)/E0)**2
                  write(6,*) X(i1),gamma
			   elseif (w.eq.3) then
			      if(X(i1-istep).gt.0) &
			      gamma=dble(0.3906) / (dble(0.056)+dble(123)/X(i1-istep)**dble(2.43))
                  write(6,*) X(i1),gamma
               endif
            else
!     EMISSION PART: 
               if(E0.NE.E2) then 
                  if (X(i1).gt.E0) then
                     gamma=W*(1-((X(i1)-E0)/(EF-E0)))**2
                  elseif (X(i1).gt.E1) then
                     gamma=W
                  else
                     gamma=W+W*(1-((X(i1)-E2)/(E1-E2)))**2
                  endif
               else
                  gamma=W*(1-((X(i1)-E0)/(EF-E0)))**2
               endif
            endif

            gamma=gamma/2
            
            do i2=1,nimax
               
               yend(i2)=yend(i2)+y(i1)/pi* &
              (atan((X(i1)-X(i2)+delta)/gamma) &
              -(atan((X(i1)-X(i2)-delta)/gamma)))

            enddo

         enddo
	  else
!     No valence broadening :
         write(6,*) 'No valence broadening; just copying spectrum.'
		 yend=y
      endif

      return
	  end
