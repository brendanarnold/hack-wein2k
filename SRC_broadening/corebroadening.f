       subroutine CoreBroadening(x,y,yend,nimax,delta,g)

      implicit none
!   IN/OUTPUT
      integer,intent(in) :: nimax       ! size of the energy mesh
      real*8, intent(in) :: x(nimax)    ! energy mesh
	  real*8, intent(in) :: delta       ! step size of energy mesh
	  real*8, intent(in) :: y(nimax)    ! unbroadened spectrum
	  real*8, intent(out):: yend(nimax) ! broadened spectrum
	  real*8, intent(in) :: g           ! Lorentzian broadening - core hole width
!   LOCALS
      integer i1,i2
      real*8,parameter ::  pi=3.14159265359	



         if (g.gt.dble(0)) then
            write(6,*) 'Core broadening with ',g, ' eV.'
            do i1=1,nimax
               do i2=1,nimax
                  yend(i2)=yend(i2)+y(i1)/pi* &
                       (atan((X(i1)-X(i2)+delta)/g) &
                       -(atan((X(i1)-X(i2)-delta)/g)))
	         enddo
            enddo
         elseif(g.eq.dble(0)) then
           write(6,*) 'No core hole broadening.'
           yend=y
		 else
		   stop "Error in CoreBroadening - width isn't supposed to be negative."
         endif


         return
         end
