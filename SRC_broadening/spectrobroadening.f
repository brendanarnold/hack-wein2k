      subroutine SpectrometerBroadening(x,y,yend,nimax,delta,S)

      implicit none
!   IN/OUTPUT
      integer,intent(in) :: nimax       ! size of the energy mesh
      real*8, intent(in) :: x(nimax)    ! energy mesh
	  real*8, intent(in) :: delta       ! step size of energy mesh
	  real*8, intent(in) :: y(nimax)    ! unbroadened spectrum
	  real*8, intent(out):: yend(nimax) ! broadened spectrum
	  real*8, intent(in) :: S           ! FWHM of Gaussian broadening
!   LOCALS
      integer i1,i2
	  real*8 sig
      real*8,parameter ::  pi=3.14159265359	

!     now we do spectrometer broadening, if selected
      if (S.gt.0) then
         write(6,*)' Gaussian spectrometer broadening with FWHM ',S,' eV,'

!     S in case.inxs is FWHM, transform S into sigma
        sig = S / dble(2.35)
		write(6,*) ' i.e., sigma ',sig,' eV.'

           do i1=1,nimax
              yend(i1)=dble(0)
              do i2=1,nimax
                 yend(i1)=yend(i1)+y(i2)* &
                      dble(2)*delta/(dsqrt(dble(2)*pi) * sig)* &
                      dexp ( ((X(i1)-X(i2))**2)/(dble(-2)*sig**2) ) 
              enddo                 
           enddo
      else
	    write(6,*) 'Spectrometer broadening set to 0.'
		write(6,*) 'No spectrometer broadening.'
      endif

	  return
	  end