!BOP
! !ROUTINE: Coreholelifetime
! !INTERFACE: 
      real*8 function CoreholeLifetime (zz,nc,lc,sc)
! !INPUT/OUTPUT PARAMETERS:
!      nc : main quantum number
!      lc : orbital quantum number
!      sc : spin quantum number
!      zz : atomic number
! !DESCRIPTION:
!     Sets core hole lifetime.  Data comes from graphs in
!     K. Rahkonen and K. Krause,
!     Atomic Data and Nuclear Data Tables, Vol 14, Number 2, 1974.
!     output gamach is in eV
! !REVISION HISTORY:
!     Actually taken from feff8.20 setgam routine
!     Need to ask permission for that.
!     Modified for TELNES2 program by Kevin Jorissen Dec 2004
!EOP

      implicit none
!  IN/OUTPUT
      integer,intent(in) :: nc,lc,sc
	real*8,intent(in)  :: zz
!  LOCALS
      real*8 zh(8,16), gamh(8,16)
	integer findihole(7,4)
      integer,external :: locat
      real*8 zk(8), gamkp(8), dy, gamach
	integer i,k,ihole,iedge

!     Note that 0.99 replaces 1.0, 95.1 replaces 95.0 to avoid roundoff trouble.
!     Gam arrays contain the gamma values.
!     We will take log10 of the gamma values so we can do linear
!     interpolation from a log plot.

      data  zh   / 0.99,  10.0, 20.0, 40.0, 50.0, 60.0, 80.0, 95.1, &
                   0.99, 18.0, 22.0, 35.0, 50.0, 52.0, 75.0,  95.1, &
                   0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1, &
                   0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1, &
                   0.99,  20.0, 28.0, 30.0, 36.0, 53.0,  80.0, 95.1, &
                   0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1, &
                   0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1, &
                   0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1, &
                   0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1, &
                   0.99,  30.0, 40.0, 47.0, 50.0, 63.0,  80.0, 95.1, &
                   0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1, &
                   0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1, &
                   0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1, &
                   0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1, &
                   0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0, &
                   0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0/

      data  gamh / 0.02,  0.28, 0.75,  4.8, 10.5, 21.0, 60.0, 105.0, &
                   0.07,  3.9,  3.8,  7.0,  6.0,  3.7,  8.0,  19.0, &
                   0.001, 0.12,  1.4,  0.8,  2.6,  4.1,   6.3, 10.5, &
                   0.001, 0.12, 0.55,  0.7,  2.1,  3.5,   5.4,  9.0, &
                   0.001,  1.0,  2.9,  2.2,  5.5, 10.0,  22.0, 22.0, &
                   0.001,0.001,  0.5,  2.0,  2.6, 11.0,  15.0, 16.0, &
                   0.001,0.001,  0.5,  2.0,  2.6, 11.0,  10.0, 10.0, &
                   0.0006,0.09, 0.07, 0.48,  1.0,  4.0,   2.7,  4.7, &
                   0.0006,0.09, 0.07, 0.48, 0.87,  2.2,   2.5,  4.3, &
                   0.001,0.001,  6.2,  7.0,  3.2, 12.0,  16.0, 13.0, &
                   0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0, &
                   0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0, &
                   0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0, &
                   0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0, &
                   0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9, &
                   0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9/


      data findihole /1, 0, 0, 0, 0, 0, 0, &
	                2, 3, 4, 0, 0, 0, 0, &
	                5, 6, 7, 8, 9, 0, 0, &
	                10,11,12,13,14,15,16/


      if(nc.le.lc) stop 'Invalid quantum numbers in CoreholeLifetime.'
	if(sc.lt.0.or.sc.gt.1) stop 'Invalid quantum numbers in CoreholeLifetime.'

	if(nc.gt.3) then
	  write(6,*) 'No core hole lifetime information for this edge.'
	  write(6,*) 'Using default value of 0.1 eV.'
        gamach = dble(-1)
	else
	    ihole=findihole(max(2*lc+sc,1),nc)
          do i=1,8
             gamkp(i) = log10 (gamh(i,ihole))
             zk(i) = zh(i,ihole)
          enddo
!         Find out between which x points x0 lies
          i = locat (zz, 8, zk)
          k = min( max(i-1/2,1) , 7 )
          call polint( zk(k), gamkp(k), 2, zz, gamach, dy)

      endif

!     Change from log10 (gamma) to gamma
      coreholelifetime = dble(10)** gamach


      return
      end


      function locat (x, n, xx)
      integer  u, m, n
      double precision x, xx(n)

!     Binary search for index of grid point immediately below x.
!     Array xx required to be monotonic increasing.
!     Returns
!     0            x <  xx(1)
!     1            x =  xx(1)
!     i            x =  xx(i)
!     n            x >= xx(n)

      locat = 0
      u = n+1

   10 if (u-locat .gt. 1)  then
         m = (u + locat) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat = m
         endif
         goto 10
      endif

      return
      end





      subroutine polint( xa, ya, n, x, y, dy)
!     draws a polynimial P(x) of order (n-1) through n points.
!     returns y = P(x) and dy - estimate of the error
!     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) pause 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end



