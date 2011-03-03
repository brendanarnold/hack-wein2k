      subroutine gaulep(x1,x2,x,w,n)
!     Routines to setup quadrature or lattice rules for theta/phi
!
!     Phi terms, just use equal sampling
      implicit real*8 (a-h,o-z)
      dimension x(n),w(n)
      logical trust
      integer igmode
      common /ldmopts/igmode,trust
!     Old version
      if(igmode .eq. 0)then
        call gauleg(x1,x2,x,w,n)
        return
      endif
      pi=acos(-1.D0)
      pi2=2.d0*pi
      top=x2
      bot=x1
10    continue
!     Try and trap multivalued possibilities
      if( top .lt. bot )then
        top=top+2.D0*pi
        goto 10
      endif
      range=top-bot
      step=(top-bot)/dble(n)
      do j=1,n
                w(j)=range/dble(n)
                x(j)=x1+(j-0.5d0)*step
      enddo
      return
      end

      SUBROUTINE gauleg(x1,x2,x,w,n)
!     Gauss-Legendre quadrature for theta values
      implicit real*8 (a-h,o-z)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS,PI
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      logical trust
      integer igmode
      common /ldmopts/igmode,trust
      if(igmode .eq. 2)then
        range=x2-x1
        step=range/dble(n)
        do j=1,n
                w(j)=step
                x(j)=x1+(j-0.5d0)*step
        enddo
         return
      endif
      pi=acos(-1.D0)      
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do i=1,m
        z=cos(PI*(i-0.25d0)/(n+.5d0))
 1      continue
        p1=1.d0
        p2=0.d0
        do j=1,n
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.EPS) goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
      end do
      return
      end
