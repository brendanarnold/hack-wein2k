      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit none
      include 'param.inc'

      INTEGER n
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
!      Given arrays x(1:n) and y(1:n) containing a tabulated function, i
!        .e., y i = f(xi), with x1 < x2 < :: : < xN , and given values
!        yp1 and ypn for the  rst derivative of the inter- polating
!        function at points 1 and n, respectively, this routine returns
!        an array y2(1:n) of length n which contains the second
!        derivatives of the interpolating function at the tabulated
!        points xi. If yp1 and/or ypn are equal to 1   10 30 or larger,
!        the routine is signaled to set the corresponding boundary
!        condition for a natural spline, with zero second derivative on
!        that boundary. Parameter: NRAD is the largest anticipated value
!        of n.
      INTEGER i,k 
      REAL*8 p,qn,sig,un,u(NRAD) 

      if (yp1.gt..99e30) then 
        y2(1)=0.D0 
        u(1)=0.D0 
      else 
        y2(1)=-0.5D0 
        u(1)=(3.D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif 
      do i=2,n-1 
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1)) 
        p=sig*y2(i-1)+2.D0
        y2(i)=(sig-1.)/p 
        u(i)=(6.D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p 
      enddo 
      if (ypn.gt..99e30) then 
        qn=0.D0 
        un=0.D0 
      else 
        qn=0.5D0 
        un=(3.D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      endif 
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.D0) 
      do k=n-1,1,-1 
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo 
      return 
      END


      SUBROUTINE splint(xa,ya,y2a,n,x,ix,y,yp,ypp,sy,syp,sypp) 
      implicit none

      INTEGER n 
      REAL*8 x,y,yp,ypp,xa(n),y2a(n),ya(n) 

      INTEGER k,khi,klo,ix
      REAL*8 a,b,h,a2,b2,inv6
      parameter (inv6=1.D0/6.D0)
      logical sy,syp,sypp
      
!      inv6=1.D0/6.D0
      if (ix.gt.0) then
        klo=ix
        khi=ix+1
      else
        klo=1 
        khi=n
 1      if (khi-klo.gt.1) then 
          k=(khi+klo)/2 
          if(xa(k).gt.x)then 
            khi=k 
          else 
            klo=k 
          endif 
          goto 1 
       endif
       write(6,*) 'klo= ',klo,' khi=',khi
      endif
      h=xa(khi)-xa(klo) 
      if (h.eq.0.) stop 'bad xa input in splint' 
      a=(xa(khi)-x)/h 
!      b=(x-xa(klo))/h 
      b=1.D0-a
      a2=a*a
      b2=b*b
      if(sy) then
        y=a*ya(klo)+b*ya(khi)+  &
           ((a2*a-a)*y2a(klo)+(b2*b-b)*y2a(khi))*(h*h)*inv6
      endif
      if(syp) then
        yp=(ya(khi)-ya(klo))/h+(-(3.D0*a2-1.D0)*y2a(klo)+ &
           (3.D0*b2-1.D0)*y2a(khi))*h*inv6
      endif
      if(sypp) then
        ypp=a*y2a(klo)+b*y2a(khi)
      endif

      return 
      END

