      subroutine losmo1(gamma,xe,eps1,eps2,negrid,deltaxy)
!ad
!ad   Calculation of Re(epsilon) with finite damping
!ad
      implicit logical (a-z)
      integer i,j,k,l,negrid
!ad
      include 'param.inc'
!ad
      real*8  eps1(maxde),eps2(maxde),xe(maxde),gamma,gamma2,y(0:2)
      real*8  pi,x,dx,xm,xp,ab6,sum,deltaxy
!ad
      parameter (pi = 3.141592654d0) 
!ad
      dx=abs(xe(3)-xe(1))  
      ab6=dx/6.d0                                           
      gamma2=gamma**2              
!ad
      do i=2,negrid-2,2
      x=xe(i)                                             
      sum=0.d0                                             
      do j=1,negrid-2,2
         do l=0,2
         xm=xe(j+l)-x
         xp=xe(j+l)+x
         y(l)=eps2(j+l)*(xm/(gamma2+(xm)**2)+xp/(gamma2+(xp)**2))
         enddo
      sum=sum+ab6*(y(0)+4.*y(1)+y(2))    
      enddo                                               
      eps1(i)=deltaxy+sum/pi                      
      enddo    
!ad                                           
      end
