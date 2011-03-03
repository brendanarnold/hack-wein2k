      subroutine  losmo(gamma,xe,epsl,eps,negrid)
      implicit logical (a-z)
      integer i,j,k,l,i12,negrid
      include 'param.inc'
      real*8  eps(maxde),epsl(maxde),xe(maxde) &
              ,gamma,gamma2,pi,x,dx,y(0:mpol-1),ab6 &
              ,sum
      parameter (pi = 3.141592654d0) 
      dx=abs(xe(3)-xe(1))       
      ab6=dx/6.                                      
      gamma2=gamma**2
      do i=2,negrid-2,2
         x=xe(i)                                             
         sum=0.                                              
         do j=1,negrid-2,2
            do l=0,2
            y(l)=eps(j+l)*(1/(gamma2+(x-xe(j+l))**2)) 
            enddo                                               
            sum=sum+(y(0)+4.*y(1)+y(2))   
         enddo                                               
         epsl(i)=gamma*ab6*sum/pi                      
      enddo                                                  
      end  
