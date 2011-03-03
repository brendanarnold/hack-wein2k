      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      implicit real*8 (a-h,o-z)
      INTEGER n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)

!     Computes all eigenvalues and eigenvectors of a real symmetric matrix a,
!     which is of size n by n, stored in a physical np by np array. On output,
!     elements of a above the diagonal are destroyed. d returns the
!     eigenvalues of a in its first n elements. v is a matrix with the same
!     logical and physical dimensions as a, whose columns contain, on output,
!     the normalized eigenvectors of a. nrot returns the number of Jacobi
!     rotations that were required.

      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do ip=1,n 
        do iq=1,n
          v(ip,iq)=0.D0
        enddo
        v(ip,ip)=1.D0
      enddo
      do ip=1,n
        b(ip)=a(ip,ip) 
        d(ip)=b(ip)
        z(ip)=0.D0
      enddo
      nrot=0
      do i=1,50
        sm=0.
        do ip=1,n-1 
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
          enddo
        enddo
        if(sm.eq.0.)return 
        if(i.lt.4)then
          tresh=0.2D0*sm/(n*n) !n**2 
        else
          tresh=0.D0 
        endif
        do ip=1,n-1
          do iq=ip+1,n
            g=100.D0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) &
               .and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.D0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h 
              else
                theta=0.5D0*h/a(ip,iq) 
                t=1.D0/(abs(theta)+sqrt(1.D0+theta**2))
                if(theta.lt.0.D0)t=-t
              endif
              c=1.D0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.D0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=iq+1,n 
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              enddo
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              enddo
              nrot=nrot+1
            endif
          enddo
        enddo
        do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip) 
          z(ip)=0.D0 
        enddo
      enddo
!      pause 'too many iterations in jacobi'
      stop 'Too many iterations in jacobi'
      return
      END

