      subroutine eshift(xe,eps,negrid,esh,maxpol,ecv)
      implicit logical (a-z)
      integer i,j,k,l,maxpol,nestart,negrid
      include 'param.inc'
      real*8  eps(maxde,mpol),xe(maxde),esh,x
      character*5 ecv
      write(6,222) esh,ecv
  222 format(' Im(epsilon) shifted by ',F8.4,1x,A5)
      do i=1,negrid
       if (xe(i).ge.esh) goto 100
      end do
  100 nestart=i
!     write(*,*) 'eshift: ',nestart,esh,xe(i)
      do i=negrid,nestart,-1
        k=i-nestart+1
        do j=1,maxpol
        eps(i,j)=eps(k,j)
        end do
      end do
      do i=1,nestart-1
        do j=1,maxpol
        eps(i,j)=0.0
        end do
      end do
      end
