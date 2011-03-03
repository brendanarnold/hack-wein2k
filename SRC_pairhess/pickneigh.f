      subroutine pickneig(distances,nuse,N2)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension distances(neigmax,natmax*48)
      integer nuse(natmax)
!
!     Just ignore those which are too small
      Mdo=0
      do j=1,N2
         M=nuse(j)
         dmin=distances(1,j)
         do n=2,M
            d=distances(n,j)
            damp=exp(-scl*(d-dmin)/dmin)
            if(damp .lt. cutoff)then
                    Mdo=n-1
                    goto 10
            endif
        enddo
10      nuse(j)=Mdo
      enddo
      return
      end
