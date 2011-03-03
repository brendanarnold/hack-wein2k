      subroutine expandset(pred,pall)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
!     Expand from reduced to full set
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax),pall(3,natmax*48)
      do j=1,index
         n=icode(1,j)
         i=icode(2,j)
!         write(88,*)j,n,i
        if(i.lt.1)then
         do jj=1,3
         pall(jj,j)=pred(jj,n)
         enddo
        else
!       Positions after applying symmetry operations
         pall(1,j) = dble(iz(1,1,i))*pred(1,n) + dble(iz(2,1,i))*pred(2,n) &
         + dble(iz(3,1,i))*pred(3,n) +  &
         tau(1,i) +dble(icode(3,j))
         pall(2,j) = dble(iz(1,2,i))*pred(1,n) + dble(iz(2,2,i))*pred(2,n) &
         + dble(iz(3,2,i))*pred(3,n) +  &
         tau(2,i) +dble(icode(4,j))
         pall(3,j) = dble(iz(1,3,i))*pred(1,n) + dble(iz(2,3,i))*pred(2,n) &
         + dble(iz(3,3,i))*pred(3,n) + &
         tau(3,i) +dble(icode(5,j))
        endif
      enddo
      return
      end

