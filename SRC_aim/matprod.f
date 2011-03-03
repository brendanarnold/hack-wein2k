      subroutine mat3prod(res,m1,m2,transp1,transp2)
      implicit none
      
      real*8 m1,m2,res
      integer i,j,k
      
      logical transp1,transp2

      dimension m1(3,3),m2(3,3),res(3,3)

      if (.not.(transp1.or.transp2)) then
        do i=1,3
          do j=1,3
            res(i,j) = 0.D0
            do k=1,3
              res(i,j) = res(i,j)+m1(i,k)*m2(k,j)
            enddo
          enddo
        enddo
      else
        if (transp1.and.(.not.transp2)) then
          do i=1,3
            do j=1,3
              res(i,j) = 0.D0
              do k=1,3
                res(i,j) = res(i,j)+m1(k,i)*m2(k,j)
              enddo
            enddo
          enddo
        else
          if (transp2.and.(.not.transp1)) then
            do i=1,3
              do j=1,3
                res(i,j) = 0.D0
                do k=1,3
                  res(i,j) = res(i,j)+m1(i,k)*m2(j,k)
                enddo
              enddo
            enddo
          else
            do i=1,3
              do j=1,3
                res(i,j) = 0.D0
                do k=1,3
                  res(i,j) = res(i,j)+m1(k,i)*m2(j,k)
                enddo
              enddo
            enddo
          endif
        endif
      endif
      
      return
      end

