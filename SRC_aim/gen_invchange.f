      subroutine gen_invchange(change,invchange)
      implicit none
      
      real*8 change(3,3),invchange(3,3),y(3,3)

      integer ipiv(3),id,info,i,j
      

      do i=1,3
        do j=1,3
          invchange(i,j)=change(i,j)
        enddo
      enddo
      call ludcmp(invchange(1,1),3,3,ipiv,id,info)
      if(info.ne.0) then
        write(6,*) 'Error inverting change:'
        do i=1,3
          write(6,*) (invchange(i,j),j=1,3)
        enddo
        stop 'ERROR INVERTING change'
      endif
      do  i=1,3
        do j=1,3
          y(i,j)=0.D0
        enddo 
        y(i,i)=1.D0
      enddo 
      do j=1,3
        call lubksb(invchange(1,1),3,3,ipiv,y(1,j))
      enddo         
      do i=1,3
        do j=1,3
          invchange(i,j)=y(i,j)
        enddo
      enddo
        
