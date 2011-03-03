      subroutine rotate_back(vt,invoiz)
      implicit real*8 (a-h,o-z)
!                                                                       
!     rotates back the real vector vt       
      real*8 vt,vn,invoiz
      integer i,j
      dimension vt(3),vn(3),invoiz(3,3)
      do i=1,3       
         vn(i)=0.0D0                                                        
         do j=1,3                                                       
            vn(i)=vn(i)+vt(j)*invoiz(j,i)
         enddo
      enddo
      do j=1,3                                                       
         vt(j)= vn(j)
      enddo
      return                                                            
      end 
                                                              
