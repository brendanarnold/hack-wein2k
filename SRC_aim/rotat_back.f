      subroutine rotat_back(vt,rotloc)                                       
!                                                                       
!     back-rotates real vector vt with local rotation matrix                 
!     uses transposed matrix, because r(-1) is read in struct           
!                                                                       
      implicit real*8 (a-h,o-z)
      dimension vt(3),rotloc(3,3),vn(3)                                 

      do jc=1,3                                                      
        dotpro=0.0D0                                                        
        do  j=1,3                                                       
          dotpro=dotpro+vt(j)*rotloc(j,jc)                                  
        enddo
        vn(jc)=dotpro
      enddo
      do j=1,3                                                       
        vt(j)= vn(j)                                                      
      enddo
      return                                                            
      end
                                                               
