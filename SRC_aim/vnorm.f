      REAL*8 FUNCTION VNORM(V)                                                 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(3)                                                    
!     Small changes by L. D. Marks, January 2006
!      VNORM=0.d0                                                         
!      do i=1,3                                                        
!        vnorm=vnorm+v(i)*v(i)
!      end do
!      VNORM=SQRT(VNORM)                                                 
      vnorm=sqrt(V(1)*V(1)+v(2)*V(2)+V(3)*V(3))
      RETURN                                                            
      END                                                               
