      SUBROUTINE MATMM2(A,B,C)   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
!
      DIMENSION        A(3,3),B(3,3),C(3,3)   
!                           
!     (MULTIPLY 3 X 3 MATRICES)                                         
!     FORM C = A B, WHERE C MAY OVERLAP WITH EITHER A OR B, OR BOTH,    
!     SINCE THE PRODUCT IS DEVELOPED IN A TEMPORARY MATRIX.             
!     (06-JUL-80)                                                       
! 
      C(1,1) = A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      C(1,2) = A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      C(1,3) = A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      C(2,1) = A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      C(2,2) = A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      C(2,3) = A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      C(3,1) = A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      C(3,2) = A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      C(3,3) = A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)                                            
!
      RETURN                                                        
      END     

