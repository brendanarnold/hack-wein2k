      SUBROUTINE inispl  
      USE param 
      USE abc                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DO 887 NUM=1,NUME
         do i=1,NDIM2
         do j=1,NDIM2
         xqtl(i,j,num)=(0.d0,0.d0)
         enddo
         enddo
 887  CONTINUE                                                          
      END                                                               
