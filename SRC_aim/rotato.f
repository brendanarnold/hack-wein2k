!$$$      SUBROUTINE ROTATO(VT,IZ,TAU)                                    
!$$$C                                                                       
!$$$C     ROTATES THE REAL VECTOR VT ,                                          
!$$$      IMPLICIT REAL*8 (A-H,O-Z)
!$$$      DIMENSION VT(3),IZ(3,3),TAU(3),VN(3)                         
!$$$      DO 10 JC=1,3                                                      
!$$$      DOTPRO=0.0                                                        
!$$$      DO 20 J=1,3                                                       
!$$$  20  DOTPRO=DOTPRO+VT(J)*IZ(J,Jc)                                      
!$$$  10  VN(JC)=DOTPRO+TAU(JC)                                       
!$$$      DO 30 J=1,3                                                       
!$$$  30  VT(J)= VN(J)                                                      
!$$$      RETURN                                                            
!$$$      END                                                               



      subroutine rotato(v,oiz,otau)
      implicit none
!                                                                       
!     ROTATES THE REAL VECTOR VT                                          
!
      real*8 v,oiz,otau,vn,dotpro
      
      integer i,j

      dimension v(3),oiz(3,3),otau(3),vn(3)                         

      do  i=1,3                                                      
        dotpro=0.0D0                                                        
        do j=1,3                                                       
          dotpro=dotpro+v(j)*oiz(j,i)
        enddo
        vn(i)=dotpro+otau(i)
      enddo
      do j=1,3                                                       
        v(j)= vn(j)
      enddo
      RETURN                                                            
      END                                                               

