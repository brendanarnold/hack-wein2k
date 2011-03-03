      SUBROUTINE ROTATE(VT,IZ,TAU,A)                                    
!                                                                       
!     ROTATES THE REAL VECTOR VT                                        
!  ... 11.1.94 for srtio3   ic(jc,j) auf iz(j,jc) changed
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VT(3),IZ(3,3),TAU(3),VN(3),A(3)                         
      DO 10 JC=1,3                                                      
      DOTPRO=0.0                                                        
      DO 20 J=1,3                                                       
  20  DOTPRO=DOTPRO+VT(J)*IZ(J,Jc)                                      
  10  VN(JC)=DOTPRO+TAU(JC)*A(JC)                                       
      DO 30 J=1,3                                                       
  30  VT(J)= VN(J)                                                      
      RETURN                                                            
      END                                                               
