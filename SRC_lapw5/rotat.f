      SUBROUTINE ROTAT(VT,ROTLOC)                                       
!                                                                       
!     ROTATES REAL VECTOR VT WITH LOCAL ROTATION MATRIX                 
!     USES TRANSPOSED MATRIX, BECAUSE R(-1) IS READ IN STRUCT           
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VT(3),ROTLOC(3,3),VN(3)                                 
      DO 10 JC=1,3                                                      
      DOTPRO=0.0                                                        
      DO 20 J=1,3                                                       
  20  DOTPRO=DOTPRO+VT(J)*ROTLOC(JC,J)                                  
  10  VN(JC)=DOTPRO                                                     
      DO 30 J=1,3                                                       
  30  VT(J)= VN(J)                                                      
      RETURN                                                            
      END                                                               
