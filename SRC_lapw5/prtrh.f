      SUBROUTINE PRTRH                                                  
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      real*8,allocatable :: CHG(:,:)                                        
      REWIND 10                                                         
       READ(10) NPX,NPY,RELX,RELY                                       
      allocate ( CHG(Npx,npy))                                        
      DO 1 I=1,NPX                                                      
      DO 1 J=1,NPY                                                      
    1  READ(10) CHG(I,J)                                                
       NPX0=1                                                           
       NPY0=1                                                           
      DO 2 J=1,NPY,NPY0                                                 
    2 WRITE(6,10) (CHG(I,J),I=1,NPX,NPX0)                               
   10 FORMAT(1H0,13(10F12.5,/))                                         
      WRITE(21,11) ((CHG(I,J),J=1,NPY),I=1,NPX)                         
 11   FORMAT(5E16.8)                                                    
      RETURN                                                            
      END                                                               
