      FUNCTION DALP (D1,D2,D3,D4)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DALP
!                                                                       
! PROCEDE DE PRATT POUR ACCELERER LA CONVERGENCE                        
! D1 INITIALE (N-1)   D2 FINALE (N-1)   D3 INITIALE (N)   D4 FINALE (N) 
! **********************************************************************
      IF ((D1+D4).EQ.(D2+D3)) GO TO 7                                   
      D=(D4-D2)/((D1+D4)-(D2+D3))                                       
      IF(D.LT.0.) GO TO 8                                               
      IF(D.LT.0.5) GO TO 9                                              
 7    D=0.5                                                             
      GO TO 9                                                           
 8    D=0.                                                              
 9    DALP=D                                                            
      RETURN                                                            
      END                                                               
