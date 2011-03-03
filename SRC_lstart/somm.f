      SUBROUTINE SOMM (DR,DP,DQ,DPAS,DA,M,NP)                           
!                                                                       
! INTEGRATION PAR LA METHODE DE SIMPSON DE (DP+DQ)*DR**M DE 0 A R=DR(NP)
! DPAS PAS EXPONENTIEL   POUR R VOISIN DE ZERO (DP+DQ)=CTE*R**DA        
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
!     IMPLICIT REAL (D)                                                 
      DIMENSION DR(*),DP(*),DQ(*)                                       
      MM=M+1                                                            
      D1=DA+MM                                                          
      DA=0.d0                                                             
      DB=0.d0                                                             
      DO 15 I=1,NP                                                      
      DL=DR(I)**MM                                                      
      IF (I.EQ.1.OR.I.EQ.NP) GO TO 5                                    
      DL=DL+DL                                                          
      IF ((I-2*(I/2)).EQ.0) DL=DL+DL                                    
 5    DC=DP(I)*DL                                                       
      IF(DC) 7,9,8                                                      
 7    DB=DB+DC                                                          
      GO TO 9                                                           
 8    DA=DA+DC                                                          
 9    DC=DQ(I)*DL                                                       
      IF(DC) 11,15,13                                                   
 11   DB=DB+DC                                                          
      GO TO 15                                                          
 13   DA=DA+DC                                                          
 15   CONTINUE                                                          
      DA=DPAS*(DA+DB)/3.d0                                                
      DC=EXP(DPAS)-1.d0                                                   
      DB=D1*(D1+1.d0)*DC*EXP((D1-1.d0)*DPAS)                                
      DB=DR(1)*(DR(2)**M)/DB                                            
      DC=(DR(1)**MM)*(1.d0+1.d0/(DC*(D1+1.d0)))/D1                            
      DA=DA+DC*(DP(1)+DQ(1))-DB*(DP(2)+DQ(2))                           
      RETURN                                                            
      END                                                               
