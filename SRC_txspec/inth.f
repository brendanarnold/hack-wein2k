      SUBROUTINE INTH (DP,DQ,DV,DR)                                     
!                                                                       
! INTEGRATION PAR LA METHODE DE ADAMS A 5 POINTS DE LA GRANDE COMPOSANTE
! DP ET DE LA PETITE COMPOSANTE DQ AU POINT DR,DV ETANT LE POTENTIEL EN 
! CE POINT                                                              
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PS1/DEP(5),DEQ(5),DB,DVC,DSAL,DK,DM                        
!                                                                       
! DEP,DEQ DERIVEES DE DP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA    
! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA             
! DM=PAS EXPONENTIEL/720.                                               
! DKOEF1=405./502., DKOEF2=27./502.                                     
! ********************************************************************* 
      DATA DKOEF1/.9462151394422310D0/,DKOEF2/.5378486055776890D-1/     
      DPR=DP+DM*((251.*DEP(1)+2616.*DEP(3)+1901.*DEP(5))-(1274.*DEP(2)+2774.*DEP(4)))                                                     
      DQR=DQ+DM*((251.*DEQ(1)+2616.*DEQ(3)+1901.*DEQ(5))-(1274.*DEQ(2)+2774.*DEQ(4)))                                                     
      DO 13 I=2,5                                                       
      DEP(I-1)=DEP(I)                                                   
 13   DEQ(I-1)=DEQ(I)                                                   
      DSUM=(DB-DV/DVC)*DR                                               
      DEP(5)=-DK*DPR+(DSAL*DR+DSUM)*DQR                                 
      DEQ(5)=DK*DQR-DSUM*DPR                                            
      DP=DP+DM*((106.*DEP(2)+646.*DEP(4)+251.*DEP(5))-(19.*DEP(1)+264.*DEP(3)))                                                           
      DQ=DQ+DM*((106.*DEQ(2)+646.*DEQ(4)+251.*DEQ(5))-(19.*DEQ(1)+264.*DEQ(3)))                                                           
      DP=DKOEF1*DP+DKOEF2*DPR                                           
      DQ=DKOEF1*DQ+DKOEF2*DQR                                           
      DEP(5)=-DK*DP+(DSAL*DR+DSUM)*DQ                                   
      DEQ(5)=DK*DQ-DSUM*DP                                              
      RETURN                                                            
      END                                                               
