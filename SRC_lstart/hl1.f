      FUNCTION HL1(SIGMA,R,IEX)                                         
!.....HEDIN-LUNDQUIST APPROXIMATION,   HARTREE UNITS |||||||            
!  IEX=2   PARAMETERS AFTER J.DE.PHYSIQUE                               
!  IEX=3   PARAMETERS AFTER J.PHYS. (OLDER PARAAM.)                     
      IMPLICIT REAL*8 (A-H,O-Z)
      HL1=0.                                                            
      IF(SIGMA.LT.1.E-30) RETURN                                        
      PI=ACOS(-1.)                                                      
      RHO=SIGMA/R/R/4./PI                                               
      RS=(0.75/PI/RHO)**(1./3.)                                         
      BETA=1.                                                           
      IF(IEX.EQ.2)BETA=1.+0.0316*RS*LOG(1.+24.3/RS)                     
      IF(IEX.EQ.3)BETA=1.+0.03683*RS*LOG(1.+21.0/RS)                    
      HL1=-BETA*1.*(3./PI*RHO)**(1./3.)                                 
      RETURN                                                            
      END                                                               
