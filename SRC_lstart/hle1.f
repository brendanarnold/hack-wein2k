      FUNCTION HLE1(SIGMA,R,IEX)                                        
!.....HEDIN-LUNDQUIST APPROXIMATION,   HARTREE UNITS |||||||            
      IMPLICIT REAL*8 (A-H,O-Z)
      HLE1=0.                                                           
      IF(SIGMA.LT.1.E-30) RETURN                                        
      PI=ACOS(-1.)                                                      
      RHO=SIGMA/R/R/4./PI                                               
      RS=(0.75/PI/RHO)**(1./3.)                                         
      EX=-1.5*(3./PI*RHO)**(1./3.)                                      
      EC=0.                                                             
      IF(IEX.EQ.2)EC=-0.112+0.0335*LOG(RS)                              
      IF(IEX.EQ.3)EC=-0.045*((1.+RS**3/21.**3)*LOG(1.+21./RS)+RS/2./21.  &
      -RS*RS/21./21.-1./3.)                                             
      HLE1=(EX+EC)/2.                                                   
      RETURN                                                            
      END                                                               
