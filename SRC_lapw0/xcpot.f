      FUNCTION XCPOT(RHO,RHOSP,PI,IEX)                          
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XCPOT 
      REAL*8 NUCE,MUCP,MUCF                                             
!     IEX=1  SPINPOLARIZED JMW PARAMETERS                               
      DATA C1,C2/-.045D0,-.0225D0/                                      
      DATA A1,A2/21.0D0,52.9167D0/                                      
      DATA T016PS,GAMMA,CMUXP/0.01899721D0,5.1297628D0,1.2217741D0/     
      DATA ONETH,FOURTH/.3333333333333333D0,1.3333333333333333D0/               
      F(X)=(1.0D0+X*X*X)*LOG(1.0D0+1.0D0/X)+0.5D0*X-X*X-ONETH           
      FOFXI(X)=(X**FOURTH+(1.0D0-X)**FOURTH-0.793700526D0)/ &
       0.206299474D0      
!                                                                       
      RS=(0.75D0/PI/RHO)**(1.D0/3.D0)                                   
      BETA=1.0d0                                            
      RHOH=2.0d0*RHOSP
      IF(IEX.EQ.1) GOTO 100                                             
      RHOH=RHO
      IF(IEX.EQ.2) THEN                                                 
      BETA=1.d0+0.0316D0*RS*LOG(1.D0+24.3D0/RS)                           
!.....HEDIN-LUNDQIST APPROXIMATION IN LINE ABOVE                        
!     BETA=1.d0+(.9604D0*RS*(RS+5.85D0)/(RS+7.8D0)**2)                    
!.....WIGNER APPROXIMATION IN LINE ABOVE                                
      ELSE IF(IEX.EQ.3) THEN                                            
      BETA=1.d0+0.03683D0*RS*LOG(1.D0+21.0D0/RS)                          
      ELSE IF(IEX.EQ.4) THEN                                            
      BETA=1.d0+0.0545D0*RS*LOG(1.D0+11.4D0/RS)                           
      ELSE IF((IEX.EQ.5) .OR. (IEX.EQ.6)) THEN
      BETA=1.0d0                                            
      RHOH=2.0d0*RHOSP
      END IF                                                            
      XCPOT=-BETA*2.D0*(3.D0/PI*RHOH)**(1.D0/3.D0)                       
      RETURN                                                            
  100 CONTINUE                                                          
!                                                                       
!.... SPINPOLARIZED POTENTIAL, JMW PARAMETERS                           
!                                                                       
      XI=RHOSP/RHO                                                      
      X=RS/A1                                                           
      MUCP=C1*LOG(1.0d0+1.0D0/X)                                          
      EFMEP=-C1*F(X)                                                    
      X=RS/A2                                                           
      MUCF=C2*LOG(1.0D0+1.0D0/X)                                        
      EFMEP=EFMEP+C2*F(X)                                               
      NUCE=GAMMA*EFMEP                                                  
      TAUCE=MUCF-MUCP-FOURTH*EFMEP                                      
      XCPOT=(-CMUXP/RS+NUCE)*((2.0D0*XI)**ONETH)-NUCE+MUCP               &
        +TAUCE*FOFXI(XI)                                                
      RETURN                                                            
      END                                                               
