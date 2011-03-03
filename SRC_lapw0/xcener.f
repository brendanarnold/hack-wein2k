      FUNCTION  XCENER(RHO,RHOSP,PI,IEX)                         
!                                                                       
!   IEX=1  SPINPOLARIZED JMW PARAMETERS                                 
!   IEX=2  HL                                                           
!   IEX=3  HL (OLD,FREEMAN)                                             
!   IEX=4  VBH (RPA)                                                    
!                                                                       
!.....SEE VASHISHTA,SINGWI PR B6,875                                    
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XCENER 
      DATA C1,C2/-.045D0,-.0225D0/                                      
      DATA A1,A2/21.0D0,52.9167D0/                                      
      DATA EPSX0,EPSX1/.9163305866D0,.23817361D0/                       
      DATA ONETH,FOURTH/.3333333333333333D0,1.3333333333333333D0/               
      F(X)=(1.0D0+X*X*X)*LOG(1.0D0+1.0D0/X)+0.5D0*X-X*X-ONETH           
      FOFXI(X)=(X**FOURTH+(1.0D0-X)**FOURTH-0.793700526D0) &
         /0.206299474D0     
!                                                                       
      RS=(0.75D0/PI/RHO)**(1.D0/3.D0)                                   
      EC=0.0d0                                            
      RHOH=2.0d0*RHOSP
         IF(IEX.EQ.1) GOTO 100                                          
      RHOH=RHO
      IF(IEX.EQ.2) THEN                                                 
      EC=-0.112D0+0.0335D0*LOG(RS)                                      
!.....HEDIN-LUNDQIST APPROXIMATION IN LINE ABOVE                        
!     EC=-0.88D0/(RS+7.8D0)                                             
!.....WIGNER APPROXIMATION IN LINE ABOVE                                
      ELSE IF(IEX.EQ.3) THEN                                            
      EC=-0.045D0*((1.D0+RS**3/21.D0**3)*LOG(1.D0+21.D0/RS)+             &
        RS/2.D0/21.D0-RS*RS/21.D0/21.D0-1.D0/3.D0)           
      ELSE IF(IEX.EQ.4) THEN                                            
      EC=-0.0504D0*((1.D0+RS**3/30.D0**3)*LOG(1.D0+30.D0/RS)+            &
        RS/2.D0/30.D0-RS*RS/30.D0/30.D0-1.D0/3.D0)        
      ELSE IF((IEX.EQ.5) .OR. (IEX.EQ.6)) THEN
      EC=0.0d0                                            
      RHOH=2.0d0*RHOSP
      END IF                                                            
      EX=-1.5D0*(3.D0/PI*RHOH)**(1.D0/3.D0)
      XCENER=EX+EC                                                      
      RETURN                                                            
 100  CONTINUE                                                          
!                                                                       
!... SPIN POLARIZED JMW PARAMETERS                                      
!                                                                       
      XI=RHOSP/RHO                                                      
      X=RS/A1                                                           
      EPSCP=C1*F(X)                                                     
      X=RS/A2                                                           
      EPSCF=C2*F(X)                                                     
      XCENER=-EPSX0/RS+EPSCP+FOFXI(XI)*(EPSCF-EPSCP-EPSX1/RS)           
      END                                                               
