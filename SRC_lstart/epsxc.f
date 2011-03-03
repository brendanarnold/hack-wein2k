      FUNCTION EPSXC(RS,XI)                                             
!   EPSXC EVALUATES THE EXCHANGE-CORRELATION ENERGY DENSITY             
!   --- SPIN POLARIZED VERSION ---                                      
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA C1,C2/-.045, -.0225/                                         
      DATA A1,A2/21.0, 52.9167/                                         
      DATA EPSX0,EPSX1/.9163305866, .23817361/                          
      DATA ONTHRD,FOTHRD/.3333333333, 1.3333333333/                     
      F(X)=(1.0+X*X*X)*LOG(1.0+1.0/X)+.5*X-X*X-ONTHRD                   
      FOFXI(X)=(X**FOTHRD+(1.0-X)**FOTHRD-.793700526)/.206299474        
      X=RS/A1                                                           
      EPSCP=C1*F(X)                                                     
      X=RS/A2                                                           
      EPSCF=C2*F(X)                                                     
      EPSXC=-EPSX0/RS+EPSCP+FOFXI(XI)*(EPSCF-EPSCP-EPSX1/RS)            
      EPSXC=EPSXC/2.                                                    
      RETURN                                                            
      END                                                               
