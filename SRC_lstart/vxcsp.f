      FUNCTION VXCSP(RHO,RHOSP,iex)                                         
!                                                                       
!     VXCSP EVALUATES THE SPIN-POLARIZED EXCHANGE CORRELATION           
!     POTENTIAL                                                         
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NUCE,MUCP,MUCF                                             
      DATA ONTHRD,FOTHRD/.333333333333D0,1.333333333333D0/
      DATA TO16PS/0.01899721d0/                                           
      DATA GAMMA,CMUXP/5.1297628d0,1.2217741d0/                             
      DATA C1,C2/-.045d0, -.0225d0/                                         
      DATA A1,A2/21.0d0, 52.9167d0/                                         
      F(X)=(1.0d0+X*X*X)*LOG(1.0d0+1.0d0/X)+.5d0*X-X*X-ONTHRD                   
      FOFXI(X)=(X**FOTHRD+(1.0d0-X)**FOTHRD-.793700526d0)/.206299474d0        
      IF(RHO.GT.1.E-18) GO TO 1                                         
      VXCSP=0.0d0                                                         
      RETURN                                                            
    1 if (iex.eq.5) then
         pi=acos(-1.d0)
         VXCSP=-(0.75D0/PI/PI*2.0d0*RHOSP)**(1.D0/3.D0)
      else
         RS=(3.d0/RHO)**ONTHRD                                               
         XI=RHOSP/RHO                                                      
         X=RS/A1                                                           
         MUCP=C1*LOG(1.0d0+1.0d0/X)                                            
         EFMEP=-C1*F(X)                                                    
         X=RS/A2                                                           
         MUCF=C2*LOG(1.0d0+1.0d0/X)                                            
         EFMEP=EFMEP+C2*F(X)                                               
         NUCE=GAMMA*EFMEP                                                  
         TAUCE=MUCF-MUCP-FOTHRD*EFMEP                                      
         VXCSP=(-CMUXP/RS+NUCE)*((2.0d0*XI)**ONTHRD)-NUCE+MUCP                &
          +TAUCE*FOFXI(XI)                                                 
         VXCSP=VXCSP/2.d0                                                    
      endif
      RETURN                                                            
      END                                                               
