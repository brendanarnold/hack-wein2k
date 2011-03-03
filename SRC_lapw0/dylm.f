      SUBROUTINE DYLM(V,LOMAX,Y,DYDT,DYDT2,DYDP,DYDP2,DYDTP,sth,cth) 
      implicit real*8 (a-h,o-z)
      COMPLEX*16 Y(*),DYDT(*),DYDT2(*),DYDP(*),DYDP2(*),DYDTP(*)   
      COMPLEX*16 IMAG                                               
      DIMENSION V(3)                                                 
      DATA ZERO/0.0D0/,SNULL/1.0E-10/,ONE/1.0D0/,FPI/12.5663706146D0/ 
      DATA IMAG/(0.D0,1.D0)/                                            
!                                                                    
!     CALCULATE SINES AND COSINES OF THE POLAR ANGLES OF THE VECTOR V
!                                                                   
      XY=V(1)**2+V(2)**2                                           
      IF(XY.lT.SNULL) then
        xy=1.1*snull   
      end if
      XYZ=XY+V(3)**2     
      IF(XY.GT.SNULL) GOTO 6   
      CTH=ONE                 
      IF(V(3).LT.ZERO) CTH=-CTH 
      STH=ZERO
!      sth=snull               
      CFI=ONE                 
      SFI=ZERO               
      GOTO 7                
    6 XY=DSQRT(XY)          
      XYZ= DSQRT(XYZ)      
      CTH=V(3)/XYZ       
      STH=XY/XYZ        
      CFI=V(1)/XY      
      SFI=V(2)/XY     
    7 CONTINUE       
                                                    
!     WRITE(6,*) 'CTH,STH,CFI,SFI',CTH,STH,CFI,SFI   
      DYDT(1)=(0.D0,0.D0)                              
      DYDT2(1)=(0.D0,0.D0)                            
      DYDP(1)=(0.D0,0.D0)                            
      DYDP2(1)=(0.D0,0.D0)                          
      DYDTP(1)=(0.D0,0.D0)                         
!                                             
      L=0         
  50  L=L+1      
      L1=L+1    
      LM1=L-1  
      DO 60 M=-L,L      
      IND=L*L1+M+1     
      DYDP(IND)=Y(IND)*IMAG*M     
      DYDP2(IND)=DYDP(IND)*IMAG*M 
      IF(IABS(M).EQ.L) THEN      
         DYDT(IND)=L*CTH/STH*Y(IND)   
         DYDT2(IND)=-L/STH/STH*Y(IND)+L*CTH/STH*DYDT(IND)  
      ELSE                                                
         IND1=LM1*L+M+1                                  
        if(ind1.lt.1) write(*,*) 'ind1'
         FACT=(2*L+1)*(L-M)*(L+M)/DFLOAT(2*L-1)          
         SQRTL=DSQRT(FACT)              
         DYDT(IND)=L*CTH/STH*Y(IND)-SQRTL/STH*Y(IND1)  
         DYDT2(IND)=-L/STH/STH*Y(IND)+L*CTH/STH*DYDT(IND)    &
            +CTH/STH/STH*SQRTL*Y(IND1)-SQRTL/STH*DYDT(IND1) 
      END IF                                               
         DYDTP(IND)=DYDT(IND)*IMAG*M                      
!     WRITE(*,100) L,M,IND,IND1,Y(IND),DYDT(IND),DYDT2(IND)   
!    *,DYDP(IND),DYDP2(IND),DYDTP(IND)                       
 100  FORMAT(4I3,6(2F8.3,2X))                               
  60  CONTINUE                                             
      IF(L.EQ.LOMAX) RETURN                               
      GOTO 50                                            
      END                                               
