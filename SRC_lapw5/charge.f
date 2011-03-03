      SUBROUTINE CHARGE(CHG,IR,R,JATOM,V,IATNR,clm,lm,lmmx,nat,rm,jri)    
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      dimension RM(NRAD,NAT),jri(nat)
      DIMENSION RHO(NCOM),ANG(NCOM),V(3)                                
      COMPLEX*16 YL((lmax2+1)*(lmax2+1)),IMAG,IMAG1                       
      REAL*8 MINU                                                        
      dimension  CLM(NRAD,NCOM,NAT),LM(2,NCOM,NAT),LMMX(NAT)           
!      DATA IMAG/(0.0,1.0)/,SQIN2/0.7071067812/                          
      parameter (imag=(0.d0,1.d0), sqin2=0.707106781186547D0)
      CALL YLM(V,lmax2,YL)                                                  
      LMMAX=LMMX(JATOM)                                                 
      DO 1 ILM=1,LMMAX                                                  
      CALL RADIAL(CLM(1,ILM,JATOM),RHO(ILM),R,IR,JATOM,rm,jri,nat)
      L=IABS(LM(1,ILM,JATOM))                                           
      M=LM(2,ILM,JATOM)                                                 
      MINU=1.d0                                                           
      IMAG1=(1.d0,0.d0)                                                     
      IF(LM(1,ILM,JATOM).LT.0) THEN                                     
      IMAG1=-IMAG                                                       
      MINU=-1.d0                                                          
      END IF                                                            
      IF(MOD(M,2).EQ.1) THEN                                            
      IMAG1=-IMAG1                                                      
      MINU=-MINU                                                        
      END IF                                                            
      IF(M.NE.0) GOTO 2                                                 
      IDX=L*(L+1)+1                                                     
      ANG(ILM)=YL(IDX)                                                  
      GOTO 1                                                            
    2 IDM=L*(L+1)+1                                                     
      IDP=IDM+M                                                         
      IDM=IDM-M                                                         
      ANG(ILM)=(YL(IDP)+MINU*YL(IDM))*SQIN2*IMAG1                       
   1  CONTINUE                                                          
!                                                                       
      IF(IATNR.GT.0) CALL SUM(RHO,ANG,CHG,LMMAX,lm,jatom)
      IF(IATNR.GT.0) RETURN                                             
!                                                                       
      CHG=0.0d0                                                           
      DO 10 ILM=1,LMMAX 
   10 CHG=CHG+RHO(ILM)*ANG(ILM)                                         
!CCC       if(jatom.eq.2) write(6,*) chg,(rho(ilm)*ang(ilm),ilm=1,lmmax)
      RETURN                                                            
      END                                                               
