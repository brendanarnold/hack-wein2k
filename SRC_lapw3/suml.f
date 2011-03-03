      SUBROUTINE SUML(LMMAX,TCC,YL,F,FA,LM,FFX)                         
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 YL(*),IMAG,IMAG1                                      
      DIMENSION TCC(*),LM(2,*),FFX(*)                          
        REAL*8 MINU                                                     
      DATA IMAG /(0.0d0,1.0d0)/                                             
      rt2i=1.d0/sqrt(2.d0)
      F=0.0d0                                                             
      FA=0.0d0                                                            
      DO 10 ILM=1,LMMAX                                                 
      IMAG1=(1.d0,0.d0)                                                     
      LL=IABS(LM(1,ILM))                                                
      MM=LM(2,ILM)                                                      
      MINU=1.d0                                                           
      IF(LM(1,ILM).LT.0) THEN                                           
      IMAG1=-IMAG                                                       
      MINU=-1.d0                                                          
      END IF                                                            
      IF(MOD(MM,2).EQ.1) THEN                                           
      IMAG1=-IMAG1                                                      
      MINU=-MINU                                                        
      END IF                                                            
      LY1=LL*(LL+1)+MM+1                                                
      IF(MM.EQ.0) GOTO 11                                               
      LY2=LL*(LL+1)-MM+1                                                
      IF(MOD(LL,2).NE.0) GOTO 12                                        
! ...SYMMETRIC PART, MM NE.0                                            
!      FFX(ILM)=TCC(ILM)*IMAG1*(YL(LY1)+YL(LY2)*MINU)/sqrt(2.d0)          
       FFX(ILM)=TCC(ILM)*IMAG1*(YL(LY1)+YL(LY2)*MINU)*rt2i
       F=F+FFX(ILM)                                                     
      GOTO 10                                                           
  11  IF(MOD(LL,2).NE.0) GOTO 13                                        
! ...SYMMETRIC PART,  MM EQ.0                                           
      FFX(ILM)=TCC(ILM)*IMAG1*YL(LY1)                                   
      F=F+FFX(ILM)                                                      
      GOTO 10                                                           
! ...ANTISYMMETRIC PART,  MM NE.0                                       
!  12  FFX(ILM)=TCC(ILM)*IMAG1*(YL(LY1)+YL(LY2)*MINU)/sqrt(2.d0)           
  12  FFX(ILM)=TCC(ILM)*IMAG1*(YL(LY1)+YL(LY2)*MINU)*rt2i
      FA=FA+FFX(ILM)                                                    
      GOTO 10                                                           
! ...ANTISYMMETRIC PART, MM EQ.0                                        
  13  FFX(ILM)=TCC(ILM)*IMAG1*YL(LY1)                                   
      FA=FA+FFX(ILM)                                                    
!                                                                       
  10  CONTINUE                                                          
      RETURN                                                            
      END                                                               
