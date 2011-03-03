      SUBROUTINE SUML (LMMAX,TCC,YL,LM)                                 
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16 YL,IMAG,IMAG1                                          
      REAL*8  MINU                                                      
      DIMENSION YL(*),LM(2,NCOM+3)                                     
      DATA IMAG /(0.d0,1.0D0)/                                           
!                                                                       
      ILM=LMMAX                                                         
      IMAG1=(1.D0,0.d0)                                                   
      LL=IABS(LM(1,ILM))                                                
      MM=LM(2,ILM)                                                      
      MINU=1.                                                           
      IF(LM(1,ILM).LT.0) THEN                                           
      IMAG1=-IMAG                                                       
      MINU=-1.D0                                                        
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
      TCC=IMAG1*(YL(LY1)+YL(LY2)*MINU)/SQRT(2.D0)                       
      GOTO 10                                                           
  11  IF(MOD(LL,2).NE.0) GOTO 13                                        
! ...SYMMETRIC PART,  MM EQ.0                                           
      TCC=IMAG1*YL(LY1)                                                 
      GOTO 10                                                           
! ...ANTISYMMETRIC PART,  MM NE.0                                       
  12  TCC=IMAG1*(YL(LY1)+YL(LY2)*MINU)/SQRT(2.D0)                       
!-------------------------------------------------------                
! 12  TCC=IMAG*( YL(LY1)+YL(LY2)*MINU ) / SQRT(2.)                      
!-------------------------------------------------------                
      GOTO 10                                                           
! ...ANTISYMMETRIC PART, MM EQ.0                                        
  13  TCC=IMAG1*YL(LY1)                                                 
!                                                                       
  10  CONTINUE                                                          
      RETURN                                                            
      END                                                               
