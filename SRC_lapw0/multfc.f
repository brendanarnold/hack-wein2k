      SUBROUTINE MULTFC(fc,JATOM,LMMAX,LM,IATNR)                        
!                                                                       
!     MULTFC COMBINES THE MULTIPOLMOMENTS QQ                            
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16   fc(NCOM+3,*),IMAG,IMAG1,A,B                          
      DIMENSION LM(2,*)                                                 
      DATA      IMAG/(0.D0,1.D0)/                                       
!-------------------------------------------------------------------    
!                                                                       
!                                                                       
!      IF(IATNR.LT.0) THEN                                               
      IF(.true.) THEN                                               
!        NONCUBIC CASE                                                  
         DO 1 LLMM=1,LMMAX                                              
            MM=LM(2,LLMM)                                               
            IF(MM.EQ.0) then 
               fc(llmm,jatom)=(1.d0,0.d0)
               GOTO 1                                          
            end if   
      IMAG1=(1.,0.)                                               
      MINU=1                                                            
            IF(LM(1,LLMM).LT.0) THEN                                    
         IMAG1=-IMAG                                                    
         MINU=-1                                                        
         END IF                                                         
            IF(MOD(MM,2).EQ.1)  THEN                                    
          IMAG1=-IMAG1                                                  
          MINU=-MINU                                                    
          END IF                                                        
      IF(MM.GT.0) MINU=1                                                
            fc(LLMM,JATOM)=IMAG1/SQRT(2.d0) * MINU         
!      IF(LM(1,LLMM).gt.0.and.MOD(MM,2).EQ.0.and.mm.gt.0) imag1=1.
!      IF(LM(1,LLMM).gt.0.and.MOD(MM,2).EQ.0.and.mm.lt.0) imag1=1.
!      IF(LM(1,LLMM).lt.0.and.MOD(MM,2).EQ.0.and.mm.gt.0) imag1=-imag
!      IF(LM(1,LLMM).lt.0.and.MOD(MM,2).EQ.0.and.mm.lt.0) imag1=imag
!      IF(LM(1,LLMM).gt.0.and.MOD(MM,2).EQ.1.and.mm.gt.0) imag1=-1.
!      IF(LM(1,LLMM).gt.0.and.MOD(MM,2).EQ.1.and.mm.lt.0) imag1=1.
!      IF(LM(1,LLMM).lt.0.and.MOD(MM,2).EQ.1.and.mm.gt.0) imag1=imag
!      IF(LM(1,LLMM).lt.0.and.MOD(MM,2).EQ.1.and.mm.lt.0) imag1=imag
!            fc(LLMM,JATOM)=IMAG1/SQRT(2.)          
 1       CONTINUE                                                       
      ELSE                                                               
!        CUBIC CASE LM= 00,40,44,60,64,-32                              
         DO llmm=1,ncom+3
            fc(llmm,jatom)=(1.d0,0.d0)
            IF(ABS(lm(1,llmm)).EQ.3.AND.lm(2,llmm).EQ.2)  &
                         fc(llmm,jatom)=(-IMAG)/SQRT(2.D0)
            IF(ABS(lm(1,llmm)).EQ.3.AND.lm(2,llmm).EQ.-2)  &
                         fc(llmm,jatom)=(IMAG)/SQRT(2.D0)
         ENDDO
      ENDIF
      RETURN                                                            
      END                                                               
