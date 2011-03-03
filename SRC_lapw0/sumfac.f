      SUBROUTINE SUMFAC(LLMM,YY,IND,LM,LMMAX)
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      REAL*8   ::  c_kub(0:10,0:10)
      COMMON /norm_kub/ c_kub
      INTEGER     :: LM(2,NCOM+3)
      DIMENSION YY(nspd,NCOM),C(7)                                       

!      sqrt2=SQRT(2.d0)

      llmm=0
      i=1
      DO
         IF(i.gt.lmmax) EXIT
         IF(lm(1,i).EQ.0.AND.lm(2,i).EQ.0) THEN
            llmm=llmm+1
            yy(ind,llmm)=yy(ind,i)
            i=i+1
         ELSEIF (lm(1,i).EQ.-3.AND.lm(2,i).EQ.2) THEN  
            llmm=llmm+1
            yy(ind,llmm)=yy(ind,i)
            i=i+1
         ELSEIF (lm(1,i).EQ.4.OR.lm(1,i).EQ.6.OR. &
                 lm(1,i).EQ.-7.OR.lm(1,i).EQ.-9 ) THEN  
            llmm=llmm+1
            yy(ind,llmm)=c_kub(abs(lm(1,i)),lm(2,i))*yy(ind,i)  &
                       + c_kub(abs(lm(1,i)),lm(2,i)+4)*yy(ind,i+1)
            i=i+2
         ELSEIF (lm(1,i).EQ.8.OR.lm(1,i).EQ.10) THEN 
            llmm=llmm+1
            yy(ind,llmm)=c_kub(lm(1,i),lm(2,i))*yy(ind,i) + &
                          c_kub(lm(1,i),lm(2,i)+4)*yy(ind,i+1) +  &
                          c_kub(lm(1,i),lm(2,i)+8)*yy(ind,i+2) 
            i=i+3
         ELSE
            WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
            STOP
         ENDIF
      END DO
      RETURN                                                            
      END                                                               
