      subroutine srolyl(clmspu,clmspd,yy,c,sqfp, &
                        iatnr,llmm,fu,fd,ft)
!
      implicit real*8 (a-h,o-z)
!
      INCLUDE 'param.inc'
!
      dimension clmspu(ncom),clmspd(ncom),yy(ncom),c(7)
!
      fu=CLMSPU(1)*YY(1)/SQFP                  
      fd=CLMSPD(1)*YY(1)/SQFP                  
!-------------
      IF(IATNR.LT.0) THEN                                        
       DO 8 LM1=2,LLMM                                                
        fu=fu+CLMSPU(lm1)*YY(lm1)                  
        fd=fd+CLMSPD(lm1)*YY(lm1)                  
 8     continue
!       IF(fu.LT.0.0d0) then
!        WRITE(*,2929) LM1,fu,CLMSPU,YY  
!        fu=1.e-6
!        end if                                      
!        IF(fd.LT.0.0d0) then
!        WRITE(*,2929) LM1,fd,CLMSPD,YY  
!        fd=1.e-6
!        end if                                      
 2929 FORMAT(' LM1,fd,CLM,YY ',i4,(4E16.5))                  
!----------------------
      ELSE IF(LLMM.GT.1)  THEN                                       
       fu=fu+(C(2)*CLMSPU(2)+C(3)*CLMSPU(3))*YY(2)+ &
             (C(4)*CLMSPU(4)+C(5)*CLMSPU(5))*YY(3)  
       fd=fd+(C(2)*CLMSPD(2)+C(3)*CLMSPD(3))*YY(2)+ &
             (C(4)*CLMSPD(4)+C(5)*CLMSPD(5))*YY(3)  
       if(llmm.eq.4) then
       fu=fu+CLMSPU(6)*YY(4)
       fd=fd+CLMSPD(6)*YY(4)
       endif 
!--------------------------------         
      ENDIF                                                          
      ft=fu+fd 
!
      return
      end
