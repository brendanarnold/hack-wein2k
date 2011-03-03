      subroutine srolyl(clmspu,clmspd,yy,sqfp, &
                        iatnr,llmm,fu,fd,ft,lm,lmmax)
!
      implicit real*8 (a-h,o-z)
!
      INCLUDE 'param.inc'
!
      REAL*8   ::  c_kub(0:10,0:10)
      COMMON /norm_kub/ c_kub
      INTEGER     :: LM(2,NCOM+3)
      dimension clmspu(ncom),clmspd(ncom),yy(ncom)
!
!      sqrt2=SQRT(2.d0)

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
! 2929 FORMAT(' LM1,fd,CLM,YY ',i4,(4E16.5))                  
!----------------------
!      ELSE IF(LLMM.GT.1)  THEN                                       

      ELSE

      iyy=0
      i=1
      DO
!         write(6,*) 'iyy',iyy

         IF(i.gt.lmmax) EXIT
         IF(lm(1,i).EQ.0.AND.lm(2,i).EQ.0) THEN
            iyy=iyy+1
            i=i+1
         ELSEIF (lm(1,i).EQ.-3.AND.lm(2,i).EQ.2) THEN  
            iyy=iyy+1
!            write(6,*) 'LM(-3,2). iyy=',iyy,'i=',i
            fu=fu+clmspu(i)*yy(iyy)
            fd=fd+clmspd(i)*yy(iyy)
            i=i+1
         ELSEIF (lm(1,i).EQ.4.OR.lm(1,i).EQ.6.OR. &
                 lm(1,i).EQ.-7.OR.lm(1,i).EQ.-9 ) THEN  
            c1=c_kub(abs(lm(1,i)),lm(2,i))
            c2=c_kub(abs(lm(1,i)),lm(2,i)+4)
            iyy=iyy+1
            fu=fu+(c1*CLMSPU(i)+c2*CLMSPU(i+1))*YY(iyy)
            fd=fd+(c1*CLMSPd(i)+c2*CLMSPd(i+1))*YY(iyy)
            i=i+2
         ELSEIF (lm(1,i).EQ.8.OR.lm(1,i).EQ.10) THEN 
            iyy=iyy+1
            c1=c_kub(abs(lm(1,i)),lm(2,i))
            c2=c_kub(abs(lm(1,i)),lm(2,i)+4)
            c3=c_kub(abs(lm(1,i)),lm(2,i)+8)
            fu=fu+(c1*CLMSPU(i)+c2*CLMSPU(i+1) +  &
                            c3*CLMSPU(i+2)  )*YY(iyy)
            fd=fd+(c1*CLMSPd(i)+c2*CLMSPd(i+1) +  &
                            c3*CLMSPd(i+2)  )*YY(iyy)
            i=i+3
         ELSE
            WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
            WRITE(6,*) 'srolyl.f'
            STOP
         ENDIF
      END DO
      IF(iyy.NE.llmm) THEN
         write(6,*) 'iyy=',iyy,' different from llmm=',llmm
         STOP
      ENDIF
      ENDIF

      ft=fu+fd 
!
      return
      end
