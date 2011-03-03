      SUBROUTINE cub_xc_back(xc,ir,llmm,lm,bb,lmmax)
      IMPLICIT NONE
      INCLUDE 'param.inc'
      REAL*8   ::  c_kub(0:10,0:10)
      COMMON /norm_kub/ c_kub
      REAL*8,INTENT(OUT)    :: xc(nrad,ncom)
      REAL*8,INTENT(IN)     :: bb(ncom) 
      INTEGER,INTENT(IN)    :: llmm,lm(2,ncom+3),ir,lmmax
      REAL*8     :: sqrt2
      INTEGER    :: i,iyy


      sqrt2=SQRT(2.d0)

      iyy=0
      i=1
      DO
         IF(i.gt.lmmax) EXIT
         IF(lm(1,i).EQ.0.AND.lm(2,i).EQ.0) THEN
            iyy=iyy+1
            xc(ir,i)=bb(iyy)
            i=i+1
         ELSEIF (lm(1,i).EQ.-3.AND.lm(2,i).EQ.2) THEN  
            iyy=iyy+1
            xc(ir,i)=bb(iyy)
            i=i+1
         ELSEIF (lm(1,i).EQ.4.OR.lm(1,i).EQ.6.OR. &
                 lm(1,i).EQ.-7.OR.lm(1,i).EQ.-9) THEN  
            iyy=iyy+1
            xc(ir,i)=bb(iyy)/c_kub(abs(lm(1,i)),lm(2,i))
            xc(ir,i+1)=0.0d0
            i=i+2
         ELSEIF (lm(1,i).EQ.8.OR.lm(1,i).EQ.10) THEN 
            iyy=iyy+1
            xc(ir,i)=bb(iyy)/c_kub(lm(1,i),lm(2,i))
            xc(ir,i+1)=0.0d0
            xc(ir,i+2)=0.0d0
            i=i+3
         ELSE
            WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
            STOP
         ENDIF
      END DO
      IF(iyy.NE.llmm) THEN
         write(6,*) 'cub_xc_back:  iyy=',iyy, &
                         ' different from llmm=',llmm
         STOP
      ENDIF
      END SUBROUTINE cub_xc_back
