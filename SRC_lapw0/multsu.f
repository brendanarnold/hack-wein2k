SUBROUTINE multsu(qq,fc,jatom,lmmax,iatnr,lm)
!                                                                       
!     MULTSU COMBINES THE MULTIPOLMOMENTS QQ                            
!
  USE struct, ONLY: nat
  IMPLICIT NONE
!  IMPLICIT REAL*8 (A-H,O-Z)
!
  INCLUDE 'param.inc'
  !
  COMPLEX*16   :: qq(ncom+3,*),fc(ncom+3,*),imag,a,b,c
  INTEGER      :: lm(2,ncom+3,nat),jatom,lmmax,iatnr,llmm,i
  REAL*8       :: sqrt2,sq1,c1,c2,c3
  REAL*8       :: c_kub(0:10,0:10)
  COMMON /norm_kub/ c_kub
!----------------------------------------------------------------------    
!
  c_kub=0.0d0
  c_kub(0,0)=1.d0
  c_kub(3,2)=1.d0
  c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
  c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
  c_kub(6,0)=.5d0*SQRT(.5d0)
  c_kub(6,2)=.25d0*SQRT(11.d0)
  c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
  c_kub(6,6)=-.25d0*SQRT(5.d0)
  c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
  c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
  c_kub(8,0)=.125d0*SQRT(33.d0)
  c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
  c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
  c_kub(9,2)=.25d0*SQRT(3.d0)
  c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
  c_kub(9,6)=-.25d0*SQRT(13.d0)
  c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
  c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
  c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
  c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
  c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
  c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
  c_kub(10,10)=-.0625d0*SQRT(85.d0)

  sqrt2=SQRT(2.d0)
  !                                                                       
  IF(IATNR.LT.0) THEN                                               
     !        NONCUBIC CASE                                                  
     DO LLMM=1,LMMAX                                              
        QQ(LLMM,JATOM)=QQ(LLMM,JATOM)*fc(llmm,jatom)         
     ENDDO
  ELSE                                                              
     !        CUBIC CASE
        
     i=1
     DO
        IMAG=CMPLX(1.D0,0.D0)
        IF(lm(1,i,jatom).LT.0) IMAG=CMPLX(0.D0,-1.D0)
        
        IF(i.gt.lmmax) EXIT
        IF(lm(1,i,jatom).EQ.0.AND.lm(2,i,jatom).EQ.0) THEN
           i=i+1
        ELSEIF (lm(1,i,jatom).EQ.-3.AND.lm(2,i,jatom).EQ.2) THEN  
           qq(i,jatom)=qq(i,jatom)*imag/sqrt2
           i=i+1
        ELSEIF (lm(1,i,jatom).EQ.4.OR.lm(1,i,jatom).EQ.6.OR. &
             lm(1,i,jatom).EQ.-7.OR.lm(1,i,jatom).EQ.-9) THEN
           IF (lm(2,i,jatom).EQ.0) THEN
              sq1=1.d0
           ELSE
              sq1=sqrt2
           ENDIF
           c1=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom))
           c2=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom)+4)
           a=qq(i,jatom)*imag
           b=qq(i+1,jatom)*imag
           qq(i,jatom)=a*c1*c1/sq1 + b*c1*c2/sq1
           qq(i+1,jatom)=a*c1*c2/sqrt2 + b*c2*c2/sqrt2        
           i=i+2
        ELSEIF (lm(1,i,jatom).EQ.8.OR.lm(1,i,jatom).EQ.10) THEN
           IF (lm(2,i,jatom).EQ.0) THEN
              sq1=1.d0
           ELSE
              sq1=sqrt2
           ENDIF
           a=qq(i,jatom)
           b=qq(i+1,jatom)
           c=qq(i+2,jatom)
           c1=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom))
           c2=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom)+4)
           c3=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom)+8)
           qq(i,jatom)=a*c1*c1/sq1 + b*c1*c2/sq1 + &
                c*c1*c3/sq1
           qq(i+1,jatom)=a*c1*c2/sqrt2 + b*c2*c2/sqrt2  &
                + c*c2*c3/sqrt2
           qq(i+2,jatom)=a*c1*c3/sqrt2 + b*c2*c3/sqrt2  &
                + c*c3*c3/sqrt2
           i=i+3
        ELSE
           WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
           WRITE(6,*) 'MULTSU.F'
           WRITE(6,'(a2,i4,a3,i4)') 'L=',lm(1,i,jatom),' M=', &
                lm(2,i,jatom)
           STOP
        ENDIF
     END DO
  ENDIF
  RETURN                                                            
END SUBROUTINE MULTsu
