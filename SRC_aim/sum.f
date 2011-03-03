       SUBROUTINE SUM(rho,ang,chg,lmmax,lm,jatom)
! Calculates cubic harmonics after Kara & Kurki-Suonio
! Acta Cryst A 1981 37 201-210
! GKHM 2/5-01
       IMPLICIT NONE
       INCLUDE 'param.inc'
       INTEGER  jatom,lmmax,lm(2,ncom,*)
       REAL*8   rho(ncom),ang(ncom)
       REAL*8   chg,c1,c2,c3
       REAL*8   c_kub(0:10,0:10)
       INTEGER  i,j
!       WRITE(6,*) 'WARNING. POS.ATOM NUMBER'
!       WRITE(6,*) 'THE LMs MUST COME IN THE ORDER GIVEN IN sum.f'

       chg=0.0                                                           
      do i=0,10
       do j=0,10
         c_kub(i,j)=0.0d0
       enddo
      enddo
      c_kub(0,0)=1.d0
      c_kub(3,2)=1.d0
      c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
      c_kub(4,4)=.5*SQRT(5.d0/3.d0)
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
      c_kub(10,0)=.125*SQRT(65.D0/6.D0)
      c_kub(10,2)=.125*SQRT(247.D0/6.D0)
      c_kub(10,4)=-.25*SQRT(11.D0/2.D0)
      c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
      c_kub(10,8)=-.125*SQRT(187.D0/6.D0)
      c_kub(10,10)=-.0625d0*SQRT(85.d0)


       i=1
 1     CONTINUE
          IF(i.gt.lmmax) GOTO 4
          IF(lm(1,i,jatom).EQ.0.AND.lm(2,i,jatom).EQ.0) THEN
             chg=chg + (rho(i) * ang(i))
             i=i+1
          ELSEIF (lm(1,i,jatom).EQ.-3.AND.lm(2,i,jatom).EQ.2) THEN  
             chg=chg+rho(i)*ang(i)
             i=i+1
          ELSEIF (lm(1,i,jatom).EQ.4.OR.lm(1,i,jatom).EQ.6.OR. &
                  lm(1,i,jatom).EQ.-7.OR.lm(1,i,jatom).EQ.-9) THEN  
            c1=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom))
            c2=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom)+4)
             chg=chg + (c1*ang(i) + c2*ang(i+1)) *  &
                 (c1*rho(i) + c2*rho(i+1))
             i=i+2
          ELSEIF (lm(1,i,jatom).EQ.8.OR.lm(1,i,jatom).EQ.10) THEN 
            c1=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom))
            c2=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom)+4)
            c3=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom)+8)
             chg=chg + (c1*ang(i) + c2*ang(i+1) + c3*ang(i+2)) &
                    * (c1*rho(i) + c2*rho(i+1) + c3*rho(i+2))
             i=i+3
          ELSE
             WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
             WRITE(6,'(a2,i4,a3,i4)') 'L=',lm(1,i,jatom), &
                   ' M=',lm(2,i,jatom)
             STOP
          ENDIF
       GOTO 1
 4     CONTINUE 
       END
