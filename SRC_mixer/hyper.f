      subroutine hyper(clmnew,r,hyperf,zz,dx,factor)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      dimension rho(NRAD),clmnew(NRAD),r(NRAD)
      pi=acos(-1.d0)
          do 4837 jrj=1,NRAD
          rho(jrj)=clmnew(jrj)*factor
!.....find Thomson radius
!.....average densities over R-Thomson (B35,3271)
          if(r(jrj).lt.ZZ*0.00005325d0) goto 4837
          if(jrj.lt.3) then
          write(6,*)' r0 .gt. r-thomson',zz,r(1),zz*0.00005325d0
          goto 4837
          end if
          CALL CHARGE(R,DX,rho,1,JRJ,QEL)
          hyperf=qel/(4.d0/3.d0*pi*(r(jrj)**3-r(1)**3))
          return
 4837     continue
      CALL OUTERR('HYPER','stop in hyper')
      STOP 'HYPER - Error'
      end
