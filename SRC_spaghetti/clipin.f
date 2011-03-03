      subroutine clipin(xa,ya,xe,ye)
      IMPLICIT REAL*8 (A-H,O-Z)
      write(11,*) 'initclip'
      write(11,*) '/windowpath'
      write(11,*) '{newpath'
      write(11,100) xa,ya,xe-xa,ye-ya,xa-xe,ya-ye
      write(6,100) xa,ya,xe-xa,ye-ya,xa-xe,ya-ye
      write(11,*) 'closepath} def'
      write(11,*) 'windowpath clip'
      write(11,*) 'newpath'
      return
 100  FORMAT(2f10.5,' M ',f10.5,' 0.0 RL 0.0 ',f10.5,' RL ' &
             ,f10.5,' 0.0 RL 0.0 ',f10.5)
      end

