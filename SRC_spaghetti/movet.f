      subroutine movet(xwert,ywert)
      IMPLICIT REAL*8 (A-H,O-Z)
! '/M { CM exch CM exch moveto} def',/
      write(11,100) xwert,ywert
      return
 100  format(f10.5,' ',f10.5,' M')
      end

