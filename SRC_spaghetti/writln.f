      subroutine writln(xwert,ywert,xwert2,ywert2)
      IMPLICIT REAL*8 (A-H,O-Z)
! '/M { CM exch CM exch moveto} def',/
      write(11,100) xwert, ywert
      write(11,101) xwert2,ywert2
      write(11,102)
      return
 100  format(f10.5,' ',f10.5,' M')
 101  format(f10.5,' ',f10.5,' L')
 102  format(' stroke ')
      end

