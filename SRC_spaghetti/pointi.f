      subroutine pointi(kcol,sizek,iprto)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer iprto(3)
!      if (kcol.eq.1) then
!         write(11,*) '1.0 setgray'
!      else
!         write(11,*) '0.0 setgray'
!      endif
      if(iprto(1).ne.3) then
        if (kcol.eq.0) then
           write (11,100) sizek
        endif
        write (11,110) sizek
      else
        write (11,100) sizek
        write(11,*) '0.0 setgray'
        write (11,110) sizek
      endif

      return
 100  format( '0 0 ',f10.5,' CM 0 360 arc fill')
 110  format( '0 0 ',f10.5,' CM 0 360 arc stroke')
      end

