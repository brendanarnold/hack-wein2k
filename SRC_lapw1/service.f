	real function dtime4(lasttime)
	implicit none
	real lasttime(3)

!
!  Block packed algorithms                                           1997
!  Based on LAPACK routines                              Dieter Kvasnicka
!                                                   Theoretical Chemistry
!                                         University of Technology Vienna
!

	double precision thetime, tarray(2), thistime(2)

!	real etime
!	external etime
	double precision cputim
	external cputim
        
	thistime(1) = cputim(tarray)
	thistime(2) = tarray(2)
	thetime = thistime(1) - lasttime(1)
	lasttime(3) = thistime(2) - lasttime(2)
	lasttime(1) = thistime(1)
	lasttime(2) = thistime(2)

	dtime4 = thetime

	end
        subroutine printtime4(text,lasttime)
	implicit none

!
!  Block packed algorithms                                           1997
!  Based on LAPACK routines                              Dieter Kvasnicka
!                                                   Theoretical Chemistry
!                                         University of Technology Vienna
!

	  character text*(*)
          double precision lasttime(3)
	  double precision thistime
	  double precision dtime4

	  thistime=dtime4(lasttime)
          write (6,*) text, thistime, lasttime(3)
!     $	  thistime/lasttime(3)

        end 
