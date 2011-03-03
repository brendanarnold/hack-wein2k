      SUBROUTINE CPUTIM(TIME)
      real*8 time                                                  
      time=second()
      end
      function second(dummy)
      real tarray(2)
      second=etime(tarray)
      end
