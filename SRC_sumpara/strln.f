      SUBROUTINE STRLN (INPUT,LENGTH)
      CHARACTER*(*) INPUT
      INTEGER LENGTH, LCOUNT
!     Evaluate the length of a string, skipping blanks
!     Note: this does not trap control characters
      DO 9998 LCOUNT=LEN(INPUT),1,-1
          IF (INPUT(LCOUNT:LCOUNT).NE.' ') GOTO 9999
 9998 CONTINUE
 9999 LENGTH=LCOUNT
      RETURN
      END
