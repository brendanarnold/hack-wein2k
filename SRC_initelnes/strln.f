

      SUBROUTINE STRLN (INPUT,LENGTH)
      CHARACTER *80 INPUT
      INTEGER LENGTH
      INTEGER LCOUNT
!     Evaluate the length of the string
!ph   Elimine the ' ' characters from the line end.
      DO 9998 LCOUNT=LEN(INPUT),1,-1
         IF (INPUT(LCOUNT-1:LCOUNT).NE.' ') GOTO 9999
 9998 CONTINUE
 9999 LENGTH=LCOUNT-1
!     Got the length of the string in LENGTH
      END
