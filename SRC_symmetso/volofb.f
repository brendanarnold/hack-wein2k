


      DOUBLE PRECISION FUNCTION VOLOFB (B)

!  CALCULATES THE VOLUME OF THE UNIT CELL, POSIBLY WITH A MINUS SIGN

      DOUBLE PRECISION B(3,3)
      VOLOFB = ( B(2,1)*B(3,2) - B(3,1)*B(2,2) ) * B(1,3) + &
               ( B(3,1)*B(1,2) - B(1,1)*B(3,2) ) * B(2,3) + &
               ( B(1,1)*B(2,2) - B(2,1)*B(1,2) ) * B(3,3)
      END
