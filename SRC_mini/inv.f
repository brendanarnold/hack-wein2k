      SUBROUTINE INV(A,B,D)
      implicit double precision( a-h, o-z )
      DIMENSION A(3,3),B(3,3)
      DO  10  I = 1,3
      I1 = I + 1
      IF(I1.GT.3)  I1 = I1 - 3
      I2 = I + 2
      IF(I2.GT.3)  I2 = I2 - 3
      DO  10  J = 1,3
      J1 = J + 1
      IF(J1.GT.3)  J1 = J1 - 3
      J2 = J + 2
      IF(J2.GT.3)  J2 = J2 - 3
   10 B(I,J) = A(J1,I1)*A(J2,I2) - A(J1,I2)*A(J2,I1)
      D = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
      IF(ABS(D).LT.1.d-6) WRITE(6,100) &
       A(1,1),B(1,1),A(1,2),B(2,1),A(1,3),B(3,1)
100   FORMAT(' DET IST NULL,DENN A(1,1) = ',F15.5/ &
             '                   B(1,1) = ',F15.5/ &
             '                   A(1,2) = ',F15.5/ &
             '                   B(2,1) = ',F15.5/ &
             '                   A(1,3) = ',F15.5/ &
             '                   B(3,1) = ',F15.5)
      DO  20  I = 1,3
      DO  20  J = 1,3
      IF(ABS(D).LT.1.d-6) WRITE(6,110) I,J,B(I,J)
110   FORMAT(' ZU TEILEN WAERE B(',I1,',',I1,') = ',F15.5)
   20 B(I,J) = B(I,J)/D
      RETURN
      END
