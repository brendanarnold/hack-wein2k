      SUBROUTINE CPUTIM(TIME)                                           
      DOUBLE PRECISION   TIME
!
%INCLUDE '/sys/ins/base.ins.ftn'
%INCLUDE '/sys/ins/proc1.ins.ftn'
!
      INTEGER            I, J, K
      INTEGER*2          CLOCK_VALUE(3), CLOCKH
      INTEGER*4          CLOCKL
      DOUBLE PRECISION   T
!
      EQUIVALENCE (CLOCKH,CLOCK_VALUE(1)),(CLOCKL,CLOCK_VALUE(2))
!
      CALL PROC1_$GET_CPUT(CLOCK_VALUE)
      T=CLOCKL
      IF (CLOCKL .LT. 0) T = 2.0D+0**32  + CLOCKL - 1
      K = 0   
      DO 1 I = 1, 10
         DO 2 J = I, 1, -1 
            K = K + 1
            IF (CLOCKH .GE. INT2(K) .AND. CLOCKH .LT. INT2(K+J)) THEN
               T = T + 2.0D+0**(31+J)
            ENDIF
   2     CONTINUE 
         IF (CLOCKH .LE. INT2(K)) GOTO 3
   1  CONTINUE 
   3  CONTINUE
      T = T*4.0D+0/1000000.0D+0
      TIME = T
!
      RETURN                                                            
!
!        End of 'CPUTIM'
!
      END 
