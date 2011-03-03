      SUBROUTINE CRWCND(ISK,ISKDEN,IZ,TAU,LKG,IKG,FGT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!
      LOGICAL         FGT
      DIMENSION       ISK(3),IZ(3,3,NSYM),TAU(3,NSYM)
      DIMENSION       LKG(NSYM),RT(3)
      DATA            TOLPH/1.0E-12/
!**********************************************************************
      FGT=.TRUE.
      DO 100 I=1,IKG
      DO 100 J=I,IKG 
         DO 120 II=1,3
 120       RT(II)= DBLE(IZ(II,1,LKG(I)))*TAU(1,LKG(J))+ &
                   DBLE(IZ(II,2,LKG(I)))*TAU(2,LKG(J))+ &
                   DBLE(IZ(II,3,LKG(I)))*TAU(3,LKG(J))-TAU(II,LKG(J))
         AG=0.D0
         DO 140 II=1,3
 140        AG=AG - 2*PI*ISK(II)*RT(II)/DBLE(ISKDEN)
         DFF=(DCOS(AG)-1.D0)**2+(DSIN(AG))**2
         IF(DFF.GT.TOLPH) FGT=.FALSE.
 100  CONTINUE
!
!.....output
      IF(.NOT.FGT) WRITE(6,510)
!
      RETURN
 510  FORMAT( &
        /,7X,'Non-symmorphic crystal and k-point at the BZ surface:', &
        /,7X,'IR of the space group for this k-point cannot simple ', &
        /,7X,'be expressed as IR of the corresponding point group  ',  &
        /,7X,'times a phase factor, since exp(-i*k(Ri*tj-tj)) not 1', &
        /,7X,'for all pair of {Ri|ti} and {Rj|tj}.',//)
      END


