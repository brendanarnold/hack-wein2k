      SUBROUTINE ROTKV(LKG,IKG,NV,ISK,ISKDEN,IZ,IIZ,TAU) 
      use felder
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      INCLUDE 'param.inc'
!
!      COMPLEX*16      PH(NSYM,NMAT)
      DIMENSION       IZ(3,3,NSYM),IIZ(3,3,NSYM),TAU(3,NSYM)
      DIMENSION       ISK(3)
      DIMENSION       KVR(3),LKG(NSYM)
!                     ,KV(3,NMAT),L(NSYM,NMAT)
!**********************************************************************
!
!.....find K' such that (k+K)inv(Ri)=k+K'
      DO 10 IG=1,IKG
      DO 10 N=1,NV
        DO 30 J=1,3
 30       KVR(J)=KV(1,N)*IIZ(1,J,LKG(IG)) &
                +KV(2,N)*IIZ(2,J,LKG(IG)) &
                +KV(3,N)*IIZ(3,J,LKG(IG))
!
        IN=1
        DO WHILE(( (ABS(KVR(1)-KV(1,IN)) &
                   +ABS(KVR(2)-KV(2,IN)) &
                   +ABS(KVR(3)-KV(3,IN))).NE.0).AND.(IN.LE.NV))
          IN=IN+1
        END DO
!
        IF(IN.EQ.NV+1)then
           WRITE(6,*) 'rotkv: cannot find (k+K)inv(Ri)'
           WRITE(6,*) IG,N, NV,(KV(J,N),J=1,3)
           WRITE(6,*) IG,IN,NV,(KVR(J), J=1,3)
           WRITE(6,*) 'All k-vectors:'
           DO 888 II=1,NV
 888          WRITE(6,*)'k+K=',(kv(J,II),J=1,3)
           STOP 'rotkv: cannot find (k+K)inv(Ri)'
        endif
!
        L(IG,N)=IN
        AG=0.D0
        DO 40 J=1,3
 40       AG=AG-2*PI*DBLE(KV(J,IN)-ISK(J))*TAU(J,LKG(IG)) &
                    /DBLE(ISKDEN) 
          PH(IG,N)=CMPLX(DCOS(AG),DSIN(AG))
 10   CONTINUE
!
      RETURN
      END


