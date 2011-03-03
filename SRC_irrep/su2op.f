      SUBROUTINE SU2OP(JIJ,IORD,DZ2,DET,SZ,ANG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!
!       see for instance ref.[2]:
!....................................................................................
!
      COMPLEX*16   SZ(2,2),IMAG,CZERO
      DIMENSION    DIX(3,3),DZ2(3,3)
      DIMENSION    ANG(3),JIJ(NSYM,NSYM)
      DATA         IMAG/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/
!**********************************************************************
!
      IS=1
      IF(DET.LT.0.D0) IS=-1 
      DO 10 I=1,3
      DO 10 J=1,3
 10      DIX(I,J)=IS*DZ2(I,J)
!
!   ..determine eulers angles
      CALL EULANG(DIX,ANG,IS)
!
!   ..u(R) for SU(2)
      SIP=DSIN(ANG(1)/2.D0)
      COP=DCOS(ANG(1)/2.D0)
      EPP=(ANG(2)+ANG(3))/2.D0
      ENN=(ANG(2)-ANG(3))/2.D0
!
      SZ(1,1)= COP*CMPLX(DCOS( EPP),DSIN( EPP))
      SZ(1,2)= SIP*CMPLX(DCOS( ENN),DSIN( ENN))
      SZ(2,1)=-SIP*CMPLX(DCOS(-ENN),DSIN(-ENN))
      SZ(2,2)= COP*CMPLX(DCOS(-EPP),DSIN(-EPP))
!
      RETURN
      END



