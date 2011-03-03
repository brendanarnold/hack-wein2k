      SUBROUTINE PWXAD4(ll1,ll2,ll3,tvec,am)
!
!
      use struct
      IMPLICIT REAL*8(A-H,O-Z)
!
      DIMENSION TVEC(10,27),RVEC(3,27),LL(3),DELT(3,3),AM(125)
!
      ll(1)=ll1
      ll(2)=ll2
      ll(3)=ll3
      DO 199 I=1,3
      DO 199 J=1,3
  199 DELT(I,J)=AVEC(I,J)/LL(J)
      DO 200 J=1,125
  200 AM(J)=0.D0
      IND=0
      DO 201 I=-1,1,1
      DO 201 J=-1,1,1
      DO 201 K=-1,1,1
      IND=IND+1
      RVEC(1,IND)=DELT(1,1)*I+DELT(1,2)*J+DELT(1,3)*K
      RVEC(2,IND)=DELT(2,1)*I+DELT(2,2)*J+DELT(2,3)*K
      RVEC(3,IND)=DELT(3,1)*I+DELT(3,2)*J+DELT(3,3)*K
  201 CONTINUE
      DO 202 I=1,27
      TVEC(1,I)=1.D0
      TVEC(2,I)=RVEC(1,I)
      TVEC(3,I)=RVEC(2,I)
      TVEC(4,I)=RVEC(3,I)
      TVEC(5,I)=RVEC(1,I)**2
      TVEC(6,I)=RVEC(2,I)**2
      TVEC(7,I)=RVEC(3,I)**2
      TVEC(8,I)=RVEC(2,I)*RVEC(3,I)
      TVEC(9,I)=RVEC(1,I)*RVEC(3,I)
      TVEC(10,I)=RVEC(1,I)*RVEC(2,I)
      JJX=0
      DO 202 J=1,10
      DO 202 K=J,10
      JJX=10*(K-1)+J
      AM(JJX)=AM(JJX)+TVEC(J,I)*TVEC(K,I)
  202 CONTINUE
      CALL CHFAC(AM,10)
!
      return
      end
