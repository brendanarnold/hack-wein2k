      subroutine iffpar(ifft1,ifft2,ifft3,iff1,iff2,iff3,nkk)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*67 ERRMSG
      complex*16 taup                                                 
      DIMENSION        KKK(3,NSYM),TAUP(NSYM)                      
      IFF1=0
      IFF2=0
      IFF3=0
      DO 40 J=2,NKK
         CALL STERN(J,IND,KKK,TAUP)
         DO 42 JJ=1,IND
            IF(KKK(1,JJ).GT.IFF1) IFF1=KKK(1,JJ)
            IF(KKK(2,JJ).GT.IFF2) IFF2=KKK(2,JJ)
            IF(KKK(3,JJ).GT.IFF3) IFF3=KKK(3,JJ)
 42      CONTINUE
 40   CONTINUE
      IFF1=(IFF1+1)*2
      IFF2=(IFF2+1)*2
      IFF3=(IFF3+1)*2
      IF(IFF1*IFF2*IFF3.GT.IFFT1*IFFT2*IFFT3) GOTO 900
      JFP=(2*IFF1+1)*(2*IFF2+1)*(2*IFF3+1)
      JFTP=(2*IFFT1+1)*(2*IFFT2+1)*(2*IFFT3+1)
      IF(JFP.GT.JFTP) GOTO 900
      WRITE(6,*) 'IFFT-PARAMETERS SET AND USED:',IFFT1,IFFT2,IFFT3, &
       IFF1,IFF2,IFF3 
      RETURN
!
!        Error messages
!
  900 CALL OUTERR('IFFPAR','ifft parameters too small')
      WRITE (ERRMSG,9000) IFF1,IFF2,IFF3
      CALL OUTERR('IFFPAR',ERRMSG)
      WRITE (ERRMSG,9010) IFFT1,IFFT2,IFFT3
      CALL OUTERR('IFFPAR',ERRMSG)
      STOP 'IFFPAR - Error'
!
!        End of 'IFFPAR'
!
 9000 FORMAT('NEEDED VALUES:',3I3)
 9010 FORMAT('PRESENT PARAMETERS:',3I3)
      END
