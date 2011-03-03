      SUBROUTINE SETFFx(nkk,ifft1,ifft2,ifft3,rhok,fft)
!
!.....sets rhok into fft-field (like setfft, but just one spin
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
!      COMMON /GENER/   BR1(3,3),br2(3,3),avec(3,3)                 
      CHARACTER*67  ERRMSG
      COMPLEX*16 FFT(IFFT1,IFFT2,IFFT3),rhok(NWAV),TAUP(NSYM)          
      dimension kkk(3,NSYM)
      DO 9 I1=1,IFFT1
      DO 9 I2=1,IFFT2
      DO 9 I3=1,IFFT3                                                   
  9   FFT(I1,I2,I3)=0.0d0                                           
      DO 2 I=1,NKK                                                      
      CALL STERN(I,NST,KKK,TAUP)                                        
      DO 3 J=1,NST                                                      
      ii1=kkk(1,j)
      ii2=kkk(2,j)
      ii3=kkk(3,j)
      IF((IABS(ii1+1)*2).GT.IFFT1) GOTO 900                       
      IF((IABS(ii2+1)*2).GT.IFFT2) GOTO 900                       
      IF((IABS(ii3+1)*2).GT.IFFT3) GOTO 900                       
      I1=ii1                                                       
        IF(I1.LT.0) I1=I1+IFFT1                                         
      I2=ii2                                                       
        IF(I2.LT.0) I2=I2+IFFT2                                         
      I3=ii3                                                       
        IF(I3.LT.0) I3=I3+IFFT3                                         
!      IF(I1.EQ.0.AND.I2.EQ.0) WRITE(*,*) I3,KKK(3,J),RHOK(I)           
      FFT(I1+1,I2+1,I3+1)=RHOK(I)*TAUP(J)/nst                      
    3 CONTINUE                                                          
    2 CONTINUE                                                          
      return
!
!        Error messages
!
  900 CALL OUTERR('SETFFX','ifft too small in setffx')
      WRITE (ERRMSG,9000)
      CALL OUTERR('SETFFX',ERRMSG)
      WRITE (ERRMSG,9010) (KKK(II,J),II=1,3)
      CALL OUTERR('SETFFX',ERRMSG)
      WRITE (ERRMSG,9020) IFFT1,IFFT2,IFFT3
      CALL OUTERR('SETFFX',ERRMSG)
      STOP 'setffx - Error'
!
!        End of 'SETFFX'
!
 9000 FORMAT('2*(KKK+1) LARGER THAN IFFT PARAMETERS IN setffx')
 9010 FORMAT('KKK=',3I4)                                       
 9020 FORMAT('IFFT=',3I4)                                       
      end
