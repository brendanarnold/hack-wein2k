      SUBROUTINE SETFFT(nkk,ifft1,ifft2,ifft3,rhoksp,fft)
!
!.....sets rhoksp into fft-field
!
      use struct
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      CHARACTER*67  ERRMSG
      COMPLEX*16 FFT(IFFT1,IFFT2,IFFT3,2),rhoksp(nkk,2),TAUP(NSYM)          
      dimension kkk(3,NSYM)

      FFT(1:IFFT1,1:IFFT2,1:IFFT3,1:2)=0.0d0

      DO 2 I=1,NKK                                                      
      CALL STERN(I,NST,KKK,TAUP)                                        
      DO 3 J=1,NST                                                      
!      IF(I.LT.6) WRITE(*,*) I,KKK(1,J),KKK(2,J),KKK(3,J),TAUP(J)       
!     transfer karth into internal coord.
!     if(ortho) then
!        ii1=nint((avec(1,1)*kkk(1,j)+avec(1,2)*kkk(2,j)+
!    *   avec(1,3)*kkk(3,j))/aa)
!        ii2=nint((avec(2,1)*kkk(1,j)+avec(2,2)*kkk(2,j)+
!    *   avec(2,3)*kkk(3,j))/bb)
!        ii3=nint((avec(3,1)*kkk(1,j)+avec(3,2)*kkk(2,j)+
!    *   avec(3,3)*kkk(3,j))/cc)
!      if(i.lt.5)write(6,*)i,j,(kkk(j1,j),j1=1,3),ii1,ii2,ii3
!     end if
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
      FFT(I1+1,I2+1,I3+1,1)=RHOKSP(I,1)*TAUP(J)/nst                      
      FFT(I1+1,I2+1,I3+1,2)=RHOKSP(I,2)*TAUP(J)/nst                      
    3 CONTINUE                                                          
    2 CONTINUE                                                          
      return
!
!        Error messages
!
  900 CALL OUTERR('SETFFT','ifft too small in xcpot3')
      WRITE (ERRMSG,9000)
      CALL OUTERR('SETFFT',ERRMSG)
      WRITE (ERRMSG,9010) (KKK(II,J),II=1,3)
      CALL OUTERR('SETFFT',ERRMSG)
      WRITE (ERRMSG,9020) IFFT1,IFFT2,IFFT3
      CALL OUTERR('SETFFT',ERRMSG)
      STOP 'XCPOT3 - Error'
!
!        End of 'SETFFT'
!
 9000 FORMAT('2*(KKK+1) LARGER THAN IFFT PARAMETERS IN XCPOT3')
 9010 FORMAT('KKK=',3I4)
 9020 FORMAT('IFFT=',3I4)
      end
