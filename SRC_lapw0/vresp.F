!@PROCESS DC(DENSIT)
      SUBROUTINE VRESP (NKK1,JSPIN,CWORK,DWORK, &
                        IFFT1,IFFT2,IFFT3,CFFT,CFFT1,FFT,VROKSP)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*19     nkktext
      INTEGER            IFFT1, IFFT2, IFFT3, JSPIN, NKK1
      DOUBLE PRECISION   DWORK(*)
      COMPLEX*16         CWORK(*)
      COMPLEX*16         CFFT(IFFT1,IFFT2,IFFT3)
      COMPLEX*16         CFFT1(IFFT1,IFFT2,IFFT3)
      COMPLEX*16         FFT(IFFT1,IFFT2,IFFT3,2)
#if defined (SCSL)
      INTEGER            ISYS(0:1)
#endif
!
!    VRESP is numerically evaluated in the unit cell.
!
      INTEGER            I1, I2, I3, IERR, II, ISIG, ISPIN
      INTEGER            J, J1, JX, NKK, NKKSP
      INTEGER            KZZSP(3)
      DOUBLE PRECISION   EXD, EXU, PI, SPINDN, SPINUP, TPI, VXD, VXU
      DOUBLE PRECISION   AM(125)
      DOUBLE PRECISION   TVEC(10,27)
!
      COMPLEX*16         VROKSP(NKK1,2)
!
!
      DOUBLE PRECISION   XCENER, XCPOT
      EXTERNAL           XCENER
      DATA THRD/0.3333333333333333D0/
!
#if defined (SCSL)
      ISYS(0) = 1
#endif

      IF(JSPIN.EQ.2) THEN
!
!.....read plane wave spin density VROKSP(K)
!
         READ(29,2060) NKKSP
!         READ(30,2060) NKKSP
  read(30,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkksp
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkksp
 6768 continue
         DO 4900 J=1,NKKSP
            READ(29,2071)  ( KZZSP(JX), JX=1,3 ) ,VROKSP(J,1)
            READ(30,2071)  ( KZZSP(JX), JX=1,3 ) ,VROKSP(J,2)
 4900    CONTINUE
      ELSE
         READ(28,2060) NKKSP
         DO 4903 J=1,NKKSP
            READ(28,2071)  ( KZZSP(JX), JX=1,3 ) ,VROKSP(J,1)
            VROKSP(J,1)=VROKSP(J,1)/2.D0
            VROKSP(J,2)=VROKSP(J,1)
 4903    CONTINUE
      END IF
 2060 FORMAT(/,13X,I6)
 2070 FORMAT(3X,3I5,2E15.7)
 2071 FORMAT(3X,3I5,2E19.12)
      PI=4.0D+0*ATAN(1.0D+0)
      TPI=2.D0*PI
      NKK=NKK1
!
!
!.....set the fft-field with rho-fourier-coefficients
!
      CALL SETFF1(NKK,IFFT1,IFFT2,IFFT3,VROKSP,FFT)
!
!.....transform rho to real space
!
      ISIG=0
#if defined (SCSL)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,0.0D0,0,1,1,0,1,1,CWORK,DWORK,&
                   ISYS)
#else
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT,IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN INITIALIZATION',IERR
#endif
      ISIG=-1
#if defined (SCSL)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT(1,1,1,1),IFFT1,     &
                   IFFT2,FFT(1,1,1,1),IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT(1,1,1,2),IFFT1,     &
                   IFFT2,FFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS)
#else
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,1),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
#endif
!
!.....calculate xcpotential and energy
!
      DO 8 I1=1,IFFT1
         DO 8 I2=1,IFFT2
            DO 8 I3=1,IFFT3
               CFFT(I1,I2,I3)=FFT(I1,I2,I3,1)
               CFFT1(I1,I2,I3)=FFT(I1,I2,I3,2)
 8    CONTINUE
      RETURN
      END
