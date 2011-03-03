      SUBROUTINE REAN4 (NKK1,IFFT1,IFFT2,IFFT3,UST,FFT,VKVXC)
      use symetr
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      INTEGER            IFFT1, IFFT2, IFFT3, NKK1
#if defined (SCSL)
      DOUBLE PRECISION   DWORK(2 * MAX(IFFT1, IFFT2, IFFT3))
      COMPLEX*16         CWORK((2*IFFT1+256)+(2*IFFT2+256)+(2*IFFT3+256))
      INTEGER            ISYS(0:1)
#else
      DOUBLE PRECISION   DWORK(4*(IFFT1+IFFT2+IFFT3)+15)
      COMPLEX*16         CWORK(IFFT1+IFFT2+IFFT3)
#endif
!lh      COMPLEX*16         CFFT(IFFT1,IFFT2,IFFT3)
      COMPLEX*16         FFT(IFFT1,IFFT2,IFFT3,2)
      COMPLEX*16         UST(ifft1,ifft2,ifft3)
      COMPLEX*16         VKVXC(nkk1,2)
!
!    REAN4 reanalysation of potential using FFT method
!
! locals
      real*8             invnorm
!
#if defined (SCSL)
      ISYS(0) = 1
#endif
      PI=4.0D+0*ATAN(1.0D+0)
      TPI=2.D0*PI
      NKK=NKK1
!
!.....transform V and Ustep to real space
!
!lh   no reordering cos of the changed structure of ust(,,,)
!lh     CALL SETFF0(IFFT1,IFFT2,IFFT3,UST,CFFT)
      do i3=1,ifft3
         do i2=1,ifft2
            do i1=1,ifft1
               ust(i1,i2,i3) = conjg(ust(i1,i2,i3))
            end do
         end do
      end do
!      write(6,929) (((i1,i2,i3,fft(i1,i2,i3,1),cfft(i1,i2,i3),i1=1,3)
!     *,i2=1,3),i3=1,ifft3)
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
!lh      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,CFFT,IFFT1,
!lh                      IFFT2,CFFT,IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,ust,IFFT1,              &
                   IFFT2,ust,IFFT1,IFFT2,CWORK,DWORK,ISYS)
#else
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,1),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
!lh      CALL C3FFT(IFFT1,IFFT2,IFFT3,CFFT,IFFT1,IFFT2,
!lh     ,           ISIG,CWORK,DWORK,IERR)
!lh      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,ust,IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
#endif
!******************************************************************
!
!.....multiply ust times v
!
      DO 7 I3=1,IFFT3
         DO 7 I2=1,IFFT2
            DO 7 I1=1,IFFT1
!lh               FFT(I1,I2,I3,1)=DREAL(FFT(I1,I2,I3,1))*
!lh     *                         DREAL(CFFT(I1,I2,I3))
!lh               FFT(I1,I2,I3,2)=DREAL(FFT(I1,I2,I3,2))*
!lh     *                         DREAL(CFFT(I1,I2,I3))
               FFT(I1,I2,I3,1)=DBLE(FFT(I1,I2,I3,1))* &
                               DBLE(ust(I1,I2,I3))
               FFT(I1,I2,I3,2)=DBLE(FFT(I1,I2,I3,2))* &
                               DBLE(ust(I1,I2,I3))
    7 CONTINUE
!
!.....backtransformation to reciprocal space
!
      ISIG=+1
#if defined (SCSL)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT,IFFT1,IFFT2,        &
                   FFT,IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT(1,1,1,2),IFFT1,     &
                   IFFT2,FFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS)
#else
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT,IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
#endif

      invnorm = 1.d0/dble(IFFT1*IFFT2*IFFT3)
      DO I3=1,IFFT3
         DO I2=1,IFFT2
            DO I1=1,IFFT1
               FFT(I1,I2,I3,1)=FFT(I1,I2,I3,1) * invnorm
               FFT(I1,I2,I3,2)=FFT(I1,I2,I3,2) * invnorm
            end do
         end do
      end do
!
!
!.....extract Fouriercoefficients from fft-field
!
      CALL GETFFT(NKK,IFFT1,IFFT2,IFFT3,VKVXC,FFT,KZZ,INST)
!
      RETURN
      END
