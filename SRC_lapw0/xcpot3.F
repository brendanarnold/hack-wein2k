      SUBROUTINE XCPOT3 (NKK1,VKEXC,VKVXC,LR2V,IEX,JSPIN, &
                    IFFT1,IFFT2,IFFT3,CFFT,UFFT,FFT,POTK,imfield,bextene)
      use defs
!     use efeld
      use symetr
      use densit
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*19     nkktext
      INTEGER            IEX, IFFT1, IFFT2, IFFT3, JSPIN, NKK1
#if defined (SCSL)
      DOUBLE PRECISION   DWORK(2 * MAX(IFFT1, IFFT2, IFFT3))
      COMPLEX*16         CWORK((2*IFFT1+256)+(2*IFFT2+256)+(2*IFFT3+256))
      INTEGER            ISYS(0:1)  
#else
      DOUBLE PRECISION   DWORK(4*(IFFT1+IFFT2+IFFT3)+15)
      COMPLEX*16         CWORK(IFFT1+IFFT2+IFFT3)
#endif

      COMPLEX*16         POTK(NKK1)
      COMPLEX*16         VKEXC(NKK1,2), VKVXC(NKK1,2)

      COMPLEX*16         CFFT(IFFT1,IFFT2,IFFT3)
!      COMPLEX*16         CFFT1(IFFT1,IFFT2,IFFT3)
      COMPLEX*16         FFT(IFFT1,IFFT2,IFFT3,2)
      COMPLEX*16         UFFT(IFFT1,IFFT2,IFFT3,2)
      COMPLEX*16,ALLOCATABLE :: CFFT1(:,:,:)
      LOGICAL            LR2V
!
!    XCPOT3 calculates the exchange-correlation potential in the
!    interstitial.
!    The charge density is numerically evaluated in the unit cell.
!    The resulting exchange-potential values are transformed to a
!    Fourier-expansion by a FFT and stored in array VKVXC.
!    The array VKEXC is used in energy.
!
!    The linear equation system of the fit is solved by the
!
      INTEGER            I1, I2, I3, IERR, II, ISIG, ISPIN
      INTEGER            J, J1, JX, NKK, NKKSP
      INTEGER            KZZSP(3)
      DOUBLE PRECISION   EXD, EXU, PI, SPINDN, SPINUP, TPI, VXD, VXU
      DOUBLE PRECISION   AM(125)
      DOUBLE PRECISION   TVEC(10,27)
!
      COMPLEX*16,allocatable ::         VROKSP(:,:)
!      COMPLEX*16         VROKSP(NKK1,2)
!
      DOUBLE PRECISION   XCENER, XCPOT
      EXTERNAL           XCENER
      DATA THRD/0.3333333333333333D0/
!
#if defined (SCSL)
      ISYS(0) = 1
#endif
      allocate (VROKSP(NKK1,2),CFFT1(IFFT1,IFFT2,IFFT3))
      IF(JSPIN.EQ.2) THEN
!
!.....read plane wave spin density RHOKSP(K)
!
         READ(12,2060) NKKSP
!         READ(13,2060) NKKSP
  read(13,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkksp
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkksp
 6768 continue
         DO 4800 J=1,NKKSP
            READ(12,2071)  ( KZZSP(JX), JX=1,3 ) ,RHOKSP(J,1)
            READ(13,2071)  ( KZZSP(JX), JX=1,3 ) ,RHOKSP(J,2)
 2060       FORMAT(/,13X,I6)
 2070       FORMAT(3X,3I5,2E15.7)
 2071       FORMAT(3X,3I5,2E19.12)
 4800    CONTINUE
      ELSE
         DO 4802 J=1,2
            DO 4801 J1=1,NKK1
               RHOKSP(J1,J)=RHOK(J1)*0.5d0
 4801       CONTINUE
 4802    CONTINUE
      END IF
      PI=4.0D+0*ATAN(1.0D+0)
      TPI=2.D0*PI
      NKK=NKK1
!
!.....generate K-vectors
!
!CCCCCCCC      DO 20 ISPIN=1,JSPIN
!
      CFFT(1:IFFT1,1:IFFT2,1:IFFT3)=zeroc
      CFFT1(1:IFFT1,1:IFFT2,1:IFFT3)=zeroc
      FFT(1:IFFT1,1:IFFT2,1:IFFT3,1:2)=zeroc
      IF (IEX.EQ.45.or.iex.eq.40.or.iex.eq.41.or.iex.eq.42.or.iex.eq.12 &
                   .or.iex.eq.27.or.iex.eq.43.or.iex.eq.44)  &
         CALL VRESP (NKK1,JSPIN,CWORK,DWORK, &
         IFFT1,IFFT2,IFFT3,CFFT,CFFT1,FFT,VROKSP)
!.....iex 41: set part of tau into ufft in real space
      if(iex.eq.41.or.iex.eq.42.or.iex.eq.12.or.iex.eq.43.or.iex.eq.44  &
                  .or.iex.eq.27) then
        do j=1,nkk
        vkvxc(j,1)=potk(j)
        vkvxc(j,2)=potk(j)
        enddo
        CALL SETFF1(NKK,IFFT1,IFFT2,IFFT3,RHOKSP,UFFT)
        CALL SETFF1(NKK,IFFT1,IFFT2,IFFT3,vkvxc,FFT)
!.....transform rho and v to real space
!
!      ISIG=0
!      CALL C3FFT(IFFT1,IFFT2,IFFT3,UFFT,IFFT1,IFFT2,
!     ,           ISIG,CWORK,DWORK,IERR)
!      IF(IERR.NE.0) WRITE(*,*) 'IERR IN INITIALIZATION',IERR
      ISIG=-1
#if defined (SCSL)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,UFFT,IFFT1,IFFT2,       &
                   UFFT,IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,UFFT(1,1,1,2),IFFT1,    &
                   IFFT2,UFFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT,IFFT1,IFFT2,        &
                   FFT,IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT(1,1,1,2),IFFT1,     &
                   IFFT2,FFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS)
#else
      CALL C3FFT(IFFT1,IFFT2,IFFT3,UFFT,IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,UFFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT,IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
#endif
!
!.....multiply rho times v
!
      DO  I3=1,IFFT3
         DO  I2=1,IFFT2
            DO  I1=1,IFFT1
               UFFT(I1,I2,I3,1)=-DBLE(FFT(I1,I2,I3,1))* &
                               DBLE(UFFT(I1,I2,I3,1))
               UFFT(I1,I2,I3,2)=-DBLE(FFT(I1,I2,I3,2))* &
                               DBLE(UFFT(I1,I2,I3,2))
            enddo
         enddo
      enddo
      CALL SETFF1(NKK,IFFT1,IFFT2,IFFT3,VROKSP,FFT)
#if defined (SCSL)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT,IFFT1,IFFT2,FFT,    &
                   IFFT1,IFFT2,CWORK,DWORK,ISYS)
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
      DO  I3=1,IFFT3
         DO  I2=1,IFFT2
            DO  I1=1,IFFT1
               UFFT(I1,I2,I3,1)=DBLE(UFFT(I1,I2,I3,1))- &
                               DBLE(FFT(I1,I2,I3,1))
               UFFT(I1,I2,I3,2)=DBLE(UFFT(I1,I2,I3,2))- &
                               DBLE(FFT(I1,I2,I3,2))
            enddo
         enddo
      enddo
      endif
!.....end of option iex 41
!

!
!.....set the fft-field with rho-fourier-coefficients
!
      CALL SETFF1(NKK,IFFT1,IFFT2,IFFT3,RHOKSP,FFT)
!CCCCC
!      ISPIN=1
!      I1=1
!      I2=1
!      DO 2829 I1=1,IFFT1
!       DO 2829 I3=1,IFFT3
! 2829  WRITE(6,*) I1-1,I2-1,I3-1,FFT(I1,I2,I3,ISPIN)
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
!******************************************************************
!
!.....calculate xcpotential and energy
!
      IF(IEX.GE.10) CALL PWXAD4(IFFT1,IFFT2,IFFT3,TVEC,AM)
!*******************************
!      DO 8887 I1=1,IFFT1
!         DO 8887 I2=1,IFFT2
!            DO 8887 I3=1,IFFT3
!               DO 8887 ISPIN=1,2
!                  FFT (I1,I2,I3,ISPIN)=2.+SIN(TPI*(I3-1+IFFT3/4)/IFFT3)
!                  IF(I1.EQ.1.AND.I2.EQ.IFFT2)
!     *               WRITE(6,8886) I1,I2,I3,FFT(I1,I2,I3,1)
!8887  CONTINUE
!*******************************
!
!CCCCC
!      WRITE(90,1818) (((FFT(I1,I2,I3,1),I1=1,IFFT1),I2=1,IFFT2),
!     ,                I3=1,IFFT3)
 1818 FORMAT(10F10.5)
      DO 8 I3=1,IFFT3
         DO 8 I2=1,IFFT2
            DO 8 I1=1,IFFT1
               IF(IEX.LT.10) THEN
                  SPINUP=FFT(I1,I2,I3,1)
                  SPINDN=FFT(I1,I2,I3,2)
                  IF(SPINUP.LT.1.D-8) SPINUP=1.D-8
                  IF(SPINDN.LT.1.D-8) SPINDN=1.D-8
                  FFT(I1,I2,I3,1)=XCPOT(SPINUP+SPINDN,SPINUP,PI,IEX)
                  FFT(I1,I2,I3,2)=XCPOT(SPINUP+SPINDN,SPINDN,PI,IEX)
                  UFFT(I1,I2,I3,1)=XCENER(SPINUP+SPINDN,SPINUP,PI,IEX)
                  UFFT(I1,I2,I3,2)=XCENER(SPINUP+SPINDN,SPINDN,PI,IEX)
       if(iex.eq.5) then
           RS=(0.75D0/(PI*(SPINUP+SPINDN)))**THRD
           ZET=(SPINUP-SPINDN)/(SPINUP+SPINDN)
           CALL CORLSD(RS,ZET,ECLSD,VLSDU,VLSDD,ECRS,ECZET,ALFC)
           FFT(I1,I2,I3,1)=FFT(I1,I2,I3,1)+2.d0*VLSDU
           FFT(I1,I2,I3,2)=FFT(I1,I2,I3,2)+2.d0*VLSDD
           UFFT(I1,I2,I3,1)=UFFT(I1,I2,I3,1)+2.d0*eclsd
           UFFT(I1,I2,I3,2)=UFFT(I1,I2,I3,2)+2.d0*eclsd
!bext
                  if (imfield.gt.0)then
                   FFT(I1,I2,I3,1)=FFT(I1,I2,I3,1)-Bextene
                   FFT(I1,I2,I3,2)=FFT(I1,I2,I3,2)+Bextene
                  endif                        
        endif
!                   IF(I1.EQ.1.AND.I2.EQ.IFFT2)
!     ,               WRITE(6,8886) I1,I2,I3,FFT(I1,I2,I3,1),vlsdu
!     ,                            ,FFT(I1,I2,I3,2),UFFT(I1,I2,I3,1),
!     ,                             UFFT(I1,I2,I3,2)
               ELSE
                  IF (IEX.EQ.45) THEN
                     VXU=CFFT (I1,I2,I3)
                     VXD=CFFT1 (I1,I2,I3)
                  ENDIF
                  tauup=ufft(i1,i2,i3,1)
                  taudn=ufft(i1,i2,i3,2)
                  CALL PWXAD5(IEX,I1,I2,I3,TVEC,AM,FFT,tauup,taudn, &
                              IFFT1,IFFT2,IFFT3,VXU,VXD,EXU,EXD)
                  CFFT (I1,I2,I3)  =VXU
                  CFFT1(I1,I2,I3)  =VXD
                  UFFT (I1,I2,I3,1)=EXU
                  UFFT (I1,I2,I3,2)=EXD
!                  IF(I1.EQ.1.AND.I2.EQ.IFFT2)
!     *               WRITE(6,8886) I1,I2,I3,VXU,VXD,EXU,EXD
! 8886             FORMAT (3(1X,I3),4(3X,F10.6,1X,F10.6))
               END IF
   8  CONTINUE
!
      IF(IEX.GT.10) THEN
         DO 8888 I3=1,IFFT3
            DO 8888 I2=1,IFFT2
               DO 8888 I1=1,IFFT1
                  FFT(I1,I2,I3,1)=CFFT (I1,I2,I3)
                  FFT(I1,I2,I3,2)=CFFT1(I1,I2,I3)
!bext
                  if (imfield.gt.0)then
                   FFT(I1,I2,I3,1)=FFT(I1,I2,I3,1)-Bextene
                   FFT(I1,I2,I3,2)=FFT(I1,I2,I3,2)+Bextene
                  endif                                                
 8888    CONTINUE
!         WRITE(91,1818) (((FFT(I1,I2,I3,1),I1=1,IFFT1),I2=1,IFFT2),
!     ,                   I3=1,IFFT3)
      ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCC
!       DO 7781 I1=1,IFFT1
!          DO 7781 I2=1,IFFT2
!             WRITE(66,7771) (FFT(I1,I2,I3,1),I3=1,IFFT3)
! 7781  CONTINUE
! 7771  FORMAT(18F8.3)
!CCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!.....backtransformation to reciprocal space
!
      ISIG=+1
#if defined (SCSL)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT(1,1,1,1),IFFT1,     &
                   IFFT2,FFT(1,1,1,1),IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT(1,1,1,2),IFFT1,     &
                   IFFT2,FFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,UFFT(1,1,1,1),IFFT1,    &
                   IFFT2,UFFT(1,1,1,1),IFFT1,IFFT2,CWORK,DWORK,ISYS)
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,UFFT(1,1,1,2),IFFT1,    &
                   IFFT2,UFFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS)
#else
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,1),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,UFFT(1,1,1,1),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,UFFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
#endif
       fftin = 1.d0/dble(IFFT1*IFFT2*IFFT3)
      DO 78 I3=1,IFFT3
         DO 78 I2=1,IFFT2
            DO 78 I1=1,IFFT1
               FFT(I1,I2,I3,1)=FFT(I1,I2,I3,1) * fftin
               FFT(I1,I2,I3,2)=FFT(I1,I2,I3,2) * fftin
               UFFT(I1,I2,I3,1)=UFFT(I1,I2,I3,1) * fftin
               UFFT(I1,I2,I3,2)=UFFT(I1,I2,I3,2) * fftin
  78  CONTINUE
      WRITE(6,2040)
      WRITE(6,2050) (0,0,I3-1,FFT(1,1,I3,1),FFT(1,1,I3,2), &
                     I3=1,IFFT3)
 2040 FORMAT(/,' SELECTED FOURIERCOEFF. OF V-XC ',I2)
 2050 FORMAT(3I3,2E15.5,5x,2e15.5)
!
!.....extract Fouriercoefficients from fft-field
!
      CALL GETFF1(NKK,IFFT1,IFFT2,IFFT3,VKVXC,FFT,KZZ,INST)
      CALL GETFF1(NKK,IFFT1,IFFT2,IFFT3,VKEXC,UFFT,KZZ,INST)
!
!!!      if(abs(iefeld).gt.0) then
!!!!     set potential to zero if iefeld negative
!!!       if(iefeld.lt.0) then
!!!          do i=1,nkk
!!!          vkvxc(i,1)=0.d0
!!!          vkvxc(i,2)=0.d0
!!!          enddo
!!!        endif
!!!        if(abs(iefeld).lt.999) then
!!!         do i=1,abs(iefeld)
!!!          vkvxc(jefeld(i),1)=vkvxc(jefeld(i),1)+4.d0/pi**2*refeld/ &
!!!                         kzz(3,jefeld(i))**2
!!!          vkvxc(jefeld(i),2)=vkvxc(jefeld(i),2)+4.d0/pi**2*refeld/ &
!!!                         kzz(3,jefeld(i))**2
!!!         enddo
!!!        else
!!!!  efeld 7:8/1:8 geteilt
!!!         iefeld1=mod(abs(iefeld),1000)
!!!         do j=1,mod(abs(iefeld),1000)
!!!          i=kzz(3,jefeld(j))
!!!          vkvxc(jefeld(j),1)=vkvxc(jefeld(j),1)+                            & 
!!!                 refeld*((8.d0*cos(i*pi/4.d0)+i*pi*sin(i*pi/4.d0)-8.d0) &
!!!                      /2.d0/pi**2/i**2 + &
!!!               (8.d0*cos(7.d0/4.d0*pi*i)+7.d0*i*pi*sin(7.d0/4.d0*pi*i)-8.d0) &
!!!                    /14.d0/i**2/pi**2 )           
!!!          vkvxc(jefeld(j),1)=vkvxc(jefeld(j),1)+                            & 
!!!               refeld*((-8.d0*sin(pi/4.d0*i)+i*pi*cos(pi/4.d0*i)+i*pi) &
!!!                      /2.d0/i**2/pi**2  &
!!!                    -(-8.d0*sin(7.d0/4.d0*pi*i)+7.d0*i*pi*cos(7.d0/4.d0*i*pi)&
!!!                    +7.d0*i*pi)/14.d0/i**2/pi**2 ) / (0.d0,1.d0)
!!!          vkvxc(jefeld(j),2)=vkvxc(jefeld(j),2)+                            & 
!!!                 refeld*((8.d0*cos(i*pi/4.d0)+i*pi*sin(i*pi/4.d0)-8.d0) &
!!!                      /2.d0/pi**2/i**2 + &
!!!               (8.d0*cos(7.d0/4.d0*pi*i)+7.d0*i*pi*sin(7.d0/4.d0*pi*i)-8.d0) &
!!!                    /14.d0/i**2/pi**2 )           
!!!          vkvxc(jefeld(j),2)=vkvxc(jefeld(j),2)+                            & 
!!!               refeld*((-8.d0*sin(pi/4.d0*i)+i*pi*cos(pi/4.d0*i)+i*pi) &
!!!                      /2.d0/i**2/pi**2  &
!!!                    -(-8.d0*sin(7.d0/4.d0*pi*i)+7.d0*i*pi*cos(7.d0/4.d0*i*pi)&
!!!                    +7.d0*i*pi)/14.d0/i**2/pi**2 ) / (0.d0,1.d0)
!!!         enddo
!!!        endif
!!!      endif
!
!.....set the fft-field with fourier-coefficients
!.....for total energy (integration of interstital)
!
      CALL SETFFT(NKK,IFFT1,IFFT2,IFFT3,RHOKSP,FFT)
!
      if(iex.eq.40) then
!.....put V-xc + V-coul into rhoksp, it will be multiplied with rho
!.....and later VKEXC(J,1)*.25d0-VROKSP(j,1)) will be added
      write(6,*) ' Local kinetic energy in xcpot3'
      write(6,*) VROKSP(1,1),RHOKSP(1,1),VROKSP(nkk,1),RHOKSP(nkk,1)
      DO  J=1,NKK
!         RHOKSP(J,1)=(-VKEXC(J,1)*.25d0-VROKSP(j,1))+
!     *                VKVXC(J,1)+POTK(J)
!         RHOKSP(J,2)=(-VKEXC(J,2)*.25d0-VROKSP(j,2))+
!     *                VKVXC(J,2)+POTK(J)
!         if(real(RHOKSP(J,1)).gt.1.d-6.and.
!     *      real(RHOKSP(J,2)).gt.1.d-6) then
         RHOKSP(J,1)=(VKVXC(J,1)+POTK(J))
         RHOKSP(J,2)=(VKVXC(J,2)+POTK(J))
!         endif
      enddo
!
      else if(iex.eq.30.or.iex.eq.31.or.iex.eq.32.or.iex.eq.33.or. &
      iex.eq.34.or.iex.eq.41.or.iex.eq.42) then
      write(6,*) ' Exchange energy only in xcpot3'
      DO 222 J=1,NKK
         RHOKSP(J,1)=-VKEXC(J,1)
         RHOKSP(J,2)=-VKEXC(J,2)
 222  CONTINUE
      else  
!
      DO 22 J=1,NKK
         RHOKSP(J,1)=VKVXC(J,1)-VKEXC(J,1)+POTK(J)*0.5D0
         RHOKSP(J,2)=VKVXC(J,2)-VKEXC(J,2)+POTK(J)*0.5D0
 22   CONTINUE
!
      endif  
!
      CALL SETFFT(NKK,IFFT1,IFFT2,IFFT3,RHOKSP,UFFT)

!      WRITE(82,*) FFT
!      WRITE(83,*) UFFT
!
!.....transform rho and v to real space
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
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT,IFFT1,IFFT2,        &
                   FFT,IFFT1,IFFT2,CWORK,DWORK,ISYS) 
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,FFT(1,1,1,2),IFFT1,     &
                   IFFT2,FFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS) 
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,UFFT,IFFT1,IFFT2,       &
                   UFFT,IFFT1,IFFT2,CWORK,DWORK,ISYS) 
      CALL ZZFFT3D(ISIG,IFFT1,IFFT2,IFFT3,1.0D0,UFFT(1,1,1,2),IFFT1,    &
                   IFFT2,UFFT(1,1,1,2),IFFT1,IFFT2,CWORK,DWORK,ISYS) 
#else
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT,IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,FFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,UFFT,IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
      CALL C3FFT(IFFT1,IFFT2,IFFT3,UFFT(1,1,1,2),IFFT1,IFFT2, &
                 ISIG,CWORK,DWORK,IERR)
      IF(IERR.NE.0) WRITE(*,*) 'IERR IN TRANSFORMATION',IERR
#endif
!
!.....multiply rho times v
!
      DO 7 I3=1,IFFT3
         DO 7 I2=1,IFFT2
            DO 7 I1=1,IFFT1
               FFT(I1,I2,I3,1)=DBLE(FFT(I1,I2,I3,1))* &
                               DBLE(UFFT(I1,I2,I3,1))
               FFT(I1,I2,I3,2)=DBLE(FFT(I1,I2,I3,2))* &
                               DBLE(UFFT(I1,I2,I3,2))
    7 CONTINUE

!
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
      DO 778 I3=1,IFFT3
         DO 778 I2=1,IFFT2
            DO 778 I1=1,IFFT1
               FFT(I1,I2,I3,1)=FFT(I1,I2,I3,1) * fftin
               FFT(I1,I2,I3,2)=FFT(I1,I2,I3,2) * fftin
!
!.....store pw-total energy contributions in cfft
!
               CFFT(I1,I2,I3)=FFT(I1,I2,I3,1)+FFT(I1,I2,I3,2)

!CCCC
!               IF(I1.EQ.1.AND.I2.EQ.1) WRITE(6,2052)
!     *            CFFT(I1,I2,I3),FFT(I1,I2,I3,1),FFT(I1,I2,I3,2)
  778 CONTINUE
      if(iex.eq.40) then
      DO  J=1,NKK
         RHOKSP(J,1)=(-VKEXC(J,1)*.25d0+VROKSP(j,1))
         RHOKSP(J,2)=(-VKEXC(J,2)*.25d0+VROKSP(j,2))
!         RHOKSP(J,1)=(0.d0)
!         RHOKSP(J,2)=(0.d0)
      enddo
      CALL SETFF2(NKK,IFFT1,IFFT2,IFFT3,RHOKSP(1,1),CFFT)
      CALL SETFF2(NKK,IFFT1,IFFT2,IFFT3,RHOKSP(1,2),CFFT)
      endif
!
!.....put total potentials in fft
!
      DO 23 J=1,NKK
         RHOKSP(J,1)=VKVXC(J,1)+POTK(J)
         RHOKSP(J,2)=VKVXC(J,2)+POTK(J)
 23   CONTINUE

!     Add E-field terms if needed now
      if(abs(iefeld).gt.0)call epot4(rhoksp,nkk)
!
      CALL SETFFT(NKK,IFFT1,IFFT2,IFFT3,RHOKSP,FFT)
!CCCCC
!      WRITE(6,2052) ((FFT(1,1,J,J1),J=1,IFFT3),J1=1,2)
 2052 FORMAT(6F15.9)

!
!.....printouts
!
!
      DEALLOCATE (vroksp,cfft1)
      RETURN
!
  998 FORMAT(3I5,2E15.7)
 1000 FORMAT(3X,'PW VOLUME IN PERCENT=',F10.5)
 1513 FORMAT(9X,I4)
 1530 FORMAT(3X,'CONDITION PARAMETER FOR INTERSTITIAL    INFO= ',I2)
      END
