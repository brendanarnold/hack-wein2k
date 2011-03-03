      SUBROUTINE WFPNT1 (sx,sy,sz,NV,NE,NUME)
      use matrices, only: eigval, kzz, Z
      use parallel, only: DESCZ, MYID, ICTXTALL, DESCLOCAL
      IMPLICIT NONE
      INTEGER            LDZ, NE, NV, NUME
      real*8             sx,sy,sz
!_REAL      DOUBLE PRECISION, allocatable ::    ZLOCAL(:,:)
!_COMPLEX      complex*16, allocatable ::       ZLOCAL(:,:)
!
!     ..................................................................
!
!        print output editor for non-relativistic wavefunctions
!        (cubic version)                     P.I.SORANTIN
!
!     ..................................................................
!
      INTEGER            NCOLS
      PARAMETER          (NCOLS = 1)
!_COMPLEX      COMPLEX*16         IMAG
!_COMPLEX      PARAMETER          (IMAG = (0.0D+0,1.0D+0))
!
      INTEGER            I, ILOWER, IUPPER, J, K
!
      IF (NE .GT. NUME) THEN
         WRITE(6,9000) NUME
         GOTO 999
      ENDIF
!
!        print out all computed Eigenvectors stored in 'Z'
!        formatted in 'NCOLS' columns
!
#ifdef Parallel
!  Note: Parallel Implementation: Possible improvement: see comment in tapewf.F
      IF(MYID.EQ.0) allocate( ZLOCAL(NV, NE) )
!_REAL      CALL PDGEMR2D( NV, NE, Z, 1, 1, DESCZ, ZLOCAL, 1, 1, DESCLOCAL, ICTXTALL)
!_COMPLEX      CALL PZGEMR2D( NV, NE, Z, 1, 1, DESCZ, ZLOCAL, 1, 1, DESCLOCAL, ICTXTALL)
#endif
      WRITE(55,6030)
      write(6,6000) sx,sy,sz,ne,nv
      write(55,6000) sx,sy,sz,ne,nv
         DO   J = 1, NV
            WRITE(55,6010) KZZ(1,J), KZZ(2,J), KZZ(3,J)
         enddo
      DO 10 ILOWER = 1, NE, NCOLS
         IUPPER = MIN(NE,ILOWER+NCOLS-1)
         WRITE(55,6005) eigval(ilower)
!         WRITE(6,6040) (I,I=ILOWER,IUPPER)
!         DO 30  J = 1, NV
!            WRITE(6,6010) KZZ(1,J), KZZ(2,J), KZZ(3,J)
#ifdef Parallel
!_REAL            WRITE(55,6020) (ZLOCAL(J,ilower),j=1,nv)
!C_COMPLEX            WRITE(6,6020) (DBLE(ZLOCAL(J,K)),K=ILOWER,IUPPER)
!C_COMPLEX            WRITE(6,6021) (DBLE(ZLOCAL(J,K)*IMAG),K=ILOWER,IUPPER)
#else
!_REAL            WRITE(55,6020) (Z(J,ilower),j=1,nv)
!_COMPLEX!
!_COMPLEX!        print the real-part of the eigenvectors
!_COMPLEX!
!C_COMPLEX            WRITE(6,6020) (DBLE(Z(J,K)),K=ILOWER,IUPPER)
!_COMPLEX!
!_COMPLEX!        print the imaginary-part of the eigenvectors
!_COMPLEX!
!C_COMPLEX            WRITE(6,6021) (DBLE(Z(J,K)*IMAG),K=ILOWER,IUPPER)
#endif
   30    CONTINUE
   10 CONTINUE
#ifdef Parallel
      IF(MYID.EQ.0) deallocate( ZLOCAL )
#endif
!
  999 RETURN
!
 6000 FORMAT(/,3f10.6,2i5)
 6010 FORMAT(1X,3i5)
6005  format(f20.8)
 6020 FORMAT(8F9.6)
!_COMPLEX 6021 FORMAT(18X,9(F9.6,2X),5X,'IMAGPART')
 6030 FORMAT(1X)
 6040 FORMAT(18X,9(I3,'.ENERGY '))
 9000 FORMAT(' NUMBER OF ENERGIES GT ',I3,' NUME IN BASISO')
!
!        End of 'WFPNT'
!
      END
