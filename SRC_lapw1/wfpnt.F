      SUBROUTINE WFPNT (NV,NE,NUME)
      use matrices, only: KZZ, Z
      use parallel, only: MYID, DESCZ, ICTXTALL, DESCLOCAL
      IMPLICIT NONE
      INTEGER            LDZ, NE, NV, NUME
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
      PARAMETER          (NCOLS = 9)
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
      WRITE(6,6030)
      DO 10 ILOWER = 1, NE, NCOLS
         IUPPER = MIN(NE,ILOWER+NCOLS-1)
         WRITE(6,6000)
         WRITE(6,6040) (I,I=ILOWER,IUPPER)
         DO 30  J = 1, NV
            WRITE(6,6010) KZZ(1,J), KZZ(2,J), KZZ(3,J)
#ifdef Parallel
!_REAL            WRITE(6,6020) (ZLOCAL(J,K),K=ILOWER,IUPPER)
!_COMPLEX            WRITE(6,6020) (DBLE(ZLOCAL(J,K)),K=ILOWER,IUPPER)
!_COMPLEX            WRITE(6,6021) (DBLE(ZLOCAL(J,K)*IMAG),K=ILOWER,IUPPER)
#else
!_REAL            WRITE(6,6020) (Z(J,K),K=ILOWER,IUPPER)
!_COMPLEX!
!_COMPLEX!        print the real-part of the eigenvectors
!_COMPLEX!
!_COMPLEX            WRITE(6,6020) (DBLE(Z(J,K)),K=ILOWER,IUPPER)
!_COMPLEX!
!_COMPLEX!        print the imaginary-part of the eigenvectors
!_COMPLEX!
!_COMPLEX            WRITE(6,6021) (DBLE(Z(J,K)*IMAG),K=ILOWER,IUPPER)
#endif
   30    CONTINUE
   10 CONTINUE
#ifdef Parallel
      IF(MYID.EQ.0) deallocate( ZLOCAL )
#endif
!
  999 RETURN
!
 6000 FORMAT(/,/,/,3X,'RECIPROCAL LATTICE VECTORS',/)
 6010 FORMAT(1X,I3,2I4)
 6020 FORMAT(18X,9(F9.6,2X),5X,'REALPART')
!_COMPLEX 6021 FORMAT(18X,9(F9.6,2X),5X,'IMAGPART')
 6030 FORMAT(1X)
 6040 FORMAT(18X,9(I3,'.ENERGY '))
 9000 FORMAT(' NUMBER OF ENERGIES GT ',I3,' NUME IN BASISO')
!
!        End of 'WFPNT'
!
      END
