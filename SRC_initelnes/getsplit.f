      SUBROUTINE GETSPLIT(NATOM, OurISPLIT)
!     Gives ISPLIT of the ato NATOM in OurISPLIT.

      INTEGER NATOM, OurISPLIT

!     lecture of case.struct:
      CHARACTER*80     TITLE
      CHARACTER*4      LATTIC
      CHARACTER*10     NAME
      CHARACTER*4      IREL
      INTEGER          INDEX, JATOM, NAT, I1, I2, M
      DOUBLE PRECISION R0, RMT, ZZ
      DOUBLE PRECISION ROTLOC(3,3)
      DOUBLE PRECISION PosX, PosY, PosZ
      DOUBLE PRECISION Al, Bl, Cl, Alpha, Beta, Gamma
      INTEGER          JRJ
      INTEGER          IATNR, MULT, ISPLIT
      
      
      REWIND(20)
!     we search the ISPLIT of the atom number: we get OurISPLIT
      READ (20,5000) TITLE
      READ (20,5010) LATTIC, NAT, IREL
      READ (20,5020) Al, Bl, Cl, Alpha, Beta, Gamma
      INDEX = 0
      DO JATOM = 1,NAT
         INDEX = INDEX + 1
         READ (20,5030) IATNR, PosX, PosY, PosZ, MULT, ISPLIT
         IF (NATOM.EQ.JATOM) THEN
            OurISPLIT = ISPLIT
            REWIND (20)
            RETURN
         ENDIF
         DO 40 M = 1, MULT - 1
            INDEX = INDEX + 1
            READ (20,5040) IATNR, PosX, PosY, PosZ
   40    CONTINUE
         READ (20,5050) NAME, JRJ, R0, RMT, ZZ
         READ (20,5060) ((ROTLOC(I1,I2),I1=1,3),I2=1,3)
      ENDDO

 5000 FORMAT(A80)
 5010 FORMAT(A4,23X,I3,/,13X,A4)
 5020 FORMAT(6F10.7)
 5030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
 5040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
 5050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 5060 FORMAT(20X,3F10.8)
 5070 FORMAT(A5)
 5080 FORMAT(F10.7,2I5)
 5090 FORMAT(20X,I1)
 5100 FORMAT(A10,4I5,3F5.2,A3)
      
      RETURN
      END
