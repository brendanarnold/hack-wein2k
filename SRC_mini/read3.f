      subroutine read3(nkknew,iterho,nat,jspin)
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!
      CHARACTER*19     nkktext
      COMPLEX*16     rokmix,ROKNEW,ROKOLD,ROKOLD1,ROKOLD2
!                                                                       
      COMMON /SYMETR/  TAU(3,NSYM),KVCNEW(3,NWAV),IZ(3,3,NSYM), &
                       INUM(NSYM),IORD
      COMMON /XA/ CLMOLD(NRAD,NCOM,NATO,2),CLMOLD1(NRAD,NCOM,NATO,2), &
                CLMOLD2(NRAD,NCOM,NATO,2),CLMNEW(NRAD,NCOM,NATO,2), &
                ROKNEW(NWAV,2),ROKOLD(NWAV,2), &
                ROKOLD1(NWAV,2),ROKOLD2(NWAV,2),rokmix(NWAV,2), &
                R(NRAD),RNOT(NATO),              &
                DX(NATO),ZZ(NATO),JRI(NATO),                       &
                LM(2,NCOM,NATO),LMMAX(NATO)   
!
!        Common blocks
!
!
!        LATTIC  - lattice type
!        NAME(i) - name of atom i in the unit-cell (e.g. 'Titanium')
!
      CHARACTER*60        test
      CHARACTER*4        LATTIC
      CHARACTER*10       NAME(NATO)
      COMMON  /CHAR/     LATTIC, NAME
      SAVE    /CHAR/
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        IATNR(i)    - atom index of (inequivalent) atom i
!                      also indicates cubic and non-cubic symmetries
!                      IATNR(i) .GT. 0 ... cubic symmetry
!                      IATNR(i) .LT. 0 ... non-cubic symmetry
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - volume of the direct lattice unit-cell
!
      INTEGER            IATNR(NATO), MULT(NATO)
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3), POS(3,NDIF)
      DOUBLE PRECISION   RMT(NATO), V(NATO)
      COMMON  /STRUK/    POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      SAVE    /STRUK/
!---------------------------------------------------------------------  
!                                                                       
!                                                                       
      iterho1=iterho
      write(6,*)'iterho1',iterho
      if (iterho.gt.0) then
      DO 60 ISPIN=1,JSPIN
         READ(13,2032)
         DO 51 JATOM=1,NAT
            READ(13,1980)
            READ(13,2000)lmmax1
            DO 50 LM1=1,lmmax1
               READ(13,2010) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
               READ(13,2021) (CLMOLD(J,LM1,JATOM,ispin),J=1,JRI(JATOM) )
  50     READ(13,2031)
  51     READ(13,2033)
         READ(13,1980)
!         READ(13,2060) NKKNEW
  read(13,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkknew
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkknew
 6768 continue
         DO 70 J=1,NKKNEW
   70    READ(13,2071) (KVCNEW(JX,J),JX=1,3),ROKOLD(J,ispin)
   60 CONTINUE
      endif    
      if (iterho.gt.1) then
      write(6,*)'iterho2',iterho
      DO 160 ISPIN=1,JSPIN
         READ(13,2032)
         DO 151 JATOM=1,NAT
            READ(13,1980)
            READ(13,2000)lmmax1
            DO 150 LM1=1,lmmax1
               READ(13,2031)
               READ(13,2021)(CLMOLD1(J,LM1,JATOM,ispin),J=1,JRI(JATOM))
  150    READ(13,2031)
  151    READ(13,2033)
         READ(13,1980)
!         READ(13,2060) NKKNEW
  read(13,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6769) nkknew
  goto 6770
 6769 read(nkktext,'(13x,i6)') nkknew
 6770 continue
         DO 170 J=1,NKKNEW
  170    READ(13,2071) (KVCNEW(JX,J),JX=1,3),ROKOLD1(J,ispin)
  160 CONTINUE
      endif    
      if (iterho.gt.2) then
      write(6,*)'iterho3',iterho
      DO 260 ISPIN=1,JSPIN
         READ(13,2032)
         DO 251 JATOM=1,NAT
            READ(13,1980)
            READ(13,2000)lmmax1
            DO 250 LM1=1,lmmax1
               READ(13,2031)
               READ(13,2021)(CLMOLD2(J,LM1,JATOM,ispin),J=1,JRI(JATOM))
  250    READ(13,2031)
  251    READ(13,2033)
         READ(13,1980)
!         READ(13,2060) NKKNEW
  read(13,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6771) nkknew
  goto 6772
 6771 read(nkktext,'(13x,i6)') nkknew
 6772 continue
         DO 270 J=1,NKKNEW
  270    READ(13,2071) (KVCNEW(JX,J),JX=1,3),ROKOLD2(J,ispin)
  260 CONTINUE
      endif    
  119 CONTINUE
      if (iterho.lt.3) iterho=iterho+1
!
! write new densities
!
      write(6,*)'bin hier',JSPIN,nat,JRI(1),LMMAX(1)
      rewind(13)
      write(13,*) iterho
      write(6,*)'bin hier',JSPIN,nat,JRI(1),LMMAX(1)
      DO 360 ISPIN=1,JSPIN
         WRITE(13,2040)-1
         DO 351 JATOM=1,NAT
            WRITE(13,1980)
            WRITE(13,2000)lmmax(jatom)
            DO 350 LM1=1,lmmax(jatom)
               WRITE(13,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
               WRITE(13,2021)(CLMNEW(J,LM1,JATOM,ispin),J=1,JRI(JATOM))
  350    WRITE(13,2031)
  351    WRITE(13,2033)
         WRITE(13,1980)
         WRITE(13,2060) NKKNEW
         DO 370 J=1,NKKNEW
  370    WRITE(13,2071) (KVCNEW(JX,J),JX=1,3),ROKNEW(J,ispin)
  360 CONTINUE
      if(iterho.gt.1) then
      DO 460 ISPIN=1,JSPIN
         WRITE(13,2040)-2
         DO 451 JATOM=1,NAT
            WRITE(13,1980)
            WRITE(13,2000)lmmax(jatom)
            DO 450 LM1=1,lmmax(jatom)
               WRITE(13,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
               WRITE(13,2021)(CLMOLD(J,LM1,JATOM,ispin),J=1,JRI(JATOM))
  450    WRITE(13,2031)
  451    WRITE(13,2033)
         WRITE(13,1980)
         WRITE(13,2060) NKKNEW
         DO 470 J=1,NKKNEW
  470    WRITE(13,2071) (KVCNEW(JX,J),JX=1,3),ROKOLD(J,ispin)
  460 CONTINUE
      endif    
      if(iterho.gt.2) then
      DO 560 ISPIN=1,JSPIN
         WRITE(13,2040)-3
         DO 551 JATOM=1,NAT
            WRITE(13,1980)
            WRITE(13,2000)lmmax(jatom)
            DO 550 LM1=1,lmmax(jatom)
               WRITE(13,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
               WRITE(13,2021)(CLMOLD1(J,LM1,JATOM,ispin),J=1,JRI(JATOM))
  550    WRITE(13,2031)
  551    WRITE(13,2033)
         WRITE(13,1980)
         WRITE(13,2060) NKKNEW
         DO 570 J=1,NKKNEW
  570    WRITE(13,2071) (KVCNEW(JX,J),JX=1,3),ROKOLD1(J,ispin)
  560 CONTINUE
      endif    
      iterho=iterho1
      RETURN
 1980 FORMAT(3X)
 2000 FORMAT(15X,I3//)
 2010 FORMAT(15X,I3,5X,I2/)
 2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)
 2021 FORMAT(3X,4E19.12)
 2031 FORMAT(/)
 2032 FORMAT(//)
 2033 FORMAT(///)
 2040 FORMAT('TOTAL CHARGE DENSITY GENERATED BY',i3,' ITERATION',//)
 2060 FORMAT(/,13X,I6)
 2071 FORMAT(3X,3I5,2E19.12)
      END
