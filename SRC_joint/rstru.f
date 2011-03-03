!
!......subroutine to read structfile :   for eps2 factor
!............                            to be able to symetrize
!...... should also be able to supply the bravaismatrices!
!
      SUBROUTINE rstru
!          FILE 20          ...   STRUCTURE
!      INCLUDE 'param.inc'
      use felder
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8       AAA(3,3),A,OPM(3,3)
      real*4       O
      DIMENSION        iii(3,3)
      CHARACTER*4      IREL,RELA
!      CHARACTER*5      MODUS
!      CHARACTER*10     ANAME
      CHARACTER*11     STATUS,FORM
      CHARACTER*77     RMARGN
      CHARACTER*80     FNAME
      LOGICAL          REL
      INTEGER          S,SI
!
!      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
!                      PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO), &
!                      JRI(NATO),R0(NATO)
!      COMMON /COM/    EMIN,EMAX,ELECN,IOUT,NSPIN,EULIMIT,EDLIMIT, &
!                      NAT,NBAND,REL,NB(NKPT),MINWAV,MAXWAV,ix,NK
!      COMMON /CHAR/   TITLE,LATTIC,MODUS,ANAME(NATO)
      COMMON /GENER/  BR1(3,3),BR2(3,3)
!      COMMON /ROTMAT/ ROTLOC(3,3,NATO),S(3,3,48),TAU(3,48), &
!                      SI(3,3,48)
!.....READ STRUCTFILE...............................................

      READ(20,1000) TITLE
      READ(20,1010) LATTIC,NAT,IREL
      allocate (ANAME(nat),RMT(nat),Vatom(nat),R0(nat),ROTLOC(3,3,NAT))
      allocate (IATNR(nat),MULT(nat),ISPLIT(nat),JRI(nat))
      allocate (POS(3,48*nat))               ! needs only ndif 
!.....READ IN LATTICE CONSTANTS
      READ(20,1020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
      IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
      IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
!C
!      if(nat.gt.nato) stop 'NATO TOO SMALL!'
      WRITE(6,800)    
      WRITE(6,805)  TITLE
      WRITE(6,810)  LATTIC
      WRITE(6,820)  AA,BB,CC,(alpha(i),i=1,3)
      WRITE(6,840)  NAT
      WRITE(6,850)  IREL
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE
!     NOTEQUIVALENT ATOMS
      INDEX=0
      DO 50 JATOM = 1,NAT
         INDEX=INDEX+1
         READ(20,1030) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM), &
                       ISPLIT(JATOM)
         IF (MULT(JATOM).EQ.0) THEN
            WRITE(6,1040) JATOM,INDEX,MULT(JATOM)
            STOP ' OPTIC: MULT EQ 0'
         ENDIF
            DO 55 M=1,MULT(JATOM)-1
               INDEX=INDEX+1
!      if(index.gt.NDIF) stop 'NDIF TOO SMALL!'
               READ(20,1031) IATNR(JATOM),( POS(J,INDEX),J=1,3)
 55         CONTINUE  
         READ(20,1050) ANAME(JATOM),JRI(JATOM),R0(JATOM),RMT(JATOM)
         READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)
 1051    FORMAT(20X,3F10.8)
 50   CONTINUE
      ndif=index
      allocate (rotij(3,3,index),TAUIJ(3,index))
!.......to take latgen.f of lapw2 package, sym is read in rotdef!.....
      goto 999
!!..........................!! skip following !!.........................
!      READ(20,10) IORD
!      DO 5 J=1,IORD   
!      READ(20,11) ( (S(J2,J1,J),J1=1,3),TAU(J2,J), J2=1,3 )
! 5    CONTINUE
!!.....READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS
! 10   FORMAT(I4)
! 11   FORMAT(3(3I2,F5.2/))
!!.....go to latgen for br1 and vol........
!!k
!!... use same as in optic but delete last liene containing ROTDEF !!
!!k
!      CALL LATGEN(NAT)
!!.....PRODUCE INVERSE SYMETRY MATRIX
!      DO J=1,IORD    
!      WRITE(6,121) J , iord
! 121  FORMAT(/,'....... OPERATION:',I2,'[',i2,']')
!        CALL INVERSSYMDEF(S(1,1,J),SI(1,1,J))
!      WRITE(6,11) ( (S(J2,J1,J),J1=1,3),TAU(J2,J), J2=1,3 )
!        WRITE(6,*) 'SYMMERTY CHECK......'
!      WRITE(6,11) ( (SI(J2,J1,J),J1=1,3),TAU(J2,J), J2=1,3 )
!        WRITE(6,*) 'SYMMERTY CHECK BY MULTIPLICATION S*SI ......'
!      do i=1,3
!      do k=1,3
!       isi=0
!         do l=1,3
!            isi=isi+s(i,l,j)*si(l,k,j)
!         end do
!       iii(i,k)=isi
!      end do
!      end do
!      WRITE(6,11) ( (III(J2,J1),J1=1,3),TAU(J2,J), J2=1,3 )
!!        WRITE(6,222) ((S(L,M,J),SI(L,M,J), M=1,3), L=1,3) ,J
! 222  FORMAT(3(5X,3I2,4X,3I2/),'   ....',I2)
!      END DO
!!......my sym skiped!!................................................
 999  continue
      CALL LATGEN(NAT)
!     write(6,*) ' bravais matrix used BR2:'
!     WRITE(6,12) ( (J2,J1,br2(J2,J1),J1=1,3), J2=1,3 )
!12   FORMAT(3(3(2I1,F12.7,2x)/))
 800  FORMAT(//,30X,50(1H-),/,33X,'S T R U C T U R A L   ', &
             'I N F O R M A T I O N',/,30X,50(1H-),//)
 805  FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)
 810  FORMAT(3X,'LATTICE',22X,'= ',A4)
 820  FORMAT(3X,'LATTICE PARAMETERS ARE',7X,'= ',3F12.7,3F7.2)
 830  FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)
 840  FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)
 850  FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)
 860  FORMAT(3X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)
 1000 FORMAT(A80)
 1010 FORMAT(A4,23X,I3,/,13X,A4)
 1020 FORMAT(6F10.5)
 1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
 1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
 1040 FORMAT(///,3X,'ERROR IN LAPW2 : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 1002 FORMAT(3F10.5,I5)
 1003 FORMAT(A5)

      RETURN
      END
