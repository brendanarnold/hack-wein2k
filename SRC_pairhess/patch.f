!     Correct for symmetry errors
!     Version 1.0, Feb 2006
!
!     Read the structure file
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pos(3,48*natmax),pred(3,natmax)
      dimension inum(48)
      CHARACTER*10    ANAME(natmax)
      dimension iatnr(natmax),ROTLOC(3,3,natmax),jri(natmax)
      dimension rnot(natmax),rmt(natmax),zatom(natmax),isplit(natmax)
!      character*1 junk
      character*80 title1
      dimension alpha(3)
      character*6 lattic
      character*12 cform
      character*80 fname,ename
!      character*23 stuff
      character*11 status,form
!      character*67 ERRMSG

!     Name
      call gtfnam(fname,ename)
!      CALL ERRFLG(Ename,'Error in PAIRHESS')
      OPEN (1,FILE=fname,STATUS='OLD')
      iopen=0
8000  CONTINUE
      READ (1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM)
      GOTO 8000
8001  continue
!.....READ STRUCT
      READ(20,1000) TITLE1
      READ(20,1010) LATTIC,NAT,cform,IREL
!
      N1=nat
!     READ IN LATTICE CONSTANTS
      READ(20,1020) alat,ALPHA
      IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
      IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
      IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
!
!     INDeX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE
!     NOTEQUIVALENT ATOMS
      INDEX=0
      TOL=1D-6
      DO 50 JATOM = 1,NAT
        INDEX=INDEX+1
        READ(20,2030) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM),  &
          ISPLIT(JATOM)
        DO J=1,3
                if(POS(J,INDEX) .lt.     -TOL)POS(J,INDEX)=POS(J,INDEX)+1.D0
                if(POS(J,INDEX) .gt. 1.D0-TOL)POS(J,INDEX)=POS(J,INDEX)-1.D0
        enddo
        DO 55 M=1,MULT(JATOM)-1
          INDEX=INDEX+1
          READ(20,2031) IATNR(JATOM),( POS(J,INDEX),J=1,3)
        DO J=1,3
                if(POS(J,INDEX) .lt.     -TOL)POS(J,INDEX)=POS(J,INDEX)+1.D0
                if(POS(J,INDEX) .gt. 1.D0-TOL)POS(J,INDEX)=POS(J,INDEX)-1.D0
        enddo
   55   CONTINUE
          READ(20,1050) ANAME(JATOM),JRI(JATOM),RNOT(JATOM),RMT(JATOM)      &
          ,zatom(jatom)     
          READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)                
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)   
 1051     FORMAT(20X,3F10.7)                              
   50 CONTINUE
!
      N2=index
!     Get symmetry
      READ(20,1151) iord
1151  FORMAT(I4)

      DO j=1,iord
         READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
1101     FORMAT(3(3I2,F10.8/),I8)
      ENDDO
      call fixup2(pos)

!     new struct-file on tape 21

      WRITE(21,1000) TITLE1
      WRITE(21,1010) LATTIC,NAT,cform,IREL
      fac=180.d0/dpi
!     WRITE IN LATTICE CONSTANTS
      WRITE(21,1020) Alat,ALPHA(1),ALPHA(2),ALPHA(3)
      INDEX=0                                                           
      DO JATOM = 1,NAT
        INDEX=INDEX+1                                                  
        do J=1,3
           if(pos(j,index) .gt. 0.999999995D0) pos(j,index)=pos(j,index)-1.D0
        enddo
        WRITE(21,1030) IATNR(JATOM),(POS(J,INDEX),J=1,3 ),MULT(JATOM) &
          ,ISPLIT(JATOM)
        DO M=1,MULT(JATOM)-1                                     
          INDEX=INDEX+1                                            
        do J=1,3
           if(pos(j,index) .gt. 0.999999995D0) pos(j,index)=pos(j,index)-1.D0
        enddo
          WRITE(21,1031)IATNR(JATOM),(POS(J,INDEX),J=1,3)
        ENDDO                                                    
        WRITE(21,1049) ANAME(JATOM),JRI(JATOM),RNOT(JATOM),RMT(JATOM), &
          zatom(jatom)     
        WRITE(21,1055) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)
        ENDDO                                                          
      write(21,10) IORD                                                  
      DO 5 J=1,IORD                                                     
    5 write(21,11) ( (IZ(J1,J2,J),J1=1,3),TAU(J2,J), J2=1,3),J          
!.....READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS         
   10 FORMAT(I4,6X,'NUMBER OF SYMMETRY OPERATIONS')    
   11 FORMAT(3(3I2,F10.7/),i4)             
      STOP

 1000 FORMAT(A)
 1010 FORMAT(A4,23X,I3,A,/,13X,A4,14X,A4)
 1020 FORMAT(6F10.6,10X,F10.7)
 2030 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)
 2031 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8)
 1030 FORMAT('ATOM',I4,': X=',F10.8,1X,'Y=',F10.8,1X,'Z=',F10.8,/, &
       10X,'MULT=',I2,10X,'ISPLIT=',I2)
 1031 FORMAT(4X,I4,': X=',F10.8,1X,'Y=',F10.8,1X,'Z=',F10.8)
 1040 FORMAT(1X,'Q(U)  :',12(2X,F7.4) )
 1049 FORMAT(A10,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f10.5)      
 1055 FORMAT('LOCAL ROT MATRIX:',3X,3F10.7,2(/,20x,3F10.7))

      END
!

