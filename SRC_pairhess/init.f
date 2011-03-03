      subroutine init(pred,pos,br1,br2,n1,n2)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
!     Read the structure file
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pos(3,48*natmax),pred(3,natmax)
      dimension inum(48),br1(3,3),br2(3,3),shift(3)
      character*1 junk
      character*80 title1
      dimension alpha(3)
      character*6 lattic,cform
      logical doshift
!.....READ STRUCT
      READ(20,1000) TITLE1
      READ(20,1010) LATTIC,NAT,cform,IREL
!
!     Shifts for CXY, CYZ, CXZ lattices
      doshift=.false.
      IF(LATTIC(1:3) .eq. 'CXY')then
        shift(1)=0.5D0
        shift(2)=0.5D0
        shift(3)=0.0D0
        doshift=.true.
      ELSE IF(LATTIC(1:3) .eq. 'CYZ')then
        shift(1)=0.0D0
        shift(2)=0.5D0
        shift(3)=0.5D0
        doshift=.true.
      ELSE IF(LATTIC(1:3) .eq. 'CXZ')then
        shift(1)=0.5D0
        shift(2)=0.0D0
        shift(3)=0.5D0
        doshift=.true.
      ENDIF
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
        READ(20,2030) IATNR,( POS(J,INDEX),J=1,3 ),MULT(JATOM),  &
          ISPLIT
          DO J=1,3
                if(POS(J,INDEX) .lt.     -TOL)POS(J,INDEX)=POS(J,INDEX)+1.D0
                if(POS(J,INDEX) .gt. 1.D0-TOL)POS(J,INDEX)=POS(J,INDEX)-1.D0
                PRED(J,JATOM)=POS(J,INDEX)
          ENDDO
!       Add for centered cells
        if(doshift)then
                index=index+1
                DO J=1,3
                POS(J,INDEX)=POS(J,INDEX-1)+shift(J)
                ENDDO
                if(POS(1,index) .gt. 0.99999D0)pos(1,index)=pos(1,index)-1
                if(POS(2,index) .gt. 0.99999D0)pos(2,index)=pos(2,index)-1
                if(POS(3,index) .gt. 0.99999D0)pos(3,index)=pos(3,index)-1
        endif
        DO 55 M=1,MULT(JATOM)-1
          INDEX=INDEX+1
          READ(20,2031) IATNR,( POS(J,INDEX),J=1,3)
          DO J=1,3
                if(POS(J,INDEX) .lt.     -TOL)POS(J,INDEX)=POS(J,INDEX)+1.D0
                if(POS(J,INDEX) .gt. 1.D0-TOL)POS(J,INDEX)=POS(J,INDEX)-1.D0
          ENDDO
!         Add for centered cells
          if(doshift)then
                index=index+1
                DO J=1,3
                POS(J,INDEX)=POS(J,INDEX-1)+shift(J)
                ENDDO
                if(POS(1,index) .gt. 0.99999D0)pos(1,index)=pos(1,index)-1
                if(POS(2,index) .gt. 0.99999D0)pos(2,index)=pos(2,index)-1
                if(POS(3,index) .gt. 0.99999D0)pos(3,index)=pos(3,index)-1
          endif
   55   CONTINUE
        if(doshift)MULT(Jatom)=MULT(Jatom)*2
          DO J=1,4
             READ(20,1000)junk
          ENDDO
   50 CONTINUE
!
      N1=nat
      N2=index
!     Get symmetry
      READ(20,1151) iord
1151  FORMAT(I4)

      DO j=1,iord
         READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
1101     FORMAT(3(3I2,F10.8/),I8)
      ENDDO
!     Average, forcing within -0.5 to 0.5
      call fixup3(pred)
!     Generate the lattice
      call gen_brav(lattic,alpha,br1,br2)

      return
 1000 FORMAT(A)
 1010 FORMAT(A4,23X,I3,A20,/,13X,A4,14X,A4)
 1020 FORMAT(6F10.6,10X,F10.7)
 2030 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)
 2031 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8)

      END
!
