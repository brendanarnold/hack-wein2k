      SUBROUTINE INILPW(NUMKPT,DEFFN,INFO,nmat_only)
!
      use loabc, only   : ELO
      use lapw_timer, only : READ_CPU_TIME, time_g1, time_atpar, &
                             START_TIMER, STOP_TIMER
      use parallel, only: MYID
      use lolog, only   : lapw
      use orb, only     : BEXT, VORB, NMOD, NSP, NATORB, IATOM, NLORB, LORB, init_orb
      use atspdt, only  : INS, NL, E
      use char, only    : LATTIC, NAME, init_char
      use comc, only    : IPGR, KNAME, init_comc
      use comi, only    : NT, LNSMAX, NAT
      use coml, only    : SPRSWF, PRNTWF, WFTAPE, REL, ITER, NOHNS
      use comr, only    : RKM, EL, EU, WEIGHT, init_comr
      use kpts, only    : SX, SY, SZ, init_kpts
      use potnlc, only  : JRJ, DH, RO, init_potnlc
      use rotmat, only  : ROTLOC, init_rotmat
      use struk, only   : POS, ALAT, ALPHA, RMT, V, VI, IATNR, MULT, multmax, init_struk
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            INFO, NUMKPT
      CHARACTER*(*)      DEFFN
!
!     ..................................................................
! 1.     PROGRAM UNIT 'INILPW'
!           Setup COMMON block data necessary for K-point processing.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Initialize the processing of K-points - setup data necessary
!           for the calculation of Hamilton and Overlap matrices and the
!           calculation of Eigenvalues and Eigenvectors).
!
! 3.     USAGE
!           INTEGER NUMKPT, INFO
!           CALL INILPW(NUMKPT,INFO)
!           CALL INIKPT(INFO)
!           DO I = 1, NUMKPT
!              CALL SETKPT(I)
!              CALL CALKPT
!              CALL PRTKPT(I,INFO)
!           ENDDO
!
!        ARGUMENT-DESCRIPTION
!           DEFFN  - CHARACTER*(*) string                        (input)
!                    name of the LAPW1 definition file 'lapw1.def'
!           INFO   - INTEGER value                              (output)
!                    returns the error condition:
!                    INFO .EQ. 0 ... no error
!                    INFO .NE. 0 ... currently not supported
!           NUMKPT - INTEGER value                              (output)
!                    returns the number of k-points to process.
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           ATPAR  - solve radial Schroedinger equation, calc. gaunts
!           CPUTIM - measure used CPU-time
!           GAUNT2 - initialize GAUNT routine
!           LATGEN - setup lattice and, if required, rotation matrices
!           NN     - checks if overlapping atomic spheres exist
!           OUTERR - output a message on the standard error device
!           RDSWAR - initialize warpin calculation
!
!        INDIRECTLY CALLED SUBROUTINES
!           ABC    - calculate local orbital coefficients A,B and C
!           CBCOMB - combine the radial total potential coefficients
!           GAUNT1 - calculate gaunt factors
!           OUTWIN - solve an IVP ODE using a Runge-Kutta method
!           RINT13 - perform radial integrals
!           ROTDEF - find symmetry operations, define rotation matrices
!           SELECT - automatic energy search
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           CALKPT - calculate Eigenvalues/-vectors for one k-point
!           INIKPT - setup values (valid for all k-points) not covered
!                    by 'INILPW' (e.g. warpin)
!           SETKPT - set current k-point
!           PRTKPT - print out/save the results obtained by 'CALKPT'
!
!        INPUT/OUTPUT (READ/WRITE)
!        Input:
!           the following files are read by 'INILPW':
!              - a file containing all files to open
!                (given as a command-line parameter)
!              - structure-file of the crystal
!              - input-file for LAPW1
!        Output:
!           'INILPW' writes to the following files:
!              - output-file for LAPW1
!              - vector-file (energy parameters for the first k-point
!                are written)
!        Error messages (to an error file):
!           - 'INILPW - can''t open definition file ' FNAME
!             where FNAME is the name of the file (given as command-line
!             parameter) which contains all files to be opened.
!           - 'INILPW - can''t open unit: ' IUNIT
!             'INILPW -         filename: ' FNAME
!             'INILPW -           status: ' STATUS '  form: ' FORM
!             where IUNIT is the intended unitnumber for the file,
!                   FNAME is the filename of the file to be opened,
!                   STATUS is the status of the file to be opened,
!                   FORM is the formatting parameter of the file.
!           - 'INILPW - MULT .EQ. 0'
!             The number of equivalent atoms has to be greater 0.
!           - 'INILPW - Maximum number of K-points exceeded.'
!             The number of K-points to be processed must be less
!             than or equal to NKPT.
!           - 'INILPW - Too many atoms (NATO too small)'
!             The number of inequivalent atoms must be less than or
!             equal to NATO.
!           - 'INILPW - Error reading file: ' FNAME
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           - COMPLEX*16 declaration for GFAC is used
!
! 4.     REMARKS
!           none
!
! 5.     METHOD
!           - setup timing
!           - open all required files
!           - initialize gaunt routine
!           - read lattice parameters
!           - read crystal structure
!           - determine output-format parameters
!           - read convergence parameters
!           - setup lattice and, if required, rotation matrices
!           - solve radial Schroedinger equation, calculate gaunts
!             and radial integrals of nonspherical potential
!           - initialize warpin calculations
!           - read first K-point
!           - write energy parameters on vector-file for the
!             first K-point
!           - read the remaining K-points
!           - perform timing
!           - return the number of read K-points to the calling
!             routine
!
!        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!     ..................................................................
!
       character*25 charop(3)
       character*5  charsp (-1:1)
!
!        Local Scalars
!
      INTEGER            KINDEX, ITAPE, IUNIT, INDEX, JATOM, J
      INTEGER            I, I1, I2, NT1, IA, M, IRECL,ios
      INTEGER            ISX, ISY, ISZ, IDV, JAT, L, M1, K
      DOUBLE PRECISION   E1N, E2N
      DOUBLE PRECISION, allocatable :: ZZ(:)
      DOUBLE PRECISION   PI, E1, E2, RNN
      CHARACTER*4        IREL
      CHARACTER*5        MODUS
      CHARACTER*11       STATUS, FORM
      CHARACTER*67       ERRMSG
      CHARACTER*80       FNAME, TITLE
      LOGICAL            RUNORB,nmat_only
      LOGICAL            FORCE,newform
      INTEGER            NK1
!
!        External Subroutines
!
      EXTERNAL           OUTERR, GAUNT2, LATGEN, NN, ATPAR
      EXTERNAL           RDSWAR
!
!        Intrinsic Functions
!
      INTRINSIC          LOG, EXP, ATAN
!
!        initialize timing for 'INILPW'
!
      CALL START_TIMER(time_g1)
!
      PI = 4.0D+0*ATAN(1.0D+0)
!
!
!        open all files listed in 'lapw1.def'
!
      nohns=.false.
      iter=.false.
      force=.false.
      runorb=.false.
      nmat_only=.false.
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         if(iunit.eq.72) nmat_only=.true.
         if(iunit.eq.7) runorb=.true.
         if(iunit.eq.71) then
!	   close(71)
           force=.true.
         end if
         if(iunit.eq.6.and.myid.ne.0) then
!	    write(*,*) 'iunit=',iunit,',  myid=',myid
            close(6)
!            close(96)
            write (fname,'(a5,i2.2)') deffn(5:9),myid
            fname="stdout." // fname
            open(6,file=fname,status=status,form=form,err=920)
         end if
         if(iunit.eq.97)   then
           close(97)
           nohns=.true.
         end if
         if(iunit.eq.98) then
           iter=.true.
         endif
      GOTO 10
   20 CONTINUE
      CLOSE (1)
!
!        calculate initialization data
!
      NT = LMAX
!
!        initialize GAUNT-routine
!
      CALL GAUNT2
!
      EU = -1.0D+0
      EL =  0.0D+0
!
      READ (20,5000) TITLE
      IF(myid.eq.0) WRITE (6,5000) TITLE
      READ (20,5010) LATTIC, NAT, IREL
      call init_orb(NAT)
      call init_char(NAT)
! read the orbital potential
       if(runorb) then
         charop(1)=' LDA+U potential'
         charop(2)=' Orbital polarization'
         charop(3)=' Interaction with Bext'
         charsp(-1)=' down'
         charsp(0)= ' para'
         charsp(1)= ' up  '
         nmod=0
         read(7,*,end=555,err=555)nmod,nsp,natorb,Bext
         DO JAT = 1, NATORB
          read(7,*)iatom(jat),nlorb(jat)
           do nl=1,nlorb(jat)
            read(7,*)lorb(nl,jat)
            l=lorb(nl,jat)
            IF(myid.eq.0) write(6,556)charop(nmod),iatom(jat),lorb(nl,jat),charsp(nsp)
            IF(myid.eq.0) write(21,556)charop(nmod),iatom(jat),lorb(nl,jat),charsp(nsp)
 556   format(a22,' added for atom type',i3,' L=',i3,' spin',a5)
              do i=-l,l
        	do j=-l,l
        	  read(7,103)vorb(jat,l,i,j)
  103   format(2f18.11)
                END DO
               enddo
 102   format(' M=',i3,7f11.5)
             IF(myid.eq.0) write(6,*)
             IF(myid.eq.0) write(6,*)' Orbital potential real part'
               do m=-l,l
        	IF(myid.eq.0) write(6,102)m,(dble(vorb(jat,l,m,m1)),m1=-l,l)
               enddo
             IF(myid.eq.0) write(6,*)
             IF(myid.eq.0) write(6,*)' Orbital potential imaginary part'
               do m=-l,l
        	IF(myid.eq.0) write(6,102)m,(dimag(vorb(jat,l,m,m1)),m1=-l,l)
               enddo
             IF(myid.eq.0) write(6,*)
          END DO
        END DO
      end if
555   continue
!
!        read in lattice constants
!
      READ (20,5020) ALAT(1), ALAT(2), ALAT(3), &
                     ALPHA(1), ALPHA(2), ALPHA(3)
      IF (ALPHA(1) .EQ. 0) ALPHA(1) = 90.0D+0
      IF (ALPHA(2) .EQ. 0) ALPHA(2) = 90.0D+0
      IF (ALPHA(3) .EQ. 0) ALPHA(3) = 90.0D+0
      ALPHA(1) = ALPHA(1)*PI/180.0D0
      ALPHA(2) = ALPHA(2)*PI/180.0D0
      ALPHA(3) = ALPHA(3)*PI/180.0D0
      IF (IREL .EQ. 'RELA') REL = .TRUE.
      IF (IREL .EQ. 'NREL') REL = .FALSE.
!
!        read crystal-structure (atompositions, symmetry-parameters,
!                                muffin-tin radius, ...)
!        'INDEX' counts all atoms in the unit cell,
!        'JATOM' counts only the notequivalent atoms
!
      call init_struk(NAT, 48*NAT, 0)
      call init_potnlc(NAT, NRAD)
      call init_rotmat(48*NAT, NAT, 0)
      allocate( ZZ(NAT) )
      INDEX = 0
      multmax=0
      DO 50 JATOM = 1,NAT
         INDEX = INDEX + 1
!
         READ (20,5030) IATNR(JATOM),(POS(J,INDEX),J=1,3),MULT(JATOM)
         IF (MULT(JATOM) .EQ. 0) THEN
!
!           illegal number of equivalent atoms
!
            WRITE (6,6000) JATOM, INDEX, MULT(JATOM)
            GOTO 930
         ENDIF
         if(multmax.lt.mult(jatom)) multmax=mult(jatom)
         DO 40 M = 1, MULT(JATOM) - 1
            INDEX = INDEX + 1
            READ (20,5040) IATNR(JATOM), (POS(J,INDEX),J=1,3)
   40    CONTINUE
         READ (20,5050) NAME(JATOM),JRJ(JATOM),RO(JATOM),RMT(JATOM) &
                       ,ZZ(JATOM)
         DH(JATOM)  = LOG(RMT(JATOM)/RO(JATOM))/(JRJ(JATOM) - 1)
         RMT(JATOM) = RO(JATOM)*EXP(DH(JATOM)*(JRJ(JATOM)-1))
         READ (20,5060) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)
   50 CONTINUE
      call init_struk(NAT, INDEX, 1)
      call init_rotmat(INDEX, NAT, 1)
!
!        read output-format parameters
!
      READ (5,5070) MODUS
      IF (MODUS .EQ. 'WFPRI') THEN
         PRNTWF = .TRUE.
         SPRSWF = .FALSE.
         WFTAPE = .TRUE.
      ELSEIF (MODUS .EQ. 'WFFIL') THEN
         PRNTWF = .FALSE.
         SPRSWF = .FALSE.
         WFTAPE = .TRUE.
      ELSEIF (MODUS .EQ. 'ENFIL') THEN
         PRNTWF = .FALSE.
         SPRSWF = .TRUE.
         WFTAPE = .FALSE.
      ELSE
         PRNTWF = .FALSE.
         SPRSWF = .TRUE.
         WFTAPE = .FALSE.
      ENDIF
      IF(myid.eq.0) WRITE (6,6010) LATTIC,IREL
!
!        read in R-MT times K-MAX (used to check convergence)
!
      READ (5,*) RKM, NT1, LNSMAX
      IF (LNSMAX .EQ. 0) LNSMAX = 2
      IF (LNSMAX .GE. NSLMAX) LNSMAX = NSLMAX - 1
      IF (NT1 .GT. 0)    NT = NT1 + 1
      IF (NT  .GT. LMAX) NT = LMAX
      IF(myid.eq.0) WRITE (6,6020) RKM, NT - 1, LNSMAX, NAT
!
!        set up lattice, and if required rotation matrices
!
      CALL LATGEN(NAT)
!
!        check if atomic spheres are non overlapping
!
      CALL NN(NAT)
!
!        determine minimum RMT
!
      RNN = RMT(1)
      DO 60 IA = 2, NAT
         IF (RMT(IA) .LT. RNN) RNN = RMT(IA)
 60   CONTINUE
!
!        solve radial Schroedinger-equation, calculate gaunts and
!        radial integrals UVU, UVUE of nonspherical potential
!
      CALL START_TIMER(time_atpar)
      CALL ATPAR(NT,REL,NAT,LNSMAX,ZZ,force)
      CALL STOP_TIMER(time_atpar)
      IF(myid.eq.0) WRITE (6,6030) READ_CPU_TIME(time_atpar)
!
!        define spherevolume/cellvolume
!
      DO 70 J = 1, NAT
         V(J) = 4.0D+0*PI*VI*RMT(J)**3/3.0D+0
 70   CONTINUE
      IF(myid.eq.0) then
      WRITE (6,6040) (RMT(J), J = 1, NAT)
      WRITE (6,6050) (V(J),   J = 1, NAT)
      WRITE (6,6060) VI
      WRITE (6,6070) ALAT
      endif
!
!        set up array for warpin
!
      IF (INS .EQ. 0) CALL RDSWAR
!
!        set up list of plane waves for hamilton-matrix
!
      RKM = RKM/RNN
      IF (NT .GT. LMAX) NT = LMAX
!
!        read in K-points and save them in an array
!
      READ (5,5090,err=71) ITAPE,e1n,e2n
!
!        read first K-point
!
 71   continue
      NK1=NKPTSTART
      call init_comc(NK1, 0)
      call init_comr(NK1, 0)
      call init_kpts(NK1, 0)
      KINDEX = 1      
      READ (ITAPE,5101,IOSTAT=ios) KNAME(KINDEX), ISX, ISY, ISZ, IDV, &
                        WEIGHT(KINDEX), E1, E2, IPGR(KINDEX)
      IF(ios==0) THEN
         newform=.TRUE.
      ELSE
         REWIND(ITAPE)
         READ (ITAPE,5100,IOSTAT=ios) KNAME(KINDEX), ISX, ISY, ISZ, IDV, &
                        WEIGHT(KINDEX), E1, E2, IPGR(KINDEX)
         newform=.FALSE.
      ENDIF
      IF (KNAME(KINDEX) .EQ. 'END       ') GOTO 510
      SX(KINDEX) = DBLE(ISX)/DBLE(IDV)
      SY(KINDEX) = DBLE(ISY)/DBLE(IDV)
      SZ(KINDEX) = DBLE(ISZ)/DBLE(IDV)
!
!        set up energy range limits for 'SECLR4'
!
      IF (E1 .NE. E2) THEN
         EL=E1
         EU=E2
      ENDIF
      IF (E1n .NE. E2n) THEN
         EL=E1n
         EU=E2n
      ENDIF
!
!        write energy parameter on vector-file for first K-point only
!
      IF(myid.eq.0) then
      DO 80 I = 1, NAT
         do j=1,lmax
            if(.not.lapw(j-1,i)) e(j,i)=e(j,i)+200.D+0
         enddo
         WRITE(10) (E(J,I),J=1,LMAX)
         WRITE(10) ((ELO(J,k,I),J=0,LOMAX),k=1,nloat)
         WRITE(11,'(100(f9.5))') (E(J,I),J=1,LMAX)
         WRITE(11,'(100(f9.5))') ((ELO(J,k,I),J=0,LOMAX),k=1,nloat)
         do j=1,lmax
            if(.not.lapw(j-1,i)) e(j,i)=e(j,i)-200.D+0
         enddo
   80 CONTINUE
      endif

      DO 
        KINDEX = KINDEX + 1
        IF(KINDEX.GE.NK1) THEN
          NK1 = NK1 + NK1 / 2
          call init_comc(NK1,1)
          call init_comr(NK1,1)
          call init_kpts(NK1,1)
        END IF
!
!        read in K-point to be calculated and energy range
!
        IF(newform) THEN
           READ (ITAPE,5101) KNAME(KINDEX), ISX, ISY, ISZ, IDV, &
                WEIGHT(KINDEX), E1, E2, IPGR(KINDEX)
        ELSE
           READ (ITAPE,5100) KNAME(KINDEX), ISX, ISY, ISZ, IDV, &
                WEIGHT(KINDEX), E1, E2, IPGR(KINDEX)
        ENDIF
         IF (KNAME(KINDEX) .EQ. 'END       ') GOTO 510
         SX(KINDEX) = DBLE(ISX)/DBLE(IDV)
         SY(KINDEX) = DBLE(ISY)/DBLE(IDV)
         SZ(KINDEX) = DBLE(ISZ)/DBLE(IDV)
      END DO
!
!        too many k-points
!
      GOTO 940
  510 CONTINUE
!
!        successful exit
!
      NUMKPT = KINDEX - 1
      INFO = 0
      CALL STOP_TIMER(time_g1)
      GOTO 998
!
!        error handling
!
  910 INFO = 1
!
!        'lapw1.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('INILPW',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('INILPW',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('INILPW',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('INILPW',ERRMSG)
      GOTO 999
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      WRITE (ERRMSG,9050) JATOM, MULT(JATOM), INDEX
      CALL OUTERR('INILPW',ERRMSG)
      CALL OUTERR('INILPW','MULT .EQ. 0')
      GOTO 999
  940 INFO = 4
      CALL OUTERR('INILPW','Maximum number of K-points exceeded.')
      GOTO 999
  950 INFO = 5
      CALL OUTERR('INILPW','Too many atoms (NATO too small).')
      GOTO 999
  955 INFO = 6
      CALL OUTERR('INILPW','Too many atoms (NDIF too small).')
      GOTO 999
  960 INFO = 7
!
!        Error reading file 'lapw1.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('INILPW',ERRMSG)
      GOTO 999
!
  998 CONTINUE
      deallocate( ZZ )
  999 CONTINUE
      RETURN
!
 5000 FORMAT(A80)
 5010 FORMAT(A4,23X,I3,/,13X,A4)
 5020 FORMAT(6F10.7)
 5030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
 5040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
 5050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 5060 FORMAT(20X,3F10.8)
 5070 FORMAT(A5)
 5080 FORMAT(F10.7,2I5)
 5090 FORMAT(20X,I1,2f10.1)
 5100 FORMAT(A10,4I5,3F5.2,A3)
 5101 FORMAT(A10,4I10,3F5.2,A3)
 6000 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
 6010 FORMAT(1H0,9X,A4,3X,'TYPE LATTICE ASSUMED'/10X,A4, &
              '-CALCULATION',/)
 6020 FORMAT(10X,'R-MT TIMES K-MAX IS',F5.2,/10X,'MAX L IS',I3,5X, &
      'MAX L IN NONSPHERICAL MATRIXELEMENTS:',I3,/ &
             ' NUMBER OF ATOMS IS',I5)
 6030 FORMAT(40X,'CPTIME ATPAR  =',F8.1)
 6040 FORMAT(///10X,'R-MT=',4F10.7/15X,4F10.7/15X,4F10.7)
 6050 FORMAT(10X,' FRACTIONAL VOLUME WITHIN MT=',4F10.7,2(/39X,4F10.7))
 6060 FORMAT(10X,' ONE/UNIT CELL VOLUME=',E16.9)
 6070 FORMAT(/10X,'LATTICE CONSTANTS ARE:',3F9.5)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
 9050 FORMAT('MULT(',I3,'=',I3,', INDEX=',I3)
!
!        End of 'INILPW'
!
      END
