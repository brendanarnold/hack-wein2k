      PROGRAM TXSPEC
      use atomgrid
      use reallocate
!
!     X-ray Spectra
!
!     TXSPEC calculates theoretical x-ray spectra
!     (either emission or absorption spectrum, 
!     according to setting in *.inxs)
!
!     for these we use:
!       emission spectra:   DOS(E.GT.EF)=0 (EMIS=.TRUE.)
!       absorption spectra: DOS(E.LT.EF)=0 (EMIS=.FALSE.)
!
!     
!     (C)1996 by Joachim Luitz
!
!
!     Splitting of L2/L3 (M5/M7) was implemented by
!     G""unter Schuster and Alexander Tomandl (1997)
! 

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'

      EXTERNAL OUTERR

!     COMMON /DIRA/ aus hfsd.f
      COMMON /DIRA/ DV(NRAD+41),DR(NRAD+41),DP(NRAD+41),DQ(NRAD+41), &
                    DPAS,Z,NSTOP,NES,TETS,NP,NUC
      dimension DPT(NRAD+41),DQT(NRAD+41)

    
      DOUBLE PRECISION   A(NRAD)
      DOUBLE PRECISION   AE(NRAD)
      DOUBLE PRECISION   AP(NRAD), B(NRAD)
      DOUBLE PRECISION   BE(NRAD)
      DOUBLE PRECISION   BP(NRAD)
      COMMON  /WORK/     A,B,AP,BP,AE,BE,VLM          
      SAVE    /WORK/
       
      real*8,allocatable :: A1(:,:), B1(:,:)
      real*8,pointer :: ENE(:)
      double precision EF
      COMMON  /ENE/      EF
      SAVE    /ENE/

      LOGICAL REL
!     Local Variables:
     
      INTEGER NATOM,NC,LC,KAPPA
      INTEGER IUNIT
      
      CHARACTER *11 STATUS,FORM
      CHARACTER*80  DEFFN, ERRFN
      CHARACTER *80 FNAME
      CHARACTER *67 ERRMSG
      CHARACTER *4  DUMMY,EMISS

!     DOS() Density of States 
!     XI()  Xray Intensity

      real*8,pointer :: DOS(:,:),XI(:,:),X(:)
      real*8,pointer :: XINTER(:),XOUT(:)
!     S Matrix-Element
      DOUBLE PRECISION S
      DOUBLE PRECISION NY

      LOGICAL TEST1, TEST2
!     EMIS: are we calculating emission spectra?
      LOGICAL EMIS
      LOGICAL INC
      LOGICAL BANDRA
      BANDRA=.FALSE.

!     we read the names or our files and open them:
      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in TXSPEC')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
 8003 CONTINUE
         READ (1,*,END=8001,ERR=911) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=912)
      GOTO 8003
 8001 CONTINUE

!     From our *.inxs we get the crucial data,
!     READ title into FNAME  and forget it
      READ(5,'(A)')FNAME
      READ(5,*)NATOM
      READ(5,*)NC
      READ(5,*)LC
      READ(5,*)SPLIT,XINT1,XINT2
    
!
!     READ(5,*)Kappa
!     just take the negative kappa

      kappa=-1-LC

!
!
      READ(5,*)EMIN,DE,EMAX
      READ(5,'(4a)')DUMMY
      IF(DUMMY.EQ.'EMIS')then
        EMIS=.TRUE.
	EMISS='EMIS'
        write(*,*)'Emis.true'
      ELSE
        EMIS=.FALSE.
	EMISS='ABS'
      ENDIF
      IF(EMIS) THEN
         READ(5,*)S
         READ(5,*)G0
         READ(5,*)W
         READ(5,'(4a)')DUMMY
         IF(DUMMY.NE.'AUTO')then
            BANDRA=.FALSE.
            write(*,*)'BANDRA.false'
         ELSE
            BANDRA=.TRUE.
         ENDIF
      endif   
!
!     We calculate our core radials using the hfsd-Sub.

      call hfsd(NATOM,NC,LC,KAPPA,REL)

!     returning from hfsd check if NP was set to zero:
!     if NP.EQ.0 then requested conditions (n,l,kappa) were not matched
!
      IF (NP.EQ.0) THEN
!        WRITE (*,*)'Error in n,l or kappa!'
         GOTO 919
      ENDIF
!
!     hfsd(n,l,kappa) returns radial function: r*P(core,n',l')
!     in DP(*) and DQ(*); 


!     now that REL is set (from HFSD): initialize constants
!     CFEIN1 and CFEIN2 are used in RINT13 to calculate the
!     integrals of radialfunctions and are derivates of the
!     Feinstruktur-constant (in Hartree)
!
      IF (REL) THEN
         CFEIN1 = 1.0D+0
! we have to play this trick because DQ was already modified
! with CFEIN2 in atpar!
         CFEIN2 = SQRT(1.d0/137.0359895d0**2)
!         CFEIN2 = 4.0D0*1.331258D-5
      ELSE
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.0D-22
      ENDIF      

!     Read partial density of states from *.dos1ev:
!     IEMAX, ENE() DOS(i,1) and, if necessary, DOS(i,2)

      READ (32,4712) IEMAX
      allocate (   A1(NRAD,IEMAX), B1(NRAD,IEMAX))
      allocate ( ENE(IEMAX))
      allocate ( DOS(IEMAX,3),XI(IEMAX,2),X(IEMAX))
      allocate ( XINTER(IEMAX),XOUT(IEMAX))
!     to be on the save side, we clear XI()
      DO I=1,IEMAX
        XI(i,1)=0
        XI(i,2)=0
      ENDDO
      I=0
 1    CONTINUE
        i=i+1
        IF (LC.EQ.0) then
          READ(32,4713,END=2,ERR=913) ENE(i),DOS(i,1),DOS(i,3)
        ELSE
          READ(32,4713,END=2,ERR=913) ENE(i),DOS(i,1),DOS(i,2),DOS(i,3)
        ENDIF
      goto 1
 2    CONTINUE
      

!     if option search for Band Ranges was set to AUTO, then search
!     for Band ranges and append them to the case.inxs file

      IF(BANDRA) THEN
         n=1
         INC=.FALSE.
         DO i=1,IEMAX
            if(INC) then
               n=n+1
               inc=.FALSE.
            endif
            
            if (n.eq.1) then
               if (DOS(i,3).gt.0) then
                  E2=ENE(i)
                  inc=.TRUE.
               endif
            endif   
            if (n.eq.2) then
               if (DOS(i,3).eq.0) then
                  E1=ENE(i)
                  INC=.TRUE.
               endif
            endif            
            if (n.eq.3) then
               if (DOS(i,3).gt.0) then
                  E0=ENE(i)
                  INC=.TRUE.
               endif
            endif   
         ENDDO   

         if(n.eq.4) then
            write(6,*)"Parameters found: E0,E1,E2"
            write(6,*)E0
            write(6,*)E1
            write(6,*)E2
            write(5,*)E0
            write(5,*)E1
            write(5,*)E2
         else
            write(6,*)"Not all parameters found!"
            write(6,*)n,i
            write(6,*)"Substituting all Ei with E2"
!     thus set all Ei to E2
            write(5,*)E2
            write(5,*)E2
            write(5,*)E2
         endif   
      endif
!     read Fermi-Energy in Ry, used for conversion
      READ(30,4714)EF


!     Test if EMIN and EMAX are within the given
!     range of *.dos1ev; else set them appropriately

!      TEST1=.FALSE.
!      TEST2=.FALSE.
!      DO 3 J=1,IEMAX
!        IF ((EMIN.LE.ENE(J)).AND.(.NOT.TEST1)) THEN
!           JEMIN=J
!           EMIN=ENE(J)
!           TEST1=.TRUE.
!        ENDIF
!
!        IF ((EMAX.LE.ENE(J)).AND.(.NOT.TEST2)) THEN
!           JEMAX=J
!           EMAX=ENE(j)
!           TEST2=.TRUE.
!        ENDIF
!
!        IF (TEST1.AND.TEST2) THEN
!           GOTO 4
!        ENDIF
! 3    CONTINUE

!     if everything else fails:
!      IF (.NOT.TEST1) THEN 
         EMIN=ENE(1)
         JEMIN=1
!      ENDIF
!      IF (.NOT.TEST2) THEN 
         EMAX=ENE(IEMAX)
         JEMAX=IEMAX
!      ENDIF
! 4    CONTINUE

      WRITE(*,*) JEMIN,JEMAX
!     Initialize processing of K-points
      call INILPW(NAT,REL)


!     ...........................................................
!     Here comes the first part of the spectrum, using
!     LL=LC+1
!     ...........................................................

      LL=LC+1

      call ATPAR(NATOM,LL,NAT,REL,JEMIN,JEMAX,a1,b1,ene)

!     Wll': for ll=lc+1!
      W=WLL(LC,LL)

      write(*,*) EMISS, ' LC=',LC,' LL=',LL
      write(*,*) 'angular multiplication factor W= ',W

!     Write out core-wf
      do i=1,jrj(natom)
         WRITE(9,4713) DR(I),DP(I)/DR(I),DQ(I)/DR(I)
!         DPT(I)=DP(I)/DR(I)
!         DQT(I)=DQ(I)/DR(I)
      enddo   
!

!      CALL RINT13(CFEIN1,CFEIN2,DPT,DQT,DPT,DQT,S,NATOM)
!      write(*,*) 'TEST......',S


      DO J=JEMIN,JEMAX
         
         
         CALL RINT13(CFEIN1,CFEIN2,A1(1,J),B1(1,J),DP,DQ,S,NATOM)

!      write(*,*) 'TEST......',S
!     ..............................................................
!     The following lines were used for debugging purposes;
!     uncomment and recompile to write out data generated.
!     Note: uncommented use xspec.def.debug as template!
!     ..............................................................
!         IF(J.EQ.jemin) THEN
!           do i=1,jrj(natom)
!           write(48,*) dr(i),a1(i,j),b1(i,j),dp(i)/dr(i),dq(i)/dr(i)
!           enddo
!         end if
!         IF(J.EQ.jemax) THEN
!           do i=1,jrj(natom)
!           write(49,*) dr(i),a1(i,j),b1(i,j),dp(i)/dr(i),dq(i)/dr(i)
!           enddo
!         end if
!     ...............................................................


         EI=ENE(J)

!     write out Matrix-Element
         WRITE(53,*)EI,(S*S)
         
!        EMISSION SPECTRA
         IF(EMIS)then
            IF(EI.GT.0)then
               XI(J,1)=0
            ELSE
               XI(J,1)=W*S*S*DOS(J,1)
            ENDIF
            
!     ABSORPTION SPECTRA
         ELSE
            IF(EI.LE.0)then
               XI(J,1)=0
            else
               XI(J,1)=W*S*S*DOS(J,1)
            ENDIF
         ENDIF
         
      ENDDO

!     ...........................................................
!     Here comes the second part of the spectrum, using
!     LL=LC-1, but only if necessary
!     ...........................................................

       IF(LC.NE.0) THEN
          LL=LC-1
          
!        Wll': for ll=lc-1
          W=WLL(LC,LL)

          write(*,*) EMISS, ' LC=',LC,' LL=',LL
          write(*,*) 'angular multiplication factor W= ',W


         call ATPAR(NATOM,LL,NAT,REL,JEMIN,JEMAX,a1,b1,ene)

         DO J=JEMIN,JEMAX
            CALL RINT13(CFEIN1,CFEIN2,A1(1,J),B1(1,J),DP,DQ,S,NATOM)
            EI=ENE(J)
            
!     here comes the 2nd Matrix-Element
            WRITE(54,*)EI,(S*S)
            
!          EMISSION SPECTRA
            IF(EMIS)then
               IF(EI.GT.0)then
                  XI(J,2)=0
               ELSE
                  XI(J,2)=W*S*S*DOS(J,2)
               ENDIF
!     ABSORPTION SPECTRA
            ELSE
               IF(EI.LE.0)then
                  XI(J,2)=0
               else
                  XI(J,2)=W*S*S*DOS(J,2)
               ENDIF
            ENDIF
!     ..............................................................
!     The following lines were used for debugging purposes;
!     uncomment and recompile to write out data generated.
!     Note: uncommented use xspec.def.debug as template!
!     ..............................................................
!     IF(J.EQ.jemin) THEN
!     do i=1,jrj(natom)
!     write(55,*) dr(i),a1(i,j),b1(i,j),dp(i)/dr(i),dq(i)/dr(i)
!     enddo
!     end if
!     IF(J.EQ.jemax) THEN
!     do i=1,jrj(natom)
!     write(56,*) dr(i),a1(i,j),b1(i,j),dp(i)/dr(i),dq(i)/dr(i)
!     enddo
!     end if
!     ..............................................................

         ENDDO
      ENDIF
      
      DO J=JEMIN,JEMAX
         X(J)=XI(J,1)+XI(J,2)
!        WRITE(47,4713)ENE(J),X(J),XI(J,1),XI(J,2)
      ENDDO
      
!     Split-Berechnung    

      ISPLIT=ABS(NINT((JEMAX-JEMIN)/(EMAX-EMIN)*SPLIT))
      DE=(EMAX-EMIN)/(JEMAX-JEMIN)
!        print*, 'isplit',isplit,jemax,jemin,emax,emin,split,de
      if(isplit.ne.0) then
         write(6,*) 'Core level splitting of ',split,' eV included'
         write(6,*) 'with intensities ',xint1,xint2
      endif
      call doreallocate (ene,jemax+isplit) 
      call doreallocate (xout,jemax+isplit) 
      call doreallocate (xinter,jemax+isplit) 
      call doreallocate (xi,jemax+isplit,2) 
      DO I=1,jemax+isplit
         XINTER(I)=0.d0
      ENDDO
      
      DO I=JEMIN,JEMAX
         XINTER(I+ISPLIT)=X(I)
      ENDDO   

      IF (SPLIT.LT.0) THEN
         EMIN=EMIN+SPLIT
      ELSE
         EMAX=EMAX+SPLIT
      ENDIF

      DO I=jemin,jemax+isplit
         ENE(I)=EMIN+(I-1)*DE
      ENDDO      

!      DO I=JEMIN,JEMAX+ISPLIT
      DO I=JEMIN,JEMAX
         XOUT(I)=XINT1*X(I)+XINT2*XINTER(I)
         WRITE(*,*)i,ene(i),XOUT(I)
      ENDDO

   
! 53   DO J=JEMIN,JEMAX+ISPLIT
! restrict to original E range for abs, otherwise a part is missing
 53   DO J=JEMIN,JEMAX
         WRITE(47,4713)ENE(J),XOUT(J),XI(J,1),XI(J,2)
      ENDDO
      
      CALL ERRCLR(ERRFN)
      STOP 'TXSPEC DONE'


!     xspec.def couldn't be opened
 910  WRITE(ERRMSG,9000) FNAME
      CALL OUTERR('TXSPEC',ERRMSG)
      GOTO 999

!     xspec.def is corrupt
 911  WRITE(ERRMSG,9001)FNAME
      CALL OUTERR('TXSPEC',ERRMSG)
      GOTO 999

!     Error opening a unit
 912  WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('TXSPEC',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('TXSPEC',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('TXSPEC',ERRMSG)
      GOTO 999

!     Error reading from case.dos1ev
 913  WRITE(ERRMSG,9050)
      CALL OUTERR('TXSPEC',ERRMSG)
      GOTO 999

 919  WRITE(ERRMSG,9090)
      CALL OUTERR('TXSPEC',ERRMSG)
      GOTO 999

 999  STOP 'TXSPEC - Error'

 4712 FORMAT(/,38x,I5,/)
 4713 FORMAT(F10.5,e14.5,e14.5,e14.5)
 4714 FORMAT(//,56x,f10.5)
 9000 FORMAT('can''t open definition file ',A40)
 9001 FORMAT('can''t read definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
 9050 FORMAT('Error reading dos1ev file')
 9090 FORMAT('Error in n,l or kappa!')
      END

      DOUBLE PRECISION FUNCTION WLL(L,LP)
      INTEGER L,LP
!      DOUBLE PRECISION WLL

!     calculate W_ll' for transition L->LP

      IF (L.EQ.LP-1) THEN
        WLL=1.0d0*(L+1)/(2*LP+1)
      ELSE
        WLL=1.0d0*L/(2*LP+1)
      ENDIF

      RETURN
      END



