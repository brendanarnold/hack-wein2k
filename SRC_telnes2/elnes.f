!BOP
! !ROUTINE: TELNES2
! !INTERFACE:
      PROGRAM TELNES2
! !USES:
	  use program_control,only : averaged,verbosity,modus,fileroot
	  use crash
! !DESCRIPTION:
!     TELNES2 calculates Electron Energy Loss Near Edge Structure.
!     (according to setting in *.innes)
!     This is the main program.
!     It opens files and calls several subroutines to do the calculation.
! !REVISION HISTORY:
! *** HISTORY OF TELNES2 ***
!     TELNES was based on TXSPEC (C)1996 by Joachim Luitz
!
!     first version :
!     TELNES (C)1997-1998 by Pierre-Henri Louf
!     Inst. f. Angew. u. Techn. Physik, TU Wien
!     Wiedner Hauptstr. 8-10, A-1140 Wien, AUSTRIA
!
!     1998,1999 TELNES modified by Michael Nelhiebel
!     implemented into WIEN package 1999 by Joachim Luitz
!
!     2004,2005 TELNES2 program, major cleanup and new functionality by Kevin Jorissen
!     (mailto kevin.jorissen@ua.ac.be for help and fan mail,
!                  or try the WIEN2k mailing list)
!EOP

      implicit none
!     LOCAL STAFF
!     Next variables are needed for opening files :
      INTEGER IUNIT
      CHARACTER *11 STATUS,FORM
      CHARACTER *80 DEFFN,FNAME


!     we read the names or our files and open them:
! windows fix :
!      ERRFN(1:12)='telnes2.error'
!      errfn(12:80)=' '
!      CALL ERRFLG(ERRFN,'Error in TELNES2')
!      OPEN (1,FILE='telnes2.def',STATUS='OLD',ERR=910)
! end of windows fix

      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in TELNES2')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)

 8003 CONTINUE
      READ (1,*,END=8001,ERR=911) IUNIT,FNAME,STATUS,FORM
      OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=912)
         READ (1,*,END=8001,ERR=911) IUNIT,FNAME,STATUS,FORM
		 if(iunit.eq.5) fileroot=fname
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=912)
      GOTO 8003
 8001 CONTINUE


      write(6,'(/,a)') 'Output from main routine :'
      write(6,'(a)') '#########################   ELNES CALCULATION ########################'
      write(6,'(a,//)') 'All files opened.  TELNES2 is ready to go!'

!     Read case.innes (principal input file) and describe the task.
      call errflg(errfn,'Error in reading case.innes')
      call ReadInnes
      
!     Read case.struct and prepare structural information.
      call errflg(errfn,'Error in readstruct')
      call ReadStruct

!     calculate the coefficients for the integral over collection and convergence angles.
      call errflg(errfn,'Error in CalculateWeights')
      call CalculateWeights

!     obtain initial state wave functions.
      call errflg(errfn,'Error while obtaining initial state wave functions.')
      call Corewavefunction

!     get DOS or cross DOS.
      call errflg(errfn,'Error in DensityOfStates')
      call DensityOvStates 


!     continue initialization of energy scale.  (started while processing DOS)
      call errflg(errfn,'Error in EnergyMesh')
      call EnergyMesh

!     calculate the radial functions for the final states.
      call errflg(errfn,'Error in the calcul of radial functions')
      call RadialFunctions
      if(verbosity.ge.2) then
!       calculate orthogonality integrales and print them in case.ortho .
        call errflg(errfn,'Err. in the calcul of orthogonality integrals')
        call Orthogonality
      endif

!     calculate the rotation matrices for the Q-vectors.
      call errflg(errfn,'Error in the calcul of the rotation matrices')
      call CalculateMatrices

      if (.not.averaged) then
!        calculate norm of radial functions for the final states.
         call errflg(errfn,'Error in the calcul of the Norme integrals')
         call NormOfRadialFunctions
      endif

      call errflg(errfn,'Main calculation')
!     calculate ELNES spectra.
      if (modus.eq.'E') then
         call CalculateEnergySpectrum
      elseif(modus.eq.'A') then
         call CalculateAngularSpectrum
	  endif
          call errflg(errfn,'End of main calculation')

      write(6,'(/,a)') 'Output from main routine :'
      write(6,'(a)') "Calculation probably finished succesfully :-)."
      write(6,'(a)') 'Main routine says goodbye!'
      call errclr(errfn)
      STOP 'TELNES2 DONE'


!     TELNES.def couldn't be opened
 910  WRITE(ERRMSG,9000) FNAME
      CALL OUTERR('TELNES2',ERRMSG)
      GOTO 999

!     TELNES.def is corrupt
 911  WRITE(ERRMSG,9001)FNAME
      CALL OUTERR('TELNES2',ERRMSG)
      GOTO 999

!     Error opening a unit
 912  WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('TELNES2',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('TELNES2',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('TELNES2',ERRMSG)
      GOTO 999

 999  STOP 'TELNES2 - Error'

 9000 FORMAT('can''t open definition file ',A40)
 9001 FORMAT('can''t read definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
       
      END


