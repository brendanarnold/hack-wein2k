!BOP
! !ROUTINE: ReadInnes
! !INTERFACE:
      SUBROUTINE ReadInnes
! !USES:
	  use input
	  use program_control
	  use constants
	  use units
	  use leftovers,only :  GammaA0
      use energygrid,only : ef
	  use dimension_constants,only : lmax, nofcross, set_nofcross
! !DESCRIPTION:
!     Read from case.innes all the input.
!
!     If a new 'qtl' type qtl-file is to be used, the file 30 (case.qtl in elnes.def) will
!     be closed, and a file case.qtl\$natom is opened instead (which is supposed to exist
!     and be the output of P. Novak's qtl program run with suitable input options).
!     This had to be include in this routine, since it will try to read the Fermi energy
!     from file 30.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!     Changes June 2005 (Kevin Jorissen)
!EOP
      implicit none

!     LOCAL VARIABLES
      real*8 ,external :: ConversionDegRad,wavelength
	  character*4 switch
      logical okay,EfIsRead
	  integer i
      character*20 ,parameter :: testselrule = '0123456789mMdDqQoOnN'
      character*90 qtlfile   ! by choosing it a little larger than fileroot, we prevent accidents.
	  character h1,h2,h3,h4

!   Set defaults.
      EfIsRead=.false.
      call init_varia
	  call init_control
      write(6,'(/,a)') 'Output from subroutine ReadInnes:'
!   First we read the basic input; these variables MUST be present.

!     read the title
      READ(5,'(a80)',err=915) title
!     WHICH ATOM (number as in case.struct) :
      READ(5,*,err=915) NATOM
!     WHICH EDGE :
!     nc : main quantum number; lc : orbital quantum number (eg . 1, 0 for K edge)
      READ(5,*,err=915) NC,LC
!     THRESHOLD - DeltaE is the energy loss of the first edge in eV.
      READ(5,*,err=915) DeltaE
!     BEAM ENERGY E0 of the incident electrons in keV :
      READ(5,*,err=915) Energy
!     ENERGY INTERVAL TO BE CALCULATED, all in eV :
      read(5,*,err=916) emin, emax,de
!     semiangles of the spectrometer aperture, and collection angle (mrad)
      READ(5,*,err=915) SpecAperture, ConvergAngle
!     CALCULATION PARAMETER - GRID IN DETECTOR PLANE :
      READ(5,*,err=915) NR, NT
!     SPECTROMETER BROADENING - GAUSSIAN FWHM
      READ(5,*,err=915) specbroad

!   Now process 'advanced' options; these variables are OPTIONAL.

      READ(5,112,end=200,err=916) switch

      do while(switch.ne.'END ')


        if (switch.eq.'DETE') then ! key DETECTOR POSITION
!          specify detector position; angles in mrad.
!          POSITION OF THE DETECTOR :
!          ThetaX (mrad) is the angle between 000 and the EELS aperture
!          in the direction of the second exited spot. ThetaY is the 
!          angle in the orthogonal direction. 
           read(5,*,err=916) ThetaX, ThetaY

        elseif (switch.eq.'ORIE') then ! key ORIENTATION SENSITIVE
!          SAMPLE TO BEAM ORIENTATION
!          read the 3 Euler angles: 
           read(5,*,err=916) Aalpha,ABeta,AGamma
           AAlpha = ConversionDegRad(AAlpha)
           ABeta = ConversionDegRad(ABeta)
           AGamma = ConversionDegRad(AGamma)
		   averaged = .false.

        elseif (switch.eq.'QGRI') then ! key QGRID
		   read(5,*) qmodus
		   if(qmodus.eq.'l') qmodus='L'
		   if(qmodus.eq.'L'.or.qmodus.eq.'1') then
		     read(5,*) Th0
			 Th0=Th0/dble(1000)   ! conversion from mrad to rad
		   elseif(qmodus.ne.'U'.and.qmodus.ne.'u') then
		     write(6,*) ':INFO - unknown qmodus given for input key QGRID.'
			 write(6,*) 'This option will be ignored and the calculation will', &
			   ' continue with a uniform grid of Q-vectors.'
		   endif

        elseif (switch.eq.'SPIN') then ! key SPIN
           nspin=2
		   stop 'Implementation of key SPIN is not finished yet.  Please try again later.'

        elseif (switch.eq.'INIT') then ! key INITIALIZATION
           read(5,*,err=916) h1,h2 !mdos,wdos
		   read(5,*,err=916) h3,h4 !mrot,wrot
		   if (h1.eq.'n'.or.h1.eq.'N') mdos=.false.  ! defaults are .true.
   		   if (h2.eq.'n'.or.h2.eq.'N') wdos=.false.
   		   if (h3.eq.'n'.or.h3.eq.'N') mrot=.false.
   		   if (h4.eq.'n'.or.h4.eq.'N') wrot=.false.
		   if(wdos.and.(.not.mdos)) then
		     write(6,*) ':INFO - you want to write (cross)dos to file without calculating it first.'
			 write(6,*) 'This makes no sense.  The option will be ignored.'
			 wdos=.false.
		   endif
		   if(wrot.and.(.not.mrot)) then
		     write(6,*) ':INFO - you want to write rotation matrices to file without calculating them first.'
			 write(6,*) 'This makes no sense.  The option will be ignored.'
			 wrot=.false.
		   endif

        elseif (switch.eq.'MODU') then ! key MODUS
		   read(5,*,err=916) modus
		   if (modus.eq.'e') modus='E'
		   if (modus.eq.'a') modus='A'
		   if (modus.ne.'E'.and.modus.ne.'A') stop 'alien modus in case.innes'

        elseif (switch.eq.'FERM') then ! key FERMI ENERGY
		   read(5,*) ef
		   EfIsRead=.true.

        elseif (switch.eq.'SELE') then ! key SELECTION RULE
		   read(5,*,err=916) selrule
		   okay=.false.
		   do i=1,20
		     if(selrule.eq.testselrule(i:i)) okay=.true.
		   enddo
		   if(.not.okay) then
		     write(6,*) ':INFO : Invalid option ',selrule,' for SELE in case.innes.'
			 write(6,*) 'This option is ignored and the calculation continues with default settings.'
           endif

        elseif (switch.eq.'LSEL') then ! key LSELECTION RULE
		   read(5,*,err=916) selrule2
		   okay=.false.
		   do i=1,20
		     if(selrule2.eq.testselrule(i:i)) okay=.true.
		   enddo
		   if(.not.okay) then
		     write(6,*) 'Invalid option ',selrule2,' for LSEL in case.innes.'
			 write(6,*) 'This option is ignored and the calculation continues with default settings.'
		   endif

        elseif(switch.eq.'ATOM') then ! key ATOMS
		  read(5,*,err=916) NEqAtDeb, NEqAt

		elseif(switch.eq.'SPLI') then ! key SPLIT
		  if(LC.eq.0) then   ! There is only one edge
		    write(6,*) ':INFO - option SPLIT found in input for edge ',NC,LC,' .'
		    write(6,*) 'Option SPLIT ignored and calculation continues.'
		    read(5,*)   ! Skip one lines
		  else   ! There are two edges
!           ENERGY SPLITTING BETWEEN eg. L3 and L2 :
!           SPLIT is the energy difference between the 2 edges of the same l'
!           there are 2 possible j, first edge: j=l'+1/2, second j=l'-1/2.
!           SPLIT is in eV.
		    read(5,*,err=916) split
		    if(split.lt.dble(0)) then
			  calcsplit = .true.
			else
			  calcsplit = .false.
			endif
          endif

        elseif(switch.eq.'BRAN') then ! key BRANCHING RATIO
		    if(LC.eq.0) then   ! There is only one edge
		    write(6,*) ':INFO - option BRANCHING RATIO found in input for edge ',NC,LC,' .'
		    write(6,*) 'Option BRANCHING RATIO ignored and calculation continues.'
		    read(5,*)   ! Skip one lines
		  else   ! There are two edges
            read(5,*,err=916) BranchingRatio
			if (BranchingRatio.lt.dble(0)) write(6,*) 'Negative Branching Ratio ',  &
			  ' - will be ignored.  Calculation continues.'
		 endif

        elseif(switch.eq.'CORE') then ! key CORE WAVEFUNCTION
		  read(5,*,err=916) corewffile
		  file2core = .true.

		elseif(switch.eq.'FINA') then ! key FINAL STATE WAVEFUNCTION
		  read(5,*,err=916) finalwffile
		  file2final = .true.

        elseif(switch.eq.'NONR') then ! key NONRELATIVISTIC
		  RelatQ=.false.

        elseif(switch.eq.'NOHE') then ! key NOHEADERS
		  headers=.false.

		elseif(switch.eq.'XQTL') then ! key XQTL FILE
!		  read(5,*,err=916) xqtlfile
		  xqtlalaNovak = .true.

	    elseif(switch.eq.'OUTP') then ! key OUTPUT
		  read(5,*,err=916) verbosity
		  if (verbosity.gt.2.or.verbosity.lt.0) then
		    write(6,*,err=916) 'Ignoring unknown verbosity value ',verbosity,' in case.innes.'
		  endif

		elseif(switch.eq.'    ') then ! empty key
!         Do nothing - blank lines are allowed.

                elseif(switch(1:1).eq.'!'.or.switch(1:1).eq.'#') then ! empty key
!         Do nothing - this is a comment line
		

        elseif(switch.eq.'END ') then ! key END
		  exit           ! step out of the loop

        else   ! key unknown
		  write(6,*) 'Unknown key ',switch,' in case.innes.'
          STOP 'Switch not recognized in LectureInnes.'

        endif       ! end of the big list of switches

        read(5,112,end=200,err=916) switch  ! next one ...
      end do

112   format(A4)
200   close(5)          ! Contrary to WIEN habits, we protect our files by closing them after use ...

      write(6,'(a)') 'ReadInnes has finished reading case.innes and will now finish data initialization.'


!     open another qtl-file if necessary
        goto 1212
        if(xqtlalaNovak) then
		  qtlfile(1:80)=fileroot
		  qtlfile(81:90)='          '
		  close(30)
		  do i=80,1,-1
		    if(qtlfile(i:i).eq.'.') exit
		  enddo
		  qtlfile(i+1:i+3)='qtl'
		  if(natom.lt.10) then
		    qtlfile(i+4:i+4)=char(48+natom)
		  elseif(natom.lt.100) then
		    qtlfile(i+4:i+4)=char(48+natom/10)
			qtlfile(i+5:i+5)=char(48+natom-10*(natom/10))
		  elseif(natom.lt.1000) then
		    qtlfile(i+4:i+4)=char(48+natom/100)
			qtlfile(i+5:i+5)=char(48+(natom-100*(natom/100))/10)
			qtlfile(i+6:i+6)=char(48+natom-10*(natom/10))
		  else
		    stop 'Adapt routine readinnes to more than 999 atoms.'
		  endif
		  close(30)
		  open(30,file=qtlfile,form='formatted',status='old')
		endif
1212   continue

      if(.not.EfIsRead) then
!       read Fermi-Energy in Ry from case.qtl .
        READ(30,4714) EF
 4714   FORMAT(//,56x,f10.5)
        REWIND(30)
      endif



!   Some post processing
      call set_nofcross(nspin)       ! number of cross-DOS terms
      call make_wi(BranchingRatio)   ! weights of the edges
      call init_UseThisLambda(lmax,lc,selrule,okay)  !  expansion of exp(iQ.r)
	  if(.not.okay) then
	    write(6,*) 'Error occurred while processing the selection rule ',selrule,' .'
		stop 'selection rule crash'
	  endif
	  call init_UseThisL(lmax,lc,selrule2,okay)      !  expansion (l,m) of final state
	  if(.not.okay) then
	    write(6,*) 'Error occurred while processing the l-selection rule ',selrule2,' .'
		stop 'l-selection rule crash'
	  endif
!     Conversions:
!     E0 in keV -> eV
!     Theta(X,Y) in mrad -> rad.
!     Spectrometer aperture and collection angle in mrad -> rad.
      Energy = Energy * dble(1000)
      ThetaX = ThetaX * dble(0.001)
      ThetaY = ThetaY * dble(0.001)
      SpecAperture = SpecAperture *dble(0.001)
      ConvergAngle = ConvergAngle * dble(0.001)
      K0Len = dble(2)*PI/WaveLength(Energy)
!     Gamma = 1 + E0/(M0 c^2); a0 sets the units (see module units).
      GammaA0 = (dble(1) + Energy/dble(511004)) * A0

      if(file2core.and.calcsplit.and.(LC.gt.0)) then
	    write(6,*) ':INFO : core wave function is taken from file and split energy is not given in input file.'
		write(6,*) 'The calculation will proceed with split energy 0 eV.'
	  endif

      IF ((SpecAperture.GT.1.0D-6).OR.(ConvergAngle.GT.1.0D-6)) THEN
         ThPart = (ConvergAngle + SpecAperture) / DBLE(2*NR)
      ELSEIF((NR+NT).gt.2) then
	     write(6,*) ':INFO - You specified (almost?) zero collection and convergence angle.'
		 write(6,*) 'Your NR and NT values were reset to 1.'
         ThPart = dble(0)
         NR = 1
         NT = 1
      ENDIF


      if(qmodus.eq.'U'.or.qmodus.eq.'L') then
	    npos=nr*nr*nt
	  elseif(qmodus.eq.'1') then
	    nt=1   ! it is not used anyway
		npos=nr
	  endif

!  Write the values of all input parameters to case.outputnes
	  call DescribeTask
      RETURN

915   write(6,*) "You have made a mistake in the basic input block of case.innes."
      write(6,*) "Please have a look at the template files and the Users Guide.  Exiting now."
	  stop
916   write(6,*) 'An error occurred while reading the advanced input options of case.innes.'
      write(6,*) 'Please check out the documentation or the Mailing List for help.  Exiting now.'
	  stop
      END


