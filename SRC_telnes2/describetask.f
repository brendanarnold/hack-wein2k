!BOP
! !ROUTINE: DescribeTask
! !INTERFACE:
      SUBROUTINE DescribeTask
! !USES:
	  use input
	  use program_control
	  use dimension_constants,only : lmax,nofcross
	  use constants
	  use units
! !DESCRIPTION:
!   Writes a precise summary to the main output file 6
!   of the tasks that will be performed, listing all input
!   variables, including those for which default settings
!   are used.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none
!  LOCAL VARIABLES :
      integer l,lamin,lamax


      write(6,'(//,a)') 'Output of subroutine DescribeTask.'
	  write(6,'(a)') '********** Description of the tasks as read from case.innes. ********'
	  write(6,'(a,a,a)') 'The name of your calculation is : ',title,' .'
	  write(6,'(a,i4,a)') 'You are calculating a spectrum of atom ',natom,' .'
	  write(6,'(a,i4,a,i4,a)') 'The contributions of its equivalent positions ',NEqAtDeb,' to ', &
	                NEqAt,' will be summed.'
	  write(6,'(a,a,i2,a,i2,a,f9.3,a)') 'The edge corresponds to excitations from an initial state with quantum numbers ', &
	                'n= ',nc,' , l= ',lc,' and has its threshold at ',DeltaE,' eV.'
	  write(6,'(a,f10.1,a)') 'You are using an electron beam of energy ',Energy,' eV'
	  write(6,'(a,f10.5,a)') '(i.e., wave vector before scattering = ',K0Len,' / atomic units).'
	  write(6,'(a,f8.3,a,f8.3,a,a,f8.3,a,i7,a,/)') 'The edge will be calculated from ',Emin,' to ',Emax,' eV above the Fermi energy', &
	             ' in steps of ',de,' eV  (i.e., ',1+int((emax-emin)/de),' energy points).'

      if(LC.gt.0) then
	    write(6,'(a)') 'Two edges will be calculated (j=l+/-1/2).'
	    if(calcsplit) then
		  write(6,'(a)') 'The energy splitting will be calculated by lcore.'
		else
		  write(6,'(a,f8.3,a)') 'The energy splitting is set to ',split,' eV.'
		endif
        write(6,'(a,f8.3,a)') 'The branching ratio is ',BranchingRatio,' .'
	  endif

	  write(6,'(/,a,f8.3,a)') 'The detector collection semiangle is ',SpecAperture*dble(1000),' mrad.'
      write(6,'(a,f8.3,a)') 'The microscope convergence semiangle is ',ConvergAngle*dble(1000),' mrad.'
	  write(6,'(a,f8.3,a,f8.3,a)') 'The detector is positioned at angles ',ThetaX*dble(1000),' , ', &
                 ThetaY*dble(1000),' in mrad w.r.t. the 000 spot.'
	  if (.not.averaged) then
	    write(6,'(a,a,f8.4,a,f8.4,a,f8.4,a)') 'The calculation is orientation resolved.  The sample to beam orientation is ', &
		 'described by the Euler angles ',AAlpha,' , ',ABeta,' , ',Agamma,' (in rad).'
	  elseif(averaged) then
	    write(6,'(a)') 'The calculation averages over all sample to beam orientations.'
	  endif
      write(6,'(a,a)') 'For sampling of the impuls transfer vectors allowed by these settings, ',&
	    'you have chosen a :'
	  if(qmodus.eq.'U') then
	     write(6,'(a,i6,a,i6,a,i3,a)') 'uniform mesh of ',nr*nr*nt,' points with nr = ',nr,' , nt = ',nt,' .'
      elseif(qmodus.eq.'1') then
	     write(6,'(a,i6,a,e10.5,a)') 'one dimensional logarithmic mesh of ',nr,' points, starting with theta = '&
                 ,Th0*dble(1000),' mrad.'
      elseif(qmodus.eq.'L') then
	     write(6,'(a,i6,a,i6,a,i3,a,e10.5,a)') 'logarithmic mesh of ',nr*nr*nt,' points with nr = ',nr,&
                 ' , nt = ',nt,' ,starting with theta = ', &
		   Th0*dble(1000),' mrad.'
      endif

	  write(6,'(/,a)') 'Allowed contributions to the spectrum :'
	  lamin=100
	  lamax=0
      do l=0,lc+lmax
	   if(UseThisLambda(l).and.iabs(l-lc).lt.lamin) lamin=iabs(l-lc)
	   if(UseThisLambda(l).and.(l+lc).gt.lamax) lamax=l+lc
       if(UseThisLambda(l).and.l.eq.0) write(6,'(a)') 'monopole term'
       if(UseThisLambda(l).and.l.eq.1) write(6,'(a)') 'dipole terms'
       if(UseThisLambda(l).and.l.eq.2) write(6,'(a)') 'quadrupole terms'
       if(UseThisLambda(l).and.l.eq.3) write(6,'(a)') 'octopole terms'
       if(UseThisLambda(l).and.l.gt.3) write(6,'(a,i3)') 'terms of order ',l
      enddo
      write(6,'(a)') 'Final states of this l-character will contribute :'
	  do l=lamin,min(lmax,lamax)
	    if(UseThisL(l)) write(6,'(x,a,i4)') 'l= ',l
	  enddo
      if(nspin.eq.2) then
	    write(6,'(a)') 'Summation over spin up and down states is done explicitly.'
	  elseif(nspin.eq.1) then
	    write(6,'(a)') 'No spin variables are taken into account.'
	  endif
      if(.not.averaged) write(6,'(a,i5,a)') 'Using ',nofcross,' cross-DOS components.'

      if(.not.RelatQ) then
	    write(6,'(a)') ':WARNING : You have switched off relativistic corrections to the scattering cross section.'
		write(6,'(a)') 'The angular behaviour of the calculated cross sections may be wrong for certain cases.'
		write(6,'(a)') 'Please refrain from using the NONRELATIVISTIC key unless you know what you are doing.'
	  endif

      write(6,*)
	  if(verbosity.eq.0) then
	    write(6,'(a)') 'You will get only basic output.'
	  elseif(verbosity.eq.1) then
	    write(6,'(a)') 'You will get basic and more advanced output.'
	  elseif(verbosity.eq.2) then
	    write(6,'(a)') 'You will get full output, including possibly less useful things.'
	  endif
	  if(modus.eq.'E') then
	    write(6,'(a)') 'The differential cross section will be given as a function of energy loss.'
	  elseif(modus.eq.'A') then
	    write(6,'(a)') 'The differential cross section will be given as a function of scattering angle.'
	  endif
	  if(file2core) then
	    write(6,'(a,a)') 'The initial state wavefunction will be read from the file ',corewffile,' .'
	  else
	    write(6,'(a)') 'The initial state wavefunction will be calculated by lcore.'
	  endif
	  if(file2final) then
	    write(6,'(a,a)') 'Final state radial functions will be read from the file ',file2final,' .'
	  else
	    write(6,'(a)') 'Final state radial functions will be calculated by atpar.'
	  endif
	  if(mrot) then 
	    write(6,'(a)') 'The program will calculate rotation matrices itself.'
	  else
	    write(6,'(a)') 'The rotation matrices will be taken from case.rotij .'  
	  endif
	  if(mdos) then
	    write(6,'(a)') 'Density of States or Cross Density of States will be calculated.'
	    if(.not.xqtlalaNovak) then
	      write(6,'(a)') 'The program will use lapw2 qtl data.'
	    else
	      write(6,'(a)') 'The program will use xqtl components as given by P. Novaks qtl program.'
	    endif
	  else
	    write(6,'(a)') 'Density of States or Cross Density of States will be taken from file.'
	  endif

	  write(6,'(/,a,a,f8.3)') 'The units in which we are working, are obtained from atomic units ',&
	     'by multiplication with ',a0,' .'


	  write(6,'(/,a)') '******* End of Task Description ******* '

	  return
	  end
