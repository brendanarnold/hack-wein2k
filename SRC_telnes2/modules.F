!BOP
! !MODULE: modules
! !DESCRIPTION:
!   Contains most of the variables used by the ELNES program.
!
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP




      module constants
! This module contains only simple constants, whose value 
!     should not change in any circumstance,
! and is not used as a dimension parameter in any way.
      real*8, parameter ::  pi = 3.14159265359D0
!   conversion from eV to Ry :
      real*8, parameter :: ev2Ry = 1.0d0/13.6058d0
!  h/2pi c in units eV a.u.
      real*8, parameter :: hbarc = 1973.2708d0/0.529177d0
! electron rest mass times c^2 in au (ie, 1 * alfa * alfa), times eV/Ha (27.2)
      real*8, parameter :: MeC2 =  511004.0d0                 
      end module constants

!*****************************************************************
      module dimension_constants
! This module absorbs the old param.inc file.  It contains compilation constants that are used
! mostly as dimension of arrays and vectors in the code.  They may have impact on memory use and speed.
! maximal number of radial points considered
      integer, parameter ::        nrad = 881
! maximal orbital quantum number l for the final states, the DOS and the spectrum	          
      integer, parameter ::        lmax = 3
! number of cross terms (reduced because of symmetry to diagonal + upper corner)
      integer                      nofcross
      CONTAINS
	   subroutine set_nofcross(nspin)
	    implicit none
	    integer, intent(in) :: nspin
        
		nofcross = nspin*(lmax+1)**2*(nspin*(lmax+1)**2+1)/2
        return
	   end subroutine set_nofcross

      end module dimension_constants

!****************************************************************
      module units
!  This module specifies the units of the results.
! We work in atomic units :
      real*8, parameter ::         a0 = 1.0d0        
      end module units

!     ***********************************************************
      module work_atpar
!   This module contains work arrays for the atpar routine.
!   They are not required in the rest of ELNES, and may be destroyed when leaving atpar.
      real*8, allocatable ::       a(:),a1(:,:),ae(:)
      real*8, allocatable ::       b(:),b1(:,:),be(:)
      real*8, allocatable ::       da1(:,:),da(:)
      real*8, allocatable ::       db1(:,:),db(:)
      
      CONTAINS
      subroutine make_abv(nrad,iemax0)
!   allocate work arrays for atpar
      implicit none
      integer nrad,iemax0
      allocate (a(nrad),a1(nrad,iemax0),ae(nrad))
      allocate (b(nrad),b1(nrad,iemax0),be(nrad))
      allocate (da1(nrad,iemax0),da(nrad))
      allocate (db1(nrad,iemax0),db(nrad))
      da1(:,:)=dble(0);da(nrad)=dble(0)
      db1(:,:)=dble(0);db(nrad)=dble(0)
      a(:)=dble(0);a1(:,:)=dble(0);ae(:)=dble(0)
      b(:)=dble(0);b1(:,:)=dble(0);be(:)=dble(0)
      end subroutine make_abv
      subroutine destroy_abv
!   deallocate work arrays for atpar
	     deallocate(a,a1,ae,b,b1,be,da,da1,db,db1)
	   end subroutine destroy_abv
      end module work_atpar
       
!     ***********************************************************
      module radial_functions
!  This module contains the APW radial basis functions for all l-values and energies.
!     Large (Uell) and small (UPell) components of the radial function :
	  real*8, allocatable ::   Uell(:,:,:), UPell(:,:,:)
!     Norm of the radial function :
	  real*8, allocatable ::   UellUell(:,:,:)
!     Derivatives of the radial function :
      real*8, allocatable ::   DUell(:,:,:),DUPell(:,:,:)
!     Radial functions divided by the radial coordinate r :
      real*8, allocatable ::   Uellr(:,:,:),UPellr(:,:,:)


	  CONTAINS
	   subroutine make_uel(nrad,iemax,lmax)
!     Allocate arrays for APW radial basis functions :
	     integer nrad,iemax,lmax
	     allocate(Uell(nrad,iemax,0:lmax),UPell(nrad,iemax,0:lmax), &
	         UellUell(iemax,0:lmax,0:lmax))
		 allocate(DUell(nrad,iemax,0:lmax),DUPell(nrad,iemax,0:lmax))
		 allocate(Uellr(nrad,iemax,0:lmax),UPellr(nrad,iemax,0:lmax))
		 uell(:,:,:)=dble(0);upell(:,:,:)=dble(0);uelluell(:,:,:)=dble(0)
		 DUell(:,:,:)=dble(0);DUPell(:,:,:)=dble(0)
		 Uellr(:,:,:)=dble(0);UPellr(:,:,:)=dble(0)
       end subroutine make_uel
	  end module radial_functions

!     ***********************************************************
      module energygrid
!  This module contains the energy grid and the Fermi energy.
!
!     The Fermi energy :      
      real*8 :: ef
!     the energy grid
	  real*8, allocatable :: ene(:)
!     indices specifying the energy grid - index of first and last energy point used; total number of points
!     If DOS is read from file, the 'available' grid may differ from that requested in case.innes.	      
      integer :: jemin, jemax, iemax
!     energy split (eg. L2-L3) in pixels
	  integer :: jdecal                
      CONTAINS
	   subroutine make_ene(iemax)
!     allocate energy grid
        integer iemax
	    allocate(ene(0:iemax))
		ene(:)=dble(0)
	   end subroutine make_ene
	  end module energygrid

!     ***********************************************************
      module potential
!  This module contains the spherical potential.
!     VR(j,i) - spherical part (l=0, m=0) of the total potential r*V at mesh point j for atom i
      real*8, allocatable ::   VR(:,:)
	  CONTAINS
	   subroutine make_pot(nato,nrad)
!     allocate the spherical potential :
	     integer nato,nrad
	     allocate(vr(nrad,nato))
		 vr(:,:)=dble(0)
	   end subroutine make_pot
	   subroutine destroy_pot
	     deallocate(vr)
	   end subroutine destroy_pot
	  end module potential

!     ***********************************************************
      module input
!     This module contains input variables from case.innes.
!     We put the 'numerical values' in this module, while the
!     module program_control groups variables which control the
!     program control.
!
!     Number of the excited atom (from case.struct)
      integer          NATOM
!     Equivalent positions NEqAtDeb to NEqAt will be summed	  
	  integer          NEqAt, NEqAtDeb
!     The initial state has main quantum number nc and orbital quantum number lc
      integer          NC, LC
!     The energy grid goes from emin to emax in steps of de (all in eV)
      real*8           EMIN, DE, EMAX
!     The edge onset (of the j=l+1/2 edge) is at deltae (in eV)
      real*8           DeltaE
!     The edge j=l-1/2 is at energy deltae + split (in eV)
      real*8           Split
!     The incoming beam has energy Energy (in eV after ReadInnes), corresponding to
!     wavevector K0Len (in au^-1)
      real*8           Energy,K0Len
!     The detector is positioned at angles ThetaX,ThetaY from the 000 spot (in rad after ReadInnes)
      real*8           ThetaX, ThetaY
!     The sample to beam orientation is described by three Euler angles (in rad after ReadInnes)
      real*8           AAlpha, ABeta, AGamma
!     The collection semiangle (in rad after ReadInnes)
      real*8           SpecAperture
!     The convergence semiangle (in rad after ReadInnes)
      real*8           ConvergAngle
!     The separation between consecutive sampling points for uniform Q-meshes (in rad after ReadInnes)
      real*8           ThPart
!     The first sampling point for logarithmic Q-meshes (in rad after ReadInnes)
      real*8           Th0
!     The selection rule for the expansion of the exponential factor
      character*1      selrule
!     The selection rule for the final states
      character*1      selrule2
!     The title of the calculation
	  character*80     title
!     Read the core wave function from this file
      character*80     corewffile
!     Read the final state radial functions from this file
      character*80     finalwffile
!     Read the xqtl matrix (qtl components) from this file
      character*80     xqtlfile
!     The Q-mesh contains NPos vectors; the mesh size is NR (radial) and NT (angular)
      integer          NR, NT, NPos
!     The spectrometer broadening (FWHM of Gaussian broadening function)
      real*8           specbroad
!     The branching ratio and the prefactors weighing the j=l+1/2 (w1) and the j=l-1/2 (w2) edge
	  real*8           BranchingRatio,w1,w2
!     Array controlling which l-values for the final states can contribute
	  logical, allocatable :: UseThisL(:)
!     Array controlling which lambda-values can contribute (ie, terms in the expansion of the exponential exp(iQ.r))
      logical, allocatable :: UseThisLambda(:)
!     The number of spin states - 1 for no spins; 2 for +/- 1/2
      integer          nspin

	  CONTAINS
	    subroutine init_varia
!     Initialize all data for the ELNES program
		  use dimension_constants,only : lmax
		  implicit none
		  NR=1; NT=1; NPos=1
          selrule='d'  ! d for dipole selection rule
		  selrule2='n' ! all final l allowed
		  SpecAperture=dble(0);ConvergAngle=dble(0);ThPart=dble(0);Th0=dble(0)
		  AAlpha=dble(0);ABeta=dble(0);AGamma=dble(0)
		  ThetaX=dble(0);ThetaY=dble(0)
		  Energy=dble(0);K0Len=dble(0);SPLIT=dble(0);DeltaE=dble(0)
		  EMIN=dble(0);EMAX=dble(15);DE=dble(0.05)
		  NATOM=0;NEqAt=0;NEqAtDeb=1
		  NC=0;LC=0
		  nspin=1
		  specbroad=dble(0.5)
		  title(:)=' ';corewffile(:)=' ';finalwffile(:)=' ';xqtlfile(:)=' '
		  BranchingRatio=dble(-1)
	    end subroutine init_varia
        subroutine make_wi(b)
!     construct weights w1,w2 that satisfy the requested branching ratio b
		  real*8 b
		  if (b.lt.dble(0)) then
		    w1=dble(2*LC+2)/dble(4*LC+2)
			w2=dble(2*LC)/dble(4*LC+2)
			b=w1/w2
		  else
		    w1=b/(b+dble(1))
			w2=dble(1)/(b+dble(1))
		  endif
		end subroutine make_wi
		subroutine init_UseThisL(lmax,lc,selrule,okay)
!    from the selection rule on final states, decide which final states may contribute
		  implicit none
		  integer lmax,lc,ldmax,i
		  character*1 selrule
		  logical okay
		  allocate(UseThisL(0:lmax))
		  UseThisL(:)=.false.
		  okay=.true.
		  if(selrule.eq.'n'.or.selrule.eq.'N') then       ! allow all transitions
		    UseThisL(:)=.true.
		  elseif(selrule.eq.'m'.or.selrule.eq.'M') then   ! allow only l=lc transition
		    UseThisL(lc)=.true.
		  elseif(selrule.eq.'d'.or.selrule.eq.'D') then   ! allow only |l-lc|=1 transition
		    if(lc.lt.lmax) UseThisL(lc+1)=.true.
			if(lc.gt.0) UseThisL(lc-1)=.true.
		  elseif(selrule.eq.'q'.or.selrule.eq.'Q') then   ! allow only |l-lc|=2 transition
		    if(lc.lt.(lmax-1)) UseThisL(lc+2)=.true.
			if(lc.gt.1) UseThisL(lc-2)=.true.
		  elseif(selrule.eq.'o'.or.selrule.eq.'O') then   ! allow only |l-lc|=3 transition
		    if(lc.lt.(lmax-2)) UseThisL(lc+3)=.true.
			if(lc.gt.2) UseThisL(lc-3)=.true.
		  else                           ! allow transitions for which |l-lc| =< ldmax
            okay=.false.
		  	do i=0,9
		      if (int(ichar(selrule)-48).eq.i) then
		        ldmax=i
		        okay=.true.
		      endif
		    enddo
            if (okay) then
		      do i=lc-ldmax,lc+ldmax
			    if(i.ge.0.and.i.le.lmax) UseThisL(i)=.true.
			  enddo
		    endif
		  endif
		end subroutine init_UseThisL
		subroutine init_UseThisLambda(lmax,lc,selrule,okay)
!       from the selection rule on the expansion of the exponential, decide which terms may contribute
		  implicit none
		  integer lmax,lc,lambdamax,i
		  character*1 selrule
		  logical okay
		  allocate(UseThisLambda(0:(lmax+lc)))            ! because of triangle, no higher lambda values are allowed
		  UseThisLambda(:)=.false.
		  okay=.true.
		  if(selrule.eq.'n'.or.selrule.eq.'N') then       ! allow all transitions
		    UseThisLambda(:)=.true.
		  elseif(selrule.eq.'m'.or.selrule.eq.'M') then   ! allow only monopole transition
		    UseThisLambda(0)=.true.
		  elseif(selrule.eq.'d'.or.selrule.eq.'D') then   ! allow only dipole transition
			UseThisLambda(1)=.true.
		  elseif(selrule.eq.'q'.or.selrule.eq.'Q') then   ! allow only quadrupole transition
			UseThisLambda(2)=.true.
		  elseif(selrule.eq.'o'.or.selrule.eq.'O') then   ! allow only octopole transition
			UseThisLambda(3)=.true.
		  else                            ! allow all transitions up to lambdamax
            okay=.false.
		  	do i=0,9
		      if (int(ichar(selrule)-48).eq.i) then
		        lambdamax=min(i,lmax+lc)
		        okay=.true.
		      endif
		    enddo
            if (okay) UseThisLambda(0:lambdamax)=.true.
		  endif

		end subroutine init_UseThisLambda
	  end module input

! ***************************************************************
	  module program_control
!     This module contains input variables from case.innes.
!     We put the 'numerical values' in module input, while this
!     module groups variables which control the program control.
!     energy resolved (modus=e) or angular resolved (modus=a) spectrum
      character*1      modus
!     what kind of q-mesh : uniform (U), logarithmic (L), or one dimensional logarithmic (1)
      character*1      qmodus
!     for making new filenames :
      character*80     fileroot
!     Calculate rotation matrices?
	  logical          mrot
!     Write rotation matrices?
      logical          wrot
!     Calculate density of states?
      logical          mdos
!     Write density of states?
      logical          wdos
!     Calculate the j=l+/-1/2 energy splitting?
	  logical          calcsplit
!     Read core wave function from file?
      logical          file2core
!     Read final state radial functions from file?
      logical          file2final
!     Read xqtl from P. Novak's qtl program?
      logical          xqtlalaNovak
!     Use relativistic corrections to the cross section formula?
	  logical          RelatQ
!     Controls the amount of output given by ELNES
	  integer          verbosity
!     Average the spectrum?  (false -> orientation sensitive spectrum)
	  logical          averaged
!     Allows to write headers to output files
      logical          headers

      CONTAINS
	    subroutine init_control
!        Initialize the variables for program flow
!        These defaults correspond more or less to 'traditional' use of telnes (except for calcsplit)
		 verbosity=0		            
		 averaged=.true.; modus='E'; qmodus='U'  ! N for no orientation; E for energy; U for uniform grid 
		 mrot=.true.;mdos=.true.;wdos=.true.;wrot=.true.
		 calcsplit=.true.;file2core=.false.;file2final=.false.;xqtlalaNovak=.false.
		 RelatQ=.true.
		 headers=.true.
        end subroutine init_control
	  end module program_control

!     ***********************************************************
      module qvectors
!     This module contains the impuls transfer vectors and the final wave vectors.

!     Weights for integration over Q-vectors
	  real*8, allocatable ::  WeightV(:)
!     Carthesian coordinates of final wave vectors in the 'detector plane' (careful : this is in a 
!     plane made by convoluting the beam cross section and the detector aperture)
      real*8, allocatable ::  ThXV(:), ThYV(:)
!     The Q-mesh for a particular energy
	  real*8, allocatable ::  QV(:,:,:)
!     The length of each Q-vector, both relativistically and 'classically'
      real*8, allocatable ::  QLenV(:,:),QLenVClas(:,:)

	  CONTAINS
	   subroutine make_Qvecs1(nposmax)
!     allocate the final wave vectors
	     integer nposmax
	     allocate(WeightV(NPOSMAX), ThXV(NPOSMAX), ThYV(NPOSMAX))
		 weightv(:)=dble(0);thxv(:)=dble(0);thyv(:)=dble(0)
	   end subroutine make_Qvecs1
	   subroutine make_Qvecs2(nposmax,ndif)
!     allocate the Q-mesh
	     integer nposmax,ndif
	     allocate(QLenV(ndif,nposmax),qv(3,ndif,nposmax),QLenVClas(ndif,nposmax))
		 qlenv(:,:)=dble(0);qv(:,:,:)=dble(0)
		 qlenvclas=dble(0)
	   end subroutine make_Qvecs2
       subroutine destroy_qvecs
!     deallocate final wave vectors and Q-mesh
	     deallocate(weightv,thxv,thyv,qlenv,qv,qlenvclas)
	   end subroutine destroy_qvecs
	  end module qvectors

!     ***********************************************************
      module initialstate1
!   This module contains the j=l+1/2 core state.
!   It also contains some work arrays for the hfsd routine.
!     The spherical potential and the radial mesh
      real*8, allocatable :: dv(:),dr(:)
!     Large and small component of the core state wave function
      real*8, allocatable :: dp(:),dq(:)
!     Step size of the exponential grid dr; atomic number; some precision parameter
	  real*8 dpas,z,tets
!     Several integers ;-)
	  integer nstop,nes,np,nuc
	  CONTAINS
	   subroutine make_dp(nrad)
!     Allocate variables for the j=l+1/2 core state
	    integer nrad
	    allocate(DV(NRAD+41),DR(NRAD+41),DP(NRAD+41),DQ(NRAD+41))
		dv(:)=dble(0);dr(:)=dble(0);dp(:)=dble(0);dq(:)=dble(0)
	   end subroutine make_dp
	  end module initialstate1

!*****************************************************************
      module initialstate2
!     This module contains the j=l-1/2 core state.
!     Large and small component of the core state wave function
      real*8, allocatable ::  dpdp(:), dqdq(:)
	  CONTAINS
	   subroutine make_dpdp(nrad)
	     integer nrad
	     allocate(dpdp(nrad+41),dqdq(nrad+41))
!     Allocate variables for the j=l-1/2 core state
	   end subroutine make_dpdp
	  end module initialstate2
      
!******************************************************************
	  module leftovers
!   This module contains a few variables that do not belong in other modules.
!     The prefactor gamma times Bohr radius for the spectrum
      real*8 GammaA0
!     Derivatives of the Feinstruktur constant, used for integration in rint13 and rint14.
      real*8 CFEIN1, CFEIN2
	  end module leftovers

!     ***********************************************************
      module densityofstates
!   This module contains the *regular* density of states (no cross terms).
!     DOS() Density of States, l=0 -> l=lmax
!     The DOS-array contains the l-partial DOS :  l=0 l=1 l=2 l=3

	  real*8, allocatable  :: DOS(:,:)

	  CONTAINS
	   subroutine make_dos(iemax,lmax)
!     allocate the DOS array
	    integer iemax,lmax
	    allocate(dos(iemax,0:lmax))
        dos(:,:)=dble(0)
	   end subroutine make_dos

	  end module densityofstates

!******************************************************************
	  module spectra_normal
!  This module contains the variables for averaged spectra.
!     spectrum(energy,one of 1 or 2 edges or the sum)
      real*8, allocatable :: x(:,:)
!     partial spectrum (energy, l-component, which edge)
	  real*8, allocatable :: sdlm(:,:,:) 
	  CONTAINS
	   subroutine makespec_normal(iemax,lmax,ispin)
!     Initialize spectrum for averaged calculation.
         integer iemax,lmax,ispin
		 allocate(x(0:iemax,ispin))
		 allocate(sdlm(0:lmax,0:iemax,ispin))
          x(:,:)=dble(0)
		  sdlm(:,:,:)=dble(0)
	   end subroutine makespec_normal
	  end module spectra_normal
! **************************************************************************
	  module spectra_hyperfine
!  This module contains the variables for orientation sensitive spectra.
!     spectrum(energy,one of 1 or 2 edges or the sum)
	  complex*16, allocatable :: x(:,:)
!     cross term contributions to the spectrum
	  complex*16, allocatable :: ctr(:,:,:)
!     partial spectrum (energy,l-component,which edge) - I expect these to be real; we may make them complex later.
	  real*8, allocatable :: sdlm(:,:,:)  
	  CONTAINS
	   subroutine makespec_hyperfine(iemax,lmax,ispin)
!     Prepare spectrum for orientation sensitive spectra
         integer iemax,lmax,ispin
		 allocate(x(0:iemax,ispin))
		 allocate(sdlm(1:10,0:iemax,ispin))   ! some function of lmax would be nicer than the 1:10, but anyway.
         allocate(ctr(1:12,0:iemax,ispin))
          x(:,:)=dcmplx(0,0)
		  sdlm(:,:,:)=dble(0)
		  ctr(:,:,:)=dcmplx(0,0)
	   end subroutine makespec_hyperfine
	  end module spectra_hyperfine

!     ***********************************************************
      module struct
!  This module contains information about the crystal structure.
!  It is initialized by the Inilpw routine.
!     number of equivalency classes (nat) and total number of atom positions in the unit cell (nats)
      integer nat,nats         ! number of equivalency classes and number of atom positions in the unit cell
!     number of symmetry operations in the space group
      integer nsym
!     true for relativistic calculations (in the sense that wave functions have large and small component)
      logical rel              
!     atomic number, for every equivalency class
	  real*8, allocatable :: zz(:)
!     useless index (iatnr), multiplicity of each equivalency class, charge splitting option (isplit)	     
      integer, allocatable  ::     IATNR(:), MULT(:),ISPLIT(:)
!     lattice parameters alat and lattice angles alpha
      real*8  alat(3), alpha(3)
!     muffin tin radius for every equivalence class
	  real*8, allocatable :: rmt(:)
!     atom position for every atom
	  real*8, allocatable :: pos(:,:)
!     lattice type
      character*4 lattic
!     label of each atom	     
      character*10,allocatable ::    NAME(:)
!     symmetry matrices
	  integer, allocatable :: iz(:,:,:)
!     symmetry translations	    
	  real*8, allocatable :: tau(:,:)
!     number of radial (logarithmic) mesh points	      
      integer, allocatable ::  JRJ(:)
!     step size of logarithmic radial mesh, first radial mesh point
      real*8, allocatable ::   DH(:), RO(:)

	  CONTAINS
	   subroutine make_struct(nato,ndif)
!     initialize arrays for crystal structure
        integer nato,ndif
	    allocate(rmt(nato),pos(3,ndif))
	    allocate(iatnr(nato),mult(nato),isplit(nato))
        allocate(zz(nato))
	    allocate(jrj(nato),dh(nato),rO(nato))
		allocate(name(nato))
		name(:)=' '
		rmt(:)=dble(0);pos(:,:)=dble(0)
		iatnr(:)=0;mult(:)=0;isplit(:)=0;zz(:)=dble(0)   ! lattic may have been read already, so we take care not to overwrite it.
		jrj(:)=0;dh(:)=dble(0);rO(:)=dble(0)
	   end subroutine make_struct
	   subroutine make_symmat(nsym)
!     initialize symmetry matrices
	    integer nsym
		allocate(iz(3,3,nsym))
		allocate(tau(3,nsym))
		iz=0;tau=dble(0)
	  end subroutine make_symmat
	  end module struct

!     ***********************************************************
      module cross_DOS
!  This module contains the cross DOS and the l,m partial DOS.
!     Cross DOS array
      complex*16, allocatable :: xdos(:,:)
!     l,m resolved DOS
	  real*8, allocatable :: doslm(:,:)
!     an index used to access the xdos matrix
      integer, allocatable :: indexlm(:,:)
	  CONTAINS
	   subroutine make_crossdos(iemax,nofcross,lmax)
!     allocate cross dos arrays
	    integer iemax,nofcross,lmax
		allocate(xdos(iemax,nofcross))
		allocate(doslm(iemax,0:12))  ! more generally, the 12 would be replaced by a function of lmax ...
		doslm(:,:)=dble(0)
		allocate(indexlm(0:lmax,-lmax:lmax))
		indexlm(0:lmax,-lmax:lmax)=0
        xdos(:,:)=dcmplx(0,0)
	   end subroutine make_crossdos
	  end module cross_DOS

!     ***********************************************************
      module rotation_matrices
!   Surprisingly, this routine contains rotation matrices !
!   They are used to transform Q-vectors from one atom to another.

!     Rotation matrices from atom to first atom in its equivalence class
      real*8, allocatable :: rotij(:,:,:)
!     Corresponding translations
	  real*8, allocatable :: tauij(:,:)
!     from first atom in the equivalency class to the global crystal coordinate system
	  real*8, allocatable :: rotloc(:,:,:)
!     from laboratory frame to frame of a particular atom	   
	  real*8, allocatable :: GeneralM(:,:,:)
!     Bravais matrix	   
      real*8                 br1(3,3)   
	  CONTAINS
	   subroutine make_rot(ndif,nato)
!     allocate rotation matrices
	    integer nato,ndif
	    allocate(rotij(3,3,ndif),tauij(3,ndif),rotloc(3,3,nato),generalm(3,3,ndif))
		rotij(:,:,:)=dble(0);rotloc(:,:,:)=dble(0)
		generalm(:,:,:)=dble(0)
	   end subroutine make_rot

	  end module rotation_matrices

!     ***********************************************************
      module crash
!  This modules enables a program unit to crash making calls like
!         CALL ERRFLG(ERRFN,'Error in this program unit')
	  character*80 errfn    ! name of the file for error messages
	  character*200 errmsg  ! error message
	  end module crash

!****************************************************************
      module tetra_params
!     This module contains array dimensions for subroutine tetraforelnes.
!     It replaces the SRC_tetra/param.inc file.
!     It also contains the idos-array which is otherwise read from case.int.
      use dimension_constants,only : lmax
!     maximum number of DOS cases to calculate :
      integer        mg
!     number of cross-qtl terms and its square
      integer        lxdos,lxdos2
!     array specifying which dos to calculate
	  integer, allocatable :: idos(:,:)
!     this array is not really used
	  character*6, allocatable :: dostyp(:)
	  contains
	    subroutine init_tetra_params(m,isplit)
!     initialize some stuff for tetraforelnes.
          use input,only : nspin
		  implicit none
		  integer m,isplit
		  if(isplit.eq.99.or.isplit.eq.88) then
		    lxdos=3
		  else
		    lxdos=1
		  endif
		  allocate(idos(m,2))
		  allocate(dostyp(m))
		  dostyp(:)='      '
		  idos(:,:)=0
		  lxdos2=nspin*(lxdos+1)**2  ! added nspin here ; lxdos2 = 4,8 / 16,32
		  mg=m
		end subroutine init_tetra_params
	  end module tetra_params


! **********************************************************************************
! Next module for tetra, taken from SRC_tetra/reallocate.f of wien2k_04.10 :
! But added extra routine because some r4 was changed to r8.

	module reallocate
	  !     nur 1 (generischer) Name wird von aussen angesprochen
	  interface doreallocate
	    module procedure doreallocate_r8_d1
	    module procedure doreallocate_r8_d2
		module procedure doreallocate_r8_d3  ! added KJ
	    module procedure doreallocate_r4_d2
	    module procedure doreallocate_r4_d3
	    module procedure doreallocate_i4_d1
	    module procedure doreallocate_i4_d2
	  end interface
	contains

	  !     leider sind mehrere subroutines notwendig fuer verschiedene Typen
	  subroutine doreallocate_r8_d1(tf, newdimension)
	    real*8, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
	    !     nur 1 mal kopieren reicht
	    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
	    real*8, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

! next subroutine added KJ
      subroutine doreallocate_r8_d3(tf,newdimension1,newdimension2,newdimension3)
	    implicit none
	    real*8, pointer :: hilsfeld(:,:,:), tf(:,:,:)
		integer newdimension1,newdimension2,newdimension3,min1,min2,min3
		allocate(hilsfeld(newdimension1,newdimension2,newdimension3))
		min1=min(newdimension1,size(tf,1))
		min2=min(newdimension2,size(tf,2))
		min3=min(newdimension3,size(tf,3))
		hilsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
		deallocate(tf)
		tf=>hilsfeld
	  end subroutine 
! end addition KJ


	  subroutine doreallocate_r4_d2(tf, newdimension1, newdimension2)
	    real*4, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_r4_d3(tf, newdimension1, newdimension2, newdimension3)
	    real*4, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2,newdimension3))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            min3=min(newdimension3,size(tf,3))
	    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d1(tf, newdimension)
	    integer*4, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
	    integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	end module 
