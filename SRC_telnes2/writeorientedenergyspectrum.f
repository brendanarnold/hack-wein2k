!BOP
! !ROUTINE: WriteOrientedEnergySpectrum
! !INTERFACE:
    subroutine WriteOrientedEnergySpectrum
! !USES:
    use spectra_hyperfine
    use dimension_constants
    use input
    use energygrid
    use leftovers
    use cross_dos,only : doslm
    use program_control,only : verbosity,fileroot,headers
    use struct,only : zz
! !DESCRIPTION:
!   The energy differential cross section (in the array X) is written to file 47.
!   Depending on verbosity, the following files may be written :
!   Its partial l-components (array SDLM) are written to file 49.
!   Some of its cross term components (array CTR) are written to file 48.
!   The matrix element (i.e., spectrum divided by density of states - array matrix) is written to file 50.
!
!   Format of the files may depend on the edge : for s-edges, there is only one spectrum,
!   while for other edges, there are two spectra (j=l+1/2 and j=l-1/2).
!
!   Additionally, a file case.inb for the broadening program 'broadening' is
!   prepared.  This file contains default parameters.
!
!   This routine differs from WriteAveragedEnergySpectra in only two points :
!   ** It writes complex variables X and CTR
!   ** It has an additional file case.ctr .

! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!   cosmetics March 2005 (Kevin Jorissen)
!EOP

	  implicit none

!   LOCAL VARIABLES
      integer L1,J,lauf,g
      real*8 matrix(0:lmax,jemax)
      real*8,external :: coreholelifetime 


      matrix(:,:)=dble(0)

!   Write case.inb for the broadening program :
      write(8,'(a80)') title
      write(8,'(a)') 'ELNES'
!!   The next instructions write 'case.elnes'.
!      do j=1,80
!        if(fileroot(j:j).eq.'.') exit
!      enddo
!      if(j.lt.76) then
!        fileroot(j+1:j+5)='elnes'
!        if(j.lt.75) fileroot(j+6:80)=' '
!      endif
!      write(8,'(a80)') fileroot
!   The next lines tell the broadening program which columns should be broadened
!   (those containing the l+1/2 and l-1/2 spectra).  This depends on the type of
!   calculation.
      if(lc.gt.0) then
        write(8,'(i3,i3,i3)') 5,3,5
        write(8,'(3(g18.10,x))') split,dble(1),dble(1)
      else
        write(8,'(i3,i3,i3)') 1,1,0
        write(8,'(3(g18.10,x))') dble(0),dble(1),dble(0)
      endif
!   Next we determine the first energy at which the spectrum differs from zero.
!   This is important for energy dependent final state broadening.
      do j=0,jemax
        if (dabs(dble(x(j,1))).gt.dble(0)) exit
      enddo
      if(j.gt.jemax) j=0 ! this only happens when the complete spectrum equals zero ...

      write(8,*) coreholelifetime(zz(natom),nc,lc,0),coreholelifetime(zz(natom),nc,lc,1)     ! core hole lifetime of two edges
      write(8,*) 1,ene(j)  ! 1 means linearly energy dependent valence broadening; ene(j) is edge offset (energy at which valence broadening should start)
      write(8,*) specbroad ! default spectrometer broadening of Gaussian FWHM 0.5 eV
      write(8,'(g18.10,/,g18.10,/,g18.10)') dble(0),dble(0),dble(0)  ! these parameters will not be used
	  close(8)



!  Write case.elnes      

      if (headers) then
!  A header containing the most important parameters of the calculation :
      write(47,'(a)') '### The orientation sensitive spectrum.  (Written as complex data)'
      write(47,'(a)') '## Parameters of the calculation :'
      write(47,4714) NATOM,NEQAT,NC,LC
      write(47,4715) DeltaE, SPLIT
      write(47,4716) Energy/1000
      write(47,4717) ThetaX*1000, ThetaY*1000
      write(47,4720) selrule,selrule2
      write(47,4724) SpecAperture*1000,ConvergAngle*1000
      write(47,4723) AAlpha*57.29577951308232D0,  &
           ABeta*57.29577951308232D0, AGamma*57.29577951308232D0
      write(47,4725) NR, NT
 4714 format('# Atom ',i3,',',i3,'edge   n=',i2,' l=',i2)
 4715 format('# Edge onset ',f6.1, ' eV;  split ',f6.2,' eV')
 4716 format('# Beam energy ',f12.0,' keV')
 4717 format('# Detector position ThetaX ',f6.2,' ThetaY ',f6.2,'  mrad')
 4720 format('# Selection rules ',a,2x,a)
 4724 format('# Collection semiangle ',f6.2,' Convergence semiangle ',f6.2,'  mrad')
 4725 format('# NR ',i3,' NT ',i3)
 4723 format('# Euler Angles: ',f8.2,f8.2,f8.2)
      endif

!  Now write the actual spectra :
      write(47,*)
      if (lc.gt.0) then
        write(47,*) '##  Energy, Re [spectrum], Im  [spectrum], Re [partial spectrum 1], Im [partial spectrum 1],',&
                    ' Re [partial spectrum 2], Im [partial spectrum 2].'
        do J=0,JEMAX
          write(47,4713) ENE(J),X(J,1),X(J,2),X(J,3)
        enddo
      else
        write(47,*) '##  Energy, Re [spectrum], Im  [spectrum].'
        do J=0,JEMAX
          write(47,4713) ENE(J),X(J,1)
        enddo
      endif


      if (verbosity.ge.1) then
	    if (headers) then
!------------header for crossterm file
        write(48,*) '### The cross terms. (Written as complex data)'
        write(48,4810)
!------------header for partial spectrum file
        write(49,*) '### The partial spectra. (Written as real data)'
        write(49,4811)
		endif
!------------header for matrix elements file
        if (lc.eq.0) then
          if (headers) write(50,*) '## Matrix Elements (ie, spectrum / dos).'
          lauf=1
        else
          if (headers) write(50,*) '## Matrix Elements (ie, spectrum / dos) for first edge.'
          lauf=3
        endif
        if (headers) write(50,4812)

        do j=1,jemax
          do l1=0,lmax
            if(doslm(j,l1).lt.0.000001) then
              matrix(l1,j)=dble(0)
            else
!  next section is clumsy, my apologies ...
              if(l1.lt.3) then
                do g=1+l1**2,(l1+1)**2
                  matrix(l1,j)=matrix(l1,j)+sdlm(g,j,lauf)/doslm(j,l1)
                enddo
              else  ! l1=3
                matrix(l1,j)=sdlm(10,j,lauf)/doslm(j,l1)
              endif
            endif
          enddo
        enddo

        do j=1,jemax
          write(49,4814) ENE(J),(SDLM(lauf,J,1),lauf=1,10)	! partial spectrum
          write(48,4813) ENE(J),(CTR(lauf,J,1),lauf=1,12)	! Crossterm
          write(50,4814) ene(j),(matrix(lauf,j),lauf=0,lmax)    ! matrix elements
        enddo

      endif  ! verbosity > 0

      RETURN


 4713 format(F10.5,6e14.5)
 4812 format('#',2x,'ENE',7x,'l=0  ',9x,'l=1  ',9x,'l=2  ',9x,'l=3   ')
 4813 format(F10.5,24e14.5)
 4814 format(F10.5,10e14.5)
 4810 format('#',2x,'ENE',7x,'00 22',23x,'00 21',23x,'00 20',23x, &
      	'11 1-1',22x,'11 10',23x,'22 2-2',22x,'22 21',23x, &
      	'22 2-1',22x,'21 2-1',22x,'20 22',23x,'20 21',23x,'all other cross term contributions')
 4811 format('#',2x,'ENE',7x,'00   ',9x,'1-1  ',9x,'10   ',9x, &
      	'11    ',8x,'2-2  ',9x,'2-1   ',8x,'20   ',9x, &
      	'21    ',8x,'22    ',9x,'sum{m} 3m')

      END



