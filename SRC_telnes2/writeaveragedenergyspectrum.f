!BOP
! !ROUTINE: WriteAveragedEnergySpectrum
! !INTERFACE:
    subroutine WriteAveragedEnergySpectrum
! !USES:
    use dimension_constants, only : lmax
    use energygrid
    use spectra_normal
    use input
    use program_control,only : verbosity,fileroot,headers
    use leftovers
    use densityofstates,only : dos
    use struct,only : zz
! !DESCRIPTION:
!   The energy differential cross section (in the array X) is written to file 47.
!   Depending on verbosity, the following files may be written :
!   Its partial l-components (array SDLM) are written to file 49.
!   The matrix element (i.e., spectrum divided by density of states - array matrix) is written to file 50.
!
!   Format of the files may depend on the edge : for s-edges, there is only one spectrum,
!   while for other edges, there are two spectra (j=l+1/2 and j=l-1/2).
!
!   A default input file case.inb for the broadening program is written.

! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!   cosmetics March 2005 (Kevin Jorissen)
!EOP

    implicit none

!   LOCAL VARIABLES
      integer L1,J,lauf
      real*8 matrix(0:jemax,0:lmax)
      real*8,external :: coreholelifetime 


      matrix(:,:)=dble(0)
!   Write case.inb for the broadening program :
      write(8,'(a80)') title
      write(8,'(a)') 'ELNES'
!!   In the next lines, we write 'case.elnes'
!      do j=1,80
!        if(fileroot(j:j).eq.'.') exit
!      enddo
!      if(j.lt.76) then
!        fileroot(j+1:j+5)='elnes'
!        if(j.lt.75) fileroot(j+6:80)=' '
!      endif

!      write(8,'(a80)') fileroot

!  The next lines tell the broadening program which columns to use for
!  broadening.  This depends on the type of calculation. 
      if(lc.gt.0) then
        write(8,'(i3,i3,i3)') 3,2,3
        write(8,'(3(g18.10,x))') split,dble(1),dble(1)
      else
        write(8,'(i3,i3,i3)') 1,1,0
        write(8,'(3(g18.10,x))') dble(0),dble(1),dble(0)
      endif
!   Next we determine the first energy at which the spectrum differs from zero.
!   This is important for energy dependent final state broadening.
      do j=0,jemax
        if (dabs(x(j,1)).gt.dble(0)) exit
      enddo
      if(j.gt.jemax) j=0  ! this only happens when the complete spectrum equals zero ...
      write(8,*) coreholelifetime(zz(natom),nc,lc,0),coreholelifetime (zz(natom),nc,lc,1)     ! core hole lifetime of two edges
      write(8,*) 1,ene(j)  ! 1 means linearly energy dependent valence broadening; ene(j) is edge offset (energy at which valence broadening should start)
      write(8,*) specbroad ! default spectrometer broadening of Gaussian FWHM 0.5 eV
      write(8,'(g18.10,/,g18.10,/,g18.10)') dble(0),dble(0),dble(0)  ! these parameters will not be used
      close(8)



!  Write case.elnes  
      if (headers) then    
!  Start with a large header containing important parameters.
      write(47,'(a)') '### The averaged spectrum.'
      write(47,'(a)') '## Parameters of the calculation :'
      write(47,4714) NATOM,NEQAT,NC,LC
      write(47,4715) DeltaE, SPLIT
      write(47,4716) Energy/1000
      write(47,4717) ThetaX*1000, ThetaY*1000
      write(47,4720) selrule,selrule2
      write(47,4724) SpecAperture*1000,ConvergAngle*1000
      write(47,4725) NR, NT
 4714 format('# Atom ',i3,',',i3,'edge   n=',i2,' l=',i2)
 4715 format('# Edge onset ',f6.1, ' eV;  split ',f6.2,' eV')
 4716 format('# Beam energy ',f12.0,' keV')
 4717 format('# Detector position ThetaX ',f6.2,' ThetaY ',f6.2,'  mrad')
 4720 format('# Selection rules ',a,2x,a)
 4724 format('# Collection semiangle ',f6.2,' Convergence semiangle ',f6.2,'  mrad')
 4725 format('# NR ',i3,' NT ',i3)
      endif

      write(47,*)
      if (lc.gt.0) then
        if (headers) write(47,*) '##  Energy, spectrum, partial spectrum 1, partial spectrum 2.'
        do J=0,JEMAX
          write(47,4713) ENE(J),X(J,1),X(J,2),X(J,3)
        enddo
      else
        if (headers) write(47,*) '##  Energy, spectrum.'
        do J=0,JEMAX
          write(47,4713) ENE(J),X(J,1)
        enddo
      endif

      if (verbosity.ge.1) then
        if (headers) write(49,*) '## Partial spectra.'
        if (lc.eq.0) then
          if (headers) write(50,*) '## Matrix Elements (ie, spectrum / dos).'
          lauf=1
        else
          if (headers) write(50,*) '## Matrix Elements (ie, spectrum / dos) for first edge.'
          lauf=3
        endif

        do j=1,jemax
          do l1=0,lmax
            if(dos(j,l1).lt.0.000001) then
              matrix(j,l1)=dble(0)
            else
              matrix(j,l1)=sdlm(l1,j,lauf)/dos(j,l1)
            endif
          enddo
        enddo

        if (headers) write(49,4811)
        if (headers) write(50,4811)
        do j=1,jemax
          write(49,4814) ENE(J),(SDLM(lauf,J,1),lauf=0,lmax)	! partial spectra
          write(50,4814) ene(j),(matrix(j,lauf),lauf=0,lmax)    ! matrix elements
        enddo

      endif   !  verbosity > 0


      RETURN
      
 4713 format(F10.5,5e14.5)
 4811 format('#',2x,'ENE',7x,'l=0  ',9x,'l=1  ',9x,'l=2  ',9x,'l=3   ')
 4814 format(F10.5,9e14.5)

      END



