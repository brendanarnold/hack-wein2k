!BOP
! !ROUTINE: WriteAngularDependence2
! !INTERFACE:
    subroutine WriteAngularDependence2
! !USES:
    use input
    use spectra_hyperfine
    use qvectors,only : WeightV,destroy_qvecs
    use program_control,only : qmodus, headers
    use constants,only : pi
! !DESCRIPTION:
!   The angular differential partial cross sections are integrated up to a variable
!   collection angle.  The resulting partial intensities are written, as a function
!   of collection semiangle, to file 59.

!   Only the first (j=l+1/2) edge is used.
!
!   Be careful when using this routine, since it will destroy the mesh of
!   Q-vectors !!!
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP

      implicit none

!  LOCAL VARIABLES
      integer ncol,j,i,nposi,nri
      real*8 WeightW(npos),betastep,beta,dx
      real*8, allocatable :: ints(:,:),sp2(:),collection(:)


      write(6,'(/,a)') 'Output from subroutine PlotAngularDependence:'	

!  There are only beta/(alfa+beta) * NR collection angles !!!
      if (qmodus.eq.'L'.or.qmodus.eq.'l'.or.qmodus.eq.'1') then
        dx=dlog((SpecAperture+ConvergAngle)/th0)/dble(NR-1)
        ncol=1+int(log(SpecAperture/th0)/dx)
      else
        betastep=(ConvergAngle+SpecAperture)/dble(NR)
        ncol=int(dble(SpecAperture/(SpecAperture+ConvergAngle))*dble(NR))
      endif

!   We first save old data :
      beta=specaperture
      nposi=npos
      nri=nr
      WeightW(:)=WeightV(:)


      write(6,'(a,i6,a)') 'Using ',ncol,' collection semi-angles.'
      allocate(ints(ncol,6),sp2(ncol),collection(ncol))
      if(qmodus.eq.'U') then
        write(6,*) 'Using betastep = ',betastep
      else
        write(6,*) 'Using dx = ',dx
      endif
      ints(:,:)=dble(0)
      sp2(:)=dble(0)
      collection(:)=dble(0)
      if (headers) write(59,'(a)') '#    beta        sp2        pi           sigmadip        total      ',&
               '   monopole     quadrupole     octopole'

      do i=1,ncol
!  Loop over collection angles :      
        if (qmodus.eq.'L'.or.qmodus.eq.'l'.or.qmodus.eq.'1') then
          if(i.eq.1) then
            specaperture=th0
          else
            specaperture=th0*dexp(dble(i-1)*dx)
          endif
        else
          specaperture=betastep*i
        endif
        collection(i)=specaperture
        if(qmodus.eq.'1') then
          npos=i
        else
          npos=i*i*nt
        endif
        nr=i
        ThPart = (ConvergAngle + SpecAperture) / DBLE(2*NR)

!  For every collection angle, new integration weights have to be calculated :        
        call destroy_qvecs
        call CalculateWeights
        if(npos.eq.1) weightv(1)=pi*(dble(1000)*(convergangle+specaperture)*min(convergangle,specaperture)/convergangle)**2
        if(npos.eq.1) write(6,*) weightv(1)
!  Normally, we don't want weights for just one position - since, then, we feel the user does not want integration, but the DOUBLE dscs.
!  However, here we do need integration.  So we make the weight ourselves.

!  We integrate all contributions to the spectrum separately :
        do j=1,npos
          ints(i,1)=ints(i,1)+weightv(j)/weightw(j)* sdlm(3,j,1)   ! pi
          ints(i,2)=ints(i,2)+weightv(j)/weightw(j)*(sdlm(4,j,1)+sdlm(2,j,1))  ! sigma dipole
          ints(i,4)=ints(i,4)+weightv(j)/weightw(j)* sdlm(1,j,1)   ! monopole
          ints(i,5)=ints(i,5)+weightv(j)/weightw(j)*(sdlm(9,j,1)+sdlm(8,j,1)+sdlm(7,j,1)+sdlm(6,j,1)+sdlm(5,j,1))  ! quadrupole
          ints(i,6)=ints(i,6)+weightv(j)/weightw(j)* sdlm(10,j,1)  ! octopole
        enddo !j
        ints(i,3)=ints(i,1)+ints(i,2)+ints(i,4)+ints(i,5)+ints(i,6) ! total
        if(dabs(ints(i,3)).gt.0) sp2(i)=ints(i,1)/ints(i,3)         ! sp2 = pi/total

        write(59,'(8(F14.9,x),I7)') collection(i),sp2(i),(ints(i,j),j=1,6),npos

      enddo  ! i


!   we restore old data :
      specaperture=beta
      npos=nposi
      nr=nri
      ThPart = (ConvergAngle + SpecAperture) / DBLE(2*NR)
      call destroy_qvecs
      call CalculateWeights
!   clean up :
      deallocate(ints,sp2,collection)


      return
      end
