!BOP
! !ROUTINE: WriteAveragedAngularSpectrum
! !INTERFACE:
    subroutine WriteAveragedAngularSpectrum
! !USES:
    use dimension_constants, only : lmax
    use rotation_matrices, only : GeneralM
    use energygrid
    use spectra_normal
    use constants
    use input
    use program_control,only : verbosity, relatQ, headers
    use leftovers
    use qvectors,only : qv
! !DESCRIPTION:
!   The angular differential cross section (in the array X) is written to file 47.
!   Depending on verbosity, the following files may be written :
!   Its partial l-components (array SDLM) are written to file 49.
!
!   Format of the files may depend on the edge : for s-edges, there is only one spectrum,
!   while for other edges, there are two spectra (j=l+1/2 and j=l-1/2).

! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP
      implicit none

      integer k,lauf
      real*8 QVs(3,NPos),QVc(3,NPos),aam(3,3),qk(3),qk2(3),beta  
! qvs : spherical coordinates of scattering vector
! qvc : cartesian
      

      if (headers) then
!  First write a header to case.elnes.
      write(47,'(a)') '### The averaged spectrum.'
      write(47,'(a)') '## Parameters of the calculation :'
      write(47,4714) NATOM,NEQAT,NC,LC
      write(47,4715) DeltaE, SPLIT
      write(47,4716) Energy/1000
      write(47,4717) ThetaX*1000, ThetaY*1000
      write(47,4720) selrule,selrule2
      write(47,4724) SpecAperture*1000,ConvergAngle*1000
      write(47,4725) NR, NT
4714  format('# Atom ',i3,',',i3,'edge   n=',i2,' l=',i2)
4715  format('# Edge onset ',f6.1, ' eV;  split ',f6.2,' eV')
4716  format('# Beam energy ',f12.0,' keV')
4717  format('# Detector position ThetaX ',f6.2,' ThetaY ',f6.2,'  mrad')
4720  format('# Selection rules ',a,2x,a)
4724  format('# Collection semiangle ',f6.2,' Convergence semiangle ',f6.2,'  mrad')
4725  format('# NR ',i3,' NT ',i3)
      endif


!  We need to transform Q back from crystal frame to laboratory frame.  Since we
!  also need to go back to k' (wave vector after scattering), we remove the
!  relativistic term which went into the calculation of the cross-section.

      call Inverse(GeneralM(1, 1,1),aam)
      do k = 1,NPos
        Qk(1:3)=qv(1:3,1,k)
        call ProductMatVect (aam, Qk,Qk2)
        qvc(1:3,k)=qk2(1:3)
        if (relatQ) then
          beta=dsqrt((2+Energy/MeC2)/(2+Energy/MeC2+MeC2/Energy))
          qvc(3,k)=qvc(3,k)-beta*(energy-ene(jemin))/hbarc
        endif
        QVs(1,k)=dsqrt(qvc(1,k)**2+qvc(2,k)**2+qvc(3,k)**2)

        if(dabs(qvc(3,k)).lt.0.00001) then
          QVs(2,k)=pi/dble(2)
        else
          QVs(2,k)=datan(dsqrt(qvc(1,k)**2+qvc(2,k)**2)/qvc(3,k))
        endif
        if (dabs(qvc(1,k)).lt.0.00001) then
          QVs(3,k)=pi/dble(2)
        else
          QVs(3,k)=datan(qvc(2,k)/qvc(1,k))
        endif
      enddo


!------------ Write case.elnes
      write(47,*)
      if (lc.gt.0) then
        if (headers) write(47,*) '##  Index, Qx, Qy, Qz, Q, Qtheta, Qfi, spectrum, partial spectrum 1, partial spectrum 2.'
        do k=1,NPos
          write(47,4713) k,QVc(1,k),QVc(2,k),QVc(3,k),QVs(1,k), &
           QVs(2,k),QVs(3,k),X(k,1),X(k,2),X(k,3)
        enddo
      else
        if (headers) write(47,*) '##  Index, Qx, Qy, Qz, Q, Qtheta, Qfi, spectrum.'
        do k=1,NPos
          write(47,4713) k,QVc(1,k),QVc(2,k),QVc(3,k),QVs(1,k), &
              QVs(2,k),QVs(3,k),X(k,1)
        enddo
      endif

      if (verbosity.ge.1) then
!------------ Write case.sdlm
        if (headers) write(49,*) '##  Partial spectra.'
        if (headers) write(49,4811)
        do k=1,NPos
          write(49,4713) k,QVc(1,k),QVc(2,k),QVc(3,k),QVs(1,k), &
           QVs(2,k),QVs(3,k),(SDLM(lauf,k,1),lauf=0,lmax)
        enddo
      endif


      RETURN
 4713 format(I6,14e14.5)
 4811 format('#','Index, Qx, Qy, Qz, Q, Qtheta, Qfi, ',7x,'l=0  ',9x,'l=1  ',9x,'l=2  ',9x,'l=3   ')

      END


