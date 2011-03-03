!BOP
! !ROUTINE: WriteOrientedAngularSpectrum
! !INTERFACE:
    subroutine WriteOrientedAngularSpectrum
! !USES:
    use spectra_hyperfine
    use dimension_constants
    use input
    use energygrid
    use leftovers
    use qvectors,only : qv
    use constants
    use rotation_matrices,only : generalm
    use program_control,only : verbosity, relatQ, headers
    use crash
! !DESCRIPTION:
!   The angular differential cross section (in the array X) is written to file 47.
!   Depending on verbosity, the following files may be written :
!   Its partial l-components (array SDLM) are written to file 49.
!   Some of its cross term components (array CTR) are written to file 48.
!   Routines Writeangulardependence 1 and 2 are called for additional output.
!
!   Format of the files may depend on the edge : for s-edges, there is only one spectrum,
!   while for other edges, there are two spectra (j=l+1/2 and j=l-1/2).
!
!   This routine differs from WriteAveragedAngularSpectra in only three points :
!   ** It writes complex variables X and CTR
!   ** It has an additional file case.ctr .
!   ** Calls to Writeangulardependence 1 and 2.

! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!   cosmetics March 2005 (Kevin Jorissen)
!EOP

      implicit none

      integer k,lauf
      real*8 QVs(3,NPos),QVc(3,NPos),aam(3,3),qk(3),qk2(3),beta
! qvs : spherical coordinates of scattering vector
! qvc : cartesian


!  First write an extensive header in case.elnes so that all the important
!  parameters are easily accessible even if case.innes got lost.

      if (headers) then
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

!  The main program uses the Q-vectors in the crystal frame, since this is how
!  the DOS is expressed.  However, if we want to study the angular behaviour of
!  scattering, we almost certainly want output in the laboratory frame.
!  Therefore, we need to transform back using the inverse of the GeneralM
!  matrix (stored in the matrix aam).  Additionally, we remove the relativistic
!  correction in qz, since we actually need output as a function of k', the wave
!  vector after scattering.

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
        if (headers) write(47,*) '##  Index, Qx, Qy, Qz, Q, Qtheta, Qfi, Re [spectrum], Im  [spectrum], Re ',&
               '[partial spectrum 1], Im [partial spectrum 1], Re [partial spectrum 2], Im [partial spectrum 2].'
        do k=1,NPos
          write(47,4713) k,QVc(1,k),QVc(2,k),QVc(3,k),QVs(1,k), &
           QVs(2,k),QVs(3,k),X(k,1),X(k,2),X(k,3)
        enddo
      else
        if (headers) write(47,*) '##  Index, Qx, Qy, Qz, Q, Qtheta, Qfi,  Re [spectrum], Im  [spectrum].'
        do k=1,NPos
          write(47,4713) k,QVc(1,k),QVc(2,k),QVc(3,k),QVs(1,k),QVs(2,k),QVs(3,k),X(k,1)
        enddo
      endif


      if (verbosity.ge.1) then
!------------header for crossterm files - case.ctr
        if (headers) write(48,*) '### The cross terms. (Written as complex data)'
        if (headers) write(48,4810)
!------------header for partial spectrum files - case.sdlm
        if (headers) write(49,*) '### The partial spectra. (Written as real data)'
        if (headers) write(49,4811)
        do k=1,npos
          write(49,4814) k,QVc(1,k),QVc(2,k),QVc(3,k),QVs(1,k), &
          QVs(2,k),QVs(3,k),(SDLM(lauf,k,1),lauf=1,10)	! partial spectrum
          write(48,4813) k,QVc(1,k),QVc(2,k),QVc(3,k),QVs(1,k), &
          QVs(2,k),QVs(3,k),(CTR(lauf,k,1),lauf=1,12)	! Crossterm
        enddo
      endif

      if (verbosity.ge.2) then
        call errflg(errfn,'Error in PlotAngularDependence')
        call WriteAngularDependence1(qvs)
        call WriteAngularDependence2
! careful ! after writeangulardependence2, the list of Q-vectors does not exist anymore ...
! fortunately, this is the end of the program.
        write(6,'(a)') 'Q-vectors destroyed by WriteAngularDependence2.'
      endif


 4713 FORMAT(I6,12e14.5)
 4813 FORMAT(I6,30e14.5)
 4814 FORMAT(I6,16e14.5)
 4810 FORMAT('#',2x,'Index, Qx, Qy, Qz, Q, Qtheta, Qfi, ',7x,'00 22',23x,'00 21',23x,'00 20',23x, &
      	'11 1-1',22x,'11 10',23x,'22 2-2',22x,'22 21',23x, &
      	'22 2-1',22x,'21 2-1',22x,'20 22',23x,'20 21',23x,'all other cross term contributions')
 4811 FORMAT('#',2x,'Index, Qx, Qy, Qz, Q, Qtheta, Qfi, ',7x,'00   ',9x,'1-1  ',9x,'10   ',9x, &
      	'11    ',8x,'2-2  ',9x,'2-1   ',8x,'20   ',9x, &
      	'21    ',8x,'22    ',9x,'sum{m} 3m')
      RETURN
      END

