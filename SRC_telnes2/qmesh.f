!BOP
! !ROUTINE: QMesh
! !INTERFACE:
      SUBROUTINE QMesh (Energy2)
! !USES:
      use program_control,only : verbosity, RelatQ
      use input, only : npos, neqatdeb, neqat, energy, deltae, natom, k0len
	  use qvectors
	  use constants
      use rotation_matrices
! !INPUT/OUTPUT PARAMETERS:
!     Energy2 : Energy of the beam electron after scattering (in eV).
! !DESCRIPTION:
!  The impuls transfer vector is calculated for one particular k0
!  and for a set of k' at energy Energy2 and corresponding to different
!  scattering directions, as expressed by ThXV and ThYV.
!  The Q-mesh will be used for evaluation and integration of the cross-section
!  over all initial and final states of the beam electron permitted by
!  collection and convergence semiangles.
!
!  A relativistic correction is added to Q (=k0-k'-correction).

!  A little output is written to the master output file 6.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none
	      
!   INPUT (energy in eV):
      real*8,intent(in) ::  Energy2
!   LOCAL variables:
      real*8  Qrot(3), QLen, Q(3), KPrLen,beta,QLenClas
      integer IRay, ITour, IAtom, IPos, NPresentTour
      real*8, external :: WaveLength

      KPrLen = 2*PI/WaveLength(Energy2)
!     the vectors K0 and KPr are:
!     K0 = K0Len*(0, 0, -1)
!     KPr = KPrLen * (sin (ThX), sin(ThY), -SQRT(1 - sin^2(ThX) - sin^2(ThY)) )
!     Q = K0 - KPr 

      do IPos=1,NPos
!     The Q vectors are given in the basis of the observator:
            Q(1) = - KPrLen * ThXV(IPos)
            Q(2) = - KPrLen * ThYV(IPos)
            Q(3) =   KPrLen * DSQRT(1-ThXV(IPos)**2-ThYV(IPos)**2) - K0Len

            QLenClas = DSQRT(Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3))
			beta=dsqrt((2+Energy/MeC2)/(2+Energy/MeC2+MeC2/Energy))
			if(RelatQ) Q(3)=Q(3)+(energy-energy2)/hbarc * beta
!  Note the sign : this is because of the stupid convention where the z-axis is antiparallel to the incident beam.

            QLen = DSQRT(Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3))

            DO IAtom = NEqAtDeb, NEqAt
!     For every equivalent atom of the unit cell, we put the Q vectors
!     in the local basis of the atom: Q -> Qrot 
!	  IEqAtom numbers the atoms in the equivalency class (as defined in case.struct).

               CALL ProductMatVect(GeneralM(1,1,IAtom),Q,Qrot)

               QV(1:3, IAtom, IPos) = Qrot(1:3)
               QLenV(IAtom, IPos)   = QLen
			   QLenVClas(IAtom,IPos) = QLenClas

            ENDDO
      enddo

      if(verbosity.ge.2.and.(dabs(energy2-(energy-deltae)).lt.dble(0.0001))) then
! print this only for the very first energy point (if present), otherwise file gets too large.
        write(6,'(/,a)') 'Output of subroutine CoeffsQCV:'
	    write(6,'(a,i4,x,i4)') '# Qx      Qy        Qz  for atom ',natom,neqatdeb
		do ipos=1,npos
		  write(6,'(3(f10.5,x))') qv(1:3,neqatdeb,ipos)
		enddo
	  endif

      RETURN
      END
