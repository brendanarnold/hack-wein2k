!BOP
! !ROUTINE: AveragedEnergySpectrum
! !INTERFACE:
      SUBROUTINE AveragedEnergySpectrum
! !USES:
      use dimension_constants
      use input
	  use spectra_normal
	  use energygrid
	  use qvectors
	  use radial_functions
      use initialstate1
	  use initialstate2
	  use densityofstates,only : dos
	  use leftovers
	  use program_control,only : relatq
	  use constants,only : hbarc
! !DESCRIPTION:
!     Calculate the part of the spectrum for every l' of the final state
!     Outputs: intensities XI(J,LL), XI2(J,LL) as a function of energy.
!
!     Depending on the core state, it calls
!     EnergyXSpectrum once or twice to calculate each separate edge, and
!     if necessary merges the two spectra into one total spectrum.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none

! LOCAL STAFF :
      integer Lambda,LL,J,IPos,IAtom,indx,signe
      real*8 EI,frac2,JQ,term,prefac


      DO LL=0, lmax
!     the up-limit (min ..) is normally LC+LambdaMax, but the DOS is
!     not calculated till such a LL (only to lmax).
      if (UseThisL(LL)) then    ! selection rules.

         DO J=JEMIN, JEMAX

!     only for positive energies:
            EI=ENE(J)
            IF (EI.GE.0) THEN
!     here we calculate the Q(E):   
               CALL QMesh (Energy-DeltaE-EI)
               DO IPos = 1, NPos
                  DO IAtom = NEqAtDeb, NEqAt
				     if(relatQ) then
				       prefac=(QLenVClas(IAtom,IPos)**2-((DeltaE+ei)/hbarc)**2)**2
					 else
					   prefac=QLenVClas(IAtom,IPos)**4
					 endif
                     DO Lambda=iabs(lc-ll),lc+ll
					 if(UseThisLambda(Lambda)) then
                        CALL threejm0(LC,Lambda,LL,Frac2,signe)
                        CALL RINT14(CFEIN1,CFEIN2,Uell(1,J,LL),UPell(1,J,LL), &
                             DP,DQ,JQ,QLenV(IAtom,IPos),Lambda,NATOM)
                        term = dble(2)*w1*gammaa0**2*DOS(J,LL)/prefac* Frac2 * &
                             WeightV(IPos)*JQ**2*DBLE((2*Lambda+1))*(2*LC+1)
                        X(J,1)=X(J,1) + term
!--------------------------different DOS contributions
                        SDLM(LL,J,1) = SDLM(LL,J,1) + term
					 endif  ! selection rule UseThisLambda
                     ENDDO  ! Lambda
                  ENDDO  ! IAtom
               ENDDO  ! IPos

!     Now for the j=l-1/2 case:    We need new Q-vectors !!
               IF (LC.GT.0.and.J.le.(jemax-jdecal)) THEN
!     here we calculate the Q(E):   
                CALL QMesh (Energy-DeltaE-EI-SPLIT)
                DO IPos = 1, NPos
                  DO IAtom = NEqAtDeb, NEqAt
				     if(relatQ) then
				       prefac=(QLenVClas(IAtom,IPos)**2-((DeltaE+Split+ei)/hbarc)**2)**2
					 else
					   prefac=QLenVClas(IAtom,IPos)**4
					 endif
                     DO Lambda=iabs(lc-ll),lc+ll
					 if(UseThisLambda(Lambda)) then
                        CALL threejm0(LC,Lambda,LL,Frac2,signe)
                           CALL RINT14(CFEIN1, CFEIN2, Uell(1,J,LL),UPell(1,J,LL), &
                                DPdp,DQdq,JQ,QLenV(IAtom,IPos),Lambda,NATOM)
                           term = dble(2)*w2*gammaa0**2*DOS(J,LL)/prefac* &
                                WeightV(IPos)*JQ**2*DBLE((2*Lambda+1))*(2*LC+1)* Frac2
                           X(J+jdecal,2)= X(J+jdecal,2) + term
!--------------------------different DOS contributions
                           SDLM(LL,J+jdecal,2) = SDLM(LL,J+jdecal,2) + term
					 endif  ! selection rule UseThisLambda
                     ENDDO  ! Lambda
                  ENDDO  ! IAtom
                ENDDO  ! IPos
               ENDIF  !  LC > 0
            ENDIF  ! triangle
         ENDDO  ! J       
	  endif ! selection rules.
      ENDDO ! LL


      IF (LC.GT.0) THEN
!        Calculate the total theoretical spectrum
!        First, we copy one partial spectrum from first to third column :
         X(:,3)=X(:,1)
		 SDLM(:,:,3)=SDLM(:,:,1)
!        now make the sum :
		 X(:,1)=X(:,2)+X(:,3)
		 SDLM(:,:,1)=SDLM(:,:,2)+SDLM(:,:,3)
      ENDIF

      
      RETURN
      END

