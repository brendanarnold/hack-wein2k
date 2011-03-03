!BOP
! !ROUTINE: OrientedAngularSpectrum
! !INTERFACE:
      SUBROUTINE OrientedAngularSpectrum
! !USES:
      use input, only : Energy, DeltaE, LC, w1, w2, split
	  use initialstate1
	  use initialstate2
	  use energygrid
	  use spectra_hyperfine
! !DESCRIPTION:
!     This is an interface for the calculation of an orientation sensitive,
!     impuls transfer resolved spectrum.  Depending on the core state, it calls
!     AngularXSpectrum once or twice to calculate each separate edge, and
!     if necessary merges the two spectra into one total spectrum.
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none

!     first contribution (j = l+1/2)
      CALL AngularXSpectrum (DP, DQ, Energy-DeltaE, w1,1)
!     eventually second contribution (j = l-1/2)
      IF (LC.NE.0) THEN
         CALL AngularXSpectrum (DPdp, DQdq, Energy-DeltaE-SPLIT, w2,2)   ! used to be +SPLIT
!        Calculate the total theoretical spectrum
!        first, we copy one partial spectrum from first to third column, and we shift the other partial spectrum :
         X(:,3)=X(:,1)
	   SDLM(:,:,3)=SDLM(:,:,1)
	   CTR(:,:,3)=CTR(:,:,1)
!        now make the sum :
	   X(:,1)=X(:,2)+X(:,3)
	   SDLM(:,:,1)=SDLM(:,:,2)+SDLM(:,:,3)
	   CTR(:,:,1)=CTR(:,:,2)+CTR(:,:,3)
      ENDIF
	  if(jemax.gt.1) then
!   Multiply by energy step to get correct result for energy integration :
        X(:,:)=X(:,:)*(ENE(jemax)-ENE(jemax-1))
	    SDLM(:,:,:)=SDLM(:,:,:)*(ENE(jemax)-ENE(jemax-1))
        CTR(:,:,:)=CTR(:,:,:)*(ENE(jemax)-ENE(jemax-1))
	  endif

      
      RETURN
      END
      

