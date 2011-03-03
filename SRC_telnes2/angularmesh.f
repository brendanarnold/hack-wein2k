!BOP
! !ROUTINE: Angularmesh
! !INTERFACE:
      SUBROUTINE AngularMesh (ThXCenter, ThYCenter)
! !USES:
      use input,only : nr,nt,npos,thpart,convergangle,specaperture,th0
  	  use constants, only : pi
	  use qvectors,only : thxv,thyv
	  use program_control,only : qmodus
! !INPUT/OUTPUT PARAMETERS:
!     ThXCenter,ThYCenter : position of the detector in mrad
! !DESCRIPTION:
!     Creates the angular mesh : the detector is placed at 'Center'.
!     Now a circle of radius (collection + convergence semiangle) around
!     this center is sampled with NPos points.
!     Precise sampling depends on modus selected - Uniform, Logarithmic
!     or 1-dimensional.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP



!     NT points at the distance ThPart of the center, 3*NT at the distance
!     3*ThPart, ... , (2*NR-1)*NT at the distance (2*NR-1)*ThPart.
!     The coordinates of the center are (ThetaX, ThetaY).
!     The coordinates of the point i are (ThXV(i), ThYV(i)).
	  implicit none
!   INPUT : Position of the center of the aperture:
      real*8,intent(in) :: ThXCenter, ThYCenter
!   LOCAL VARIABLES :
      integer IRay, ITour, IndexPos, NPresentTour
      real*8 InterAngle,dxx
      

      if(qmodus.eq.'L'.or.qmodus.eq.'1') then
          dxx= dlog((convergangle+specaperture)/th0)/dble(NR-1)
      endif
      IF (npos.gt.1) THEN
         IndexPos = 0
         DO IRay = 1, NR
            NPresentTour = NT*(2*IRay-1)
            if(qmodus.eq.'1') NPresentTour=1
            InterAngle = dble(2) * PI / DBLE (NPresentTour)
            DO ITour = 1, NPresentTour
               IndexPos = IndexPos + 1
               if (qmodus.eq.'L'.or.qmodus.eq.'1') then
                 if(IRay.eq.1) then
                    ThXV(IndexPos) = ThXCenter+Th0*DCOS(InterAngle*ITour) &
                      /dble(2)
                    ThYV(IndexPos) = ThYCenter+Th0*DSIN(InterAngle*ITour) &
                      /dble(2)
                 else
                    ThXV(IndexPos) = ThXCenter+Th0*DCOS(InterAngle*ITour)* &
                      dexp(dxx*dble(IRay-2))*(dble(1)+dexp(dxx))/dble(2)
                    ThYV(IndexPos) = ThYCenter+Th0*DSIN(InterAngle*ITour)* &
                      dexp(dxx*dble(IRay-2))*(dble(1)+dexp(dxx))/dble(2)
                 endif
               else   
                 ThXV(IndexPos) = ThXCenter +  &
                    DCOS(InterAngle*ITour) * dsin(Thpart*(2*IRay-1))
                 ThYV(IndexPos) = ThYCenter +  &
                    DSIN(InterAngle*ITour) * dsin(Thpart*(2*IRay-1))
               endif
            ENDDO
         ENDDO      
      ELSE
         ThXV(1) = ThXCenter
         ThYV(1) = ThYCenter 
      ENDIF
      RETURN
      END


