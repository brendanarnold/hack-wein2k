!BOP
! !ROUTINE: CalculateWeights
! !INTERFACE:
      SUBROUTINE CalculateWeights
! !USES:
      use qvectors
	  use constants,only : pi
	  use input
	  use struct,only : nats,nat,mult
	  use program_control,only : qmodus
! !DESCRIPTION:
!  Set up the mesh of (k,k') - pairs for which the double differential scattering cross-section
!  will be calculated.
!  The mesh is determined by collection and convergence semiangle, and by the type of mesh
!  chosen (Uniform, Logarithmic, 1-dimensional) and its specified size (Nr, Nt, NPos).
!
!  A note : since we do not change the microscope and detector semiangles while we measure
!  a full edge, also the definition of the 'pixels' in our calculation, i.e. the *angular* mesh
!  used for sampling the Q vectors, should not change.
!  Since there is a relation connecting energy loss, the magnitude of impuls transfer, and the
!  scattering angle, keeping the latter fixed does imply that the vector Q and its length change
!  when energy changes.
!  Therefore, we run calculweight and newthxthy only once at the beginning of the calculation,
!  but we run coeffsqcf for every energy in the edge.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP



	  implicit none
!   LOCAL variables
      integer IRay, ITour, IndexPos, NPresentTour
      real*8 ConvolValue, Lfactor, dxx
      real*8 Theta,sa,ca,p

	  write(6,'(/,a)') 'Output from subroutine CalculWeight:'
	  write(6,*) specaperture,convergangle
	  write(6,'(a)') 'ThXV          ThYV       weight       overlap_factor'
      
      call make_Qvecs1(NPos)
      call make_Qvecs2(NPos,mult(natom))
      if(qmodus.eq.'L'.or.qmodus.eq.'1') then
          dxx= dlog((convergangle+specaperture)/th0)/dble(NR-1)
      endif
      CALL AngularMesh (dble(0),dble(0))
      IF (npos.gt.1) THEN

	     sa=specaperture    ! for convenience, we use these abbreviations :
		 ca=convergangle
         IndexPos = 0
         DO IRay = 1, NR
            if(qmodus.eq.'1') then
              NPresentTour=1
            else
              NPresentTour = NT*(2*IRay-1)
            endif
            Theta = ThXV(IndexPos+NPresentTour) ! The last position in a ring has cos(2 pi) = 1; and here ThXCenter = 0

            IF     (Theta.LE.dabs(sa-ca)) THEN
			   if (ca.gt.dble(0.000001).and.sa.gt.dble(0.000001)) then
                 ConvolValue = pi* min(ca,sa)**2  / (pi*ca**2)
			   else
			     ConvolValue = dble(1)
			   endif
            ELSEIF (Theta.GE.(sa+ca)) THEN
                  ConvolValue = dble(0)
            ELSE  ! |a-b| < Theta < a+b
                  p=(Theta**2+ca**2-sa**2)/(dble(2)*Theta)
				  convolvalue=pi/dble(2)*(ca**2+sa**2)-p*dsqrt(ca**2-p**2) &
				   -(Theta-p)*dsqrt(sa**2-(Theta-p)**2)-sa**2*dasin((Theta-p)/sa) &
				   -ca**2*dasin(p/ca)
				  convolvalue = convolvalue / (pi*ca**2)   ! if ca=0, we wouldn't be here !
            ENDIF

            DO ITour = 1, NPresentTour
              IndexPos = IndexPos + 1
!             Weight is the surface around a point multiplied with the value at the
!             point of the convolution of SpecAperture and Convergence distribution.
              WeightV(IndexPos) = (1.D6*ThPart**2)/ DBLE(NPresentTour) &
                    * PI * dble(4) * (2 * IRay - 1) * ConvolValue
!             Now correct the weights for nonuniform Q-meshes :
              if(qmodus.eq.'L'.or.qmodus.eq.'1') then
                if (Iray.eq.1) then
                  Lfactor=(dble(NR)*Th0/(sa+ca))**2*NT/NPresentTour
                else
                  Lfactor=(dble(NR)*Th0*dexp(dxx*(iray-2)) &
                    /(sa+ca))**2*(dexp(2*dxx)-dble(1))*NT/NPresentTour
                endif
              else
                Lfactor=dble(1)
              endif
              WeightV(IndexPos)=WeightV(IndexPos)*Lfactor
              WRITE(6,'(4(f10.5,2x))') ThXV(IndexPos), ThYV(IndexPos),  &
                    WeightV(IndexPos), Convolvalue

            ENDDO
         ENDDO
      ELSE
!     Only one spectrum
         WeightV(1) = dble(1)
         write(6,'(4(f10.5,2x))') ThXV(1), ThYV(1),WeightV(1),dble(1)
      ENDIF

!     for the rest of the calculation:
      call AngularMesh(thetax,thetay)  ! makes thx and thy for given alfa, beta
      RETURN
      END
      

