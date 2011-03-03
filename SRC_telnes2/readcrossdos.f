!BOP
! !ROUTINE: ReadCrossDOS
! !INTERFACE:
      SUBROUTINE ReadCrossDOS

! !USES:
      use dimension_constants,only : lmax,nofcross
      use energygrid
	  use cross_DOS
! !DESCRIPTION:
!     READ the DD coefficients from case.xdos.
!     Calculate l,m DOS.
!
!     xdos(nofcross) : Dlm DLM*
!     $I1 = (l+1)^2 - l + m$
!     $I2 = (L+1)^2 - L + M$
!     it means that:
!     
!     l m     I
!     ---------
!     0 0  -> 1
!     1-1  -> 2
!     1 0  -> 3
!     1 1  -> 4
!     2-2  -> 5
!          :
!     
!     so that:
!     Dlm DLM* = xdos( I1*(I1-1)/2 + I2 )
!      
!       I2 1  2  3  4  5  6 
!     I1 
!       1  1
!       2  2  3
!       3  4  5  6 
!       4  7  8  9  10
!       :  11 12 ...
! !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!EOP

      implicit none
!  LOCAL VARIABLES
      CHARACTER*2 CH
      CHARACTER*7 FNAME
      INTEGER I, J, ITAP, IDD
      INTEGER el, em, I1

!     Read partial (cross) density of states from *.xdos :
      READ (58,4712) IEMAX

	  call make_ene(iemax)
      call make_crossdos(iemax,nofcross,lmax)

      do i=1,iemax
	    read(58,*,end=913,err=913) ene(i),(xdos(i,j),j=1,nofcross)
	  enddo

!     Calculate the partial DOS:
      CALL DefineIndex
      DO I=1, IEMAX
         DO el = 0, lmax  ! used to be 3
            DOSlm(I, el) = dble(0)         !!! doslm used to be called dos
            DO em = -el, el
               I1 = (el+1)**2-el+em
               DOSlm(I, IndexLM(el, em)) = dble(xdos(i, ((I1+1)*I1)/2))  != DDR(I, ((I1+1)*I1)/2)
               IF (el.NE.0) THEN
!     if el = 0, IndexLM (0,0) = 0 = el ...
                  DOSlm(I, el) = DOSlm(I, el) + DOSlm(I, IndexLM(el,em))
               ENDIF
            ENDDO
            IF (el.EQ.3) THEN   ! since apparently we still do not have charge splitting for f ??
			   do em=-el,el
                 DOSlm(I, IndexLM(3, em)) = DOSlm(I, 3) / dble(7)
			   enddo
            ENDIF
         ENDDO
      ENDDO
      

      RETURN
 911  WRITE(6,*) 'ERROR IN OPENING ', FNAME
      stop
 912  WRITE(6,*) 'ERROR IN READING ', FNAME
      stop
 913  WRITE(6,*) 'ERROR WHILE READING XDOS FILE.'
      stop

 4712 FORMAT(/,38x,I5,/)
 4713 FORMAT(F10.5,7e14.5)
      
      END
