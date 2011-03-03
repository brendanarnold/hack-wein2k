      DOUBLE PRECISION FUNCTION D1MACH(I)
      implicit real*8 (a-h,o-z)
!C
!C  DOUBLE-PRECISION MACHINE CONSTANTS
!C
!C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!C  D1MACH( 5) = LOG10(B)
!C
!c     use values in float.h to define these
      save dmach
!C
      DOUBLE PRECISION DMACH(5),values(10)
!C
      call mdble(values)
      dmach(1)=values(2)
      dmach(2)=values(3)
      dmach(3)=values(1)
!c     take largest and smallest spacings as equal
      dmach(4)=values(1)
!c     we do not use dmach(5) -- so ignore
      if(i.eq.5)then
           write(6,*)'Need to set dimach(5) !'
           write(6,*)'Dying....'
           call exit(0)
      endif
!C
      IF (I .LT. 1  .OR.  I .GT. 5) Write(6,*)'D1MACH - I OUT OF BOUNDS'
      D1MACH = DMACH(I)
      RETURN
!C
      END
