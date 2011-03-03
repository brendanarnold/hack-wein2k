!********************************************************
!******** INPUT: Cubic lattice parameter a                                   
!********        Strain factor eps (ch=ch0(1+eps))
!********        Lattice type lat (P,B or F)
!******** OUTPUT: hexagonal lattice parameters for 
!********         the R (strained) lattice
!********
!********  ah=ar-br     ch=ar+br+cr
!********




      SUBROUTINE c2rh(a,eps,lat,ah,ch)
      CHARACTER*1 lat
      REAL*8 a,eps
      REAL*8 ah,ch,ch0




      IF (lat.eq.'F') THEN
!****** ar=1/2(b+c) br=1/2(a+c) cr=1/2(a+b) alphar=60
         ah=a/(2.d0**(0.5d0))
         ch0=a*(3.d0**(0.5d0))
         ch=ch0*(1.d0+eps)

         GOTO 5000
      ENDIF


      IF (lat.eq.'B') THEN
!****** ar=1/2(-a+b+c) br=1/2(a-b+c) cr=1/2(a+b-c) alphar=109°28'
         ah=a*(2.d0**(0.5d0))
         ch0=a*(3.d0**(0.5d0))/2.d0
         ch=ch0*(1.d0+eps)

         GOTO 5000
      ENDIF


      IF (lat.eq.'P') THEN
!****** ar=a br=b cr=c alphr=90
         ah=a*(2.d0**(0.5d0))
         ch0=a*(3.d0**(0.5d0))
         ch=ch0*(1.d0+eps)

         GOTO 5000
      ENDIF



      
      WRITE(*,*) 'Problem: your structure is neither P, B or F'
      STOP

 5000 RETURN
      END
