!*********************************************************
!********* INPUT: Cubic coordinates of atom x,y,z
!*********        Lattice type lat (P,B or F)
!********* OUTPUT: Rhombohedral coordinates xr,yr,zr 
!*********         
!*********
!*********  
!*********




      SUBROUTINE xc2rh(x,y,z,lat,xr,yr,zr)
      CHARACTER*1 lat
      REAL*8 x,y,z,xr,yr,zr



      IF (lat.eq.'F') THEN

!********
!******** a     -1  1  1    ar
!******** b  =   1 -1  1 .  br
!******** c      1  1 -1    cr
!********
         
         xr=-x+y+z
         yr=x-y+z
         zr=x+y-z
         GOTO 5000
      ENDIF


      IF (lat.eq.'B') THEN

!********
!******** a      0  1  1    ar
!******** b  =   1  0  1 .  br
!******** c      1  1  0    cr
!********

         xr=y+z
         yr=x+z
         zr=x+y

         GOTO 5000
      ENDIF


      IF (lat.eq.'P') THEN
!******** ar=a br=b cr=c

         xr=x
         yr=y
         zr=z

         GOTO 5000
      ENDIF



      
      WRITE(*,*) 'Problem: your structure is neither P, B or F'
      STOP

 5000 xr=xr-int(xr)
      yr=yr-int(yr)
      zr=zr-int(zr)
      IF (xr.lt.0.d0) xr=xr+1.d0
      IF (yr.lt.0.d0) yr=yr+1.d0
      IF (zr.lt.0.d0) zr=zr+1.d0

      RETURN
      END
