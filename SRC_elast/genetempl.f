!******* Generates rhomb.templ, tetra.templ from init.struct.
!*******  Those files will be used for initializing rhombohedral
!*******  and tetragonal strain calculations (i.e. set up proper 
!*******  symetries, kmesh...). They are 0.5% strained.
!*******  
!*******   Author: CHARPIN Thomas
!*******           Laboratoire des Geomateriaux de l'IPGP
!*******           4,pl Jussieu
!*******           F-75252 Paris Cedex 05
!*******           France
!*******           charpin@ipgp.jussieu.fr   
!*******


      PROGRAM genetempl
      CHARACTER*79 row
      CHARACTER*23 goodi
      CHARACTER*4 latt
      CHARACTER*1 lat
      INTEGER noa,i,indic,j
      REAL*8 a,b,c,alp,bet,gam
      REAL*8 ah,ch,at,ct,eps
      REAL*8 x,y,z,xr,yr,zr
      EXTERNAL c2rh,c2te,xc2rh

      PARAMETER(eps=5.d-3)

      OPEN(unit=10,file='init.struct',status='old',err=8000)
      OPEN(unit=20,file='rhomb.templ')
      OPEN(unit=30,file='tetra.templ')

      READ(10,1000) row
      WRITE(20,1000)row
      WRITE(30,1000)row

      READ(10,3010) latt,goodi,noa
      WRITE(20,1010) 'R   ',goodi,noa
      WRITE(30,1010) latt,goodi,noa
      lat=latt(1:1)

      READ(10,1000) row
      WRITE(20,1000)row
      WRITE(30,1000)row

      READ(10,1020) a,b,c,alp,bet,gam
      CALL c2rh(a,eps,lat,ah,ch)
      CALL c2te(a,eps,at,ct)
      WRITE(20,1020) ah,ah,ch,alp,bet,gam
      WRITE(30,1020) at,at,ct,alp,bet,gam

      DO i=1,noa
         READ(10,1030) indic,x,y,z
         WRITE(30,1040) indic,x,y,z
         CALL xc2rh(x,y,z,lat,xr,yr,zr)
         WRITE(20,1040) indic,xr,yr,zr
         read(10,'(15x,i2)') mult
         write(30,'(15x,i2)') mult
         write(20,'(15x,i2)') mult
          do j=1,mult-1
         READ(10,1030) indic,x,y,z
         WRITE(30,1040) indic,x,y,z
         CALL xc2rh(x,y,z,lat,xr,yr,zr)
         WRITE(20,1040) indic,xr,yr,zr
          enddo 
         DO j=1,4
            READ(10,1000) row
            write(20,1000) row
            write(30,1000) row
         ENDDO

      ENDDO
      WRITE(20,"('   0 SYMMETRY OPERATIONS:')")
      WRITE(30,"('   0 SYMMETRY OPERATIONS:')")
      

      CLOSE(30)
      CLOSE(20)
      CLOSE(10)


 1000 FORMAT(A79)
 3010 FORMAT(A4,A23,I3)
 1010 FORMAT(A4,A23,I3)
 1020 FORMAT(6F10.6)
 1030 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8)
 1040 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)


      GOTO 8001
 8000 write(*,*) 'No valid file init.struct'

      STOP

 8001 write(*,*) 'Template files generated...'

      END
