MODULE param
  INTEGER,PARAMETER :: IBLOCK= 128
  !c.....Optimize IBLOCK for your hardware (32-255)
  INTEGER,PARAMETER :: LMAX2=    8
  INTEGER,PARAMETER :: LOMAX=    3
  INTEGER,PARAMETER :: NCOM=   121
  ! for ncom parameter check format 1003 in l2main.frc
  INTEGER,PARAMETER :: NLOAT=    3
  INTEGER,PARAMETER :: NRAD=  1481
  INTEGER,PARAMETER :: NGAU=  2350
  ! for x-dos set lxdos to 3
  INTEGER,PARAMETER :: LXDOS=    1

  INTEGER           :: NMAT=     0
  INTEGER           :: NUME=     0
  INTEGER           :: NSYM=     0 
  INTEGER           :: NKPT=     0
END MODULE param


MODULE defs
  REAL*8,PARAMETER       :: CLIGHT= 137.0359895d0
  REAL*8,PARAMETER       :: PI=     3.1415926535897932d0
  REAL*8,PARAMETER       :: PI2A=     180.0d0/3.1415926535897932d0
  REAL*8,PARAMETER       :: TEST=   1.D-12
  REAL*8,PARAMETER       :: ZERO=   0.0d0
  REAL*8,PARAMETER       :: TWO=    2.0d0
  REAL*8,PARAMETER       :: NINETY=  90.0d0
  COMPLEX*16,PARAMETER   :: ZEROC=  (0.0d0,0.0d0)
  COMPLEX*16,PARAMETER   :: IMAG=   (0.0D0,1.0D0)
END MODULE defs

	module reallocate
	  !     nur 1 (generischer) Name wird von aussen angesprochen
	  interface doreallocate
	    module procedure doreallocate_r8_d1
	    module procedure doreallocate_r8_d2
	    module procedure doreallocate_c16_d1
	    module procedure doreallocate_c16_d2
	    module procedure doreallocate_i4_d1
	    module procedure doreallocate_i4_d2
	    module procedure hugo     !   ;)
	  end interface
	contains

	  !     leider sind mehrere subroutines notwendig fuer verschiedene Typen
	  subroutine doreallocate_r8_d1(tf, newdimension)
	    real*8, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
	    !     nur 1 mal kopieren reicht
	    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
	    real*8, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d1(tf, newdimension)
	    complex*16, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d2(tf, newdimension1, newdimension2)
	    complex*16, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d1(tf, newdimension)
	    integer*4, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
	    integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	!     Es gibt auch Methoden, um das Programm unleserlich zu machen :-)
	!     das sollten wir besser vermeiden ;-)
	  subroutine hugo
	    write(*,*) " Hier ist Hugo"
	  end subroutine
	end module 

MODULE struct

  INTEGER                  :: nato, ndif, nsym
  LOGICAL                  :: rel
  REAL*8                   :: AA,BB,CC,VOL,pia(3),alpha(3)
  REAL*8,ALLOCATABLE       :: R0(:),DX(:),RMT(:),zz(:),rotloc(:,:,:),v(:)
  REAL*8,ALLOCATABLE       :: tau(:,:)
  REAL*8,POINTER           :: pos(:,:)
  CHARACTER*4              :: lattic,irel,cform
  CHARACTER*80             :: title
  CHARACTER*10,ALLOCATABLE :: aname(:)
  INTEGER                  :: nat,iord
  INTEGER,ALLOCATABLE      :: mult(:),jrj(:),iatnr(:),isplit(:)
  INTEGER,ALLOCATABLE      :: iz(:,:,:),inum(:)
  REAL*8                   :: br1_rec(3,3),br2_rec(3,3)
  REAL*8                   :: br1_dir(3,3),br2_dir(3,3)
  LOGICAL                  :: ortho
  character*20             :: sgroup
  REAL*8,ALLOCATABLE         :: rotij(:,:,:),tauij(:,:) 
  integer, allocatable:: tinv(:)


 CONTAINS
  SUBROUTINE init_struct
    USE reallocate
    IMPLICIT NONE

    INTEGER                :: ios
    REAL*8                 :: test,ninety
!loop indexs
    INTEGER                :: index,i,j,j1,j2,m,jatom
    logical lmult

    test=1.D-5
    ninety=90.0D0

    read (20,1000) title

    lmult=.false.
    if (trim(title).eq.'octavetmp') lmult=.true.
    
    read (20,1010)  lattic,nat,sgroup,cform,irel
    nato=nat
    REL=.TRUE.                                     
    IF(IREL.EQ.'NREL') REL=.FALSE.                                    
    ALLOCATE(aname(nato),mult(nato),jrj(nato),r0(nato),dx(nato),rmt(nato))
    allocate(zz(nato),rotloc(3,3,nato),iatnr(nato),isplit(nato),v(nato))
    v=0.0d0
    ALLOCATE (pos(3,10000*nat))
    read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    IF(ABS(ALPHA(1)).LT.test) ALPHA(1)=ninety
    IF(ABS(ALPHA(2)).LT.test) ALPHA(2)=ninety
    IF(ABS(ALPHA(3)).LT.test) ALPHA(3)=ninety
    INDEX=0
    DO jatom=1,NAT
       INDEX=INDEX+1
       if (lmult) then
          READ(20,1040,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), &
               mult(jatom),isplit(jatom) 
       else
          READ(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), &
               mult(jatom),isplit(jatom) 
       endif
       WRITE(6,*) index,mult(jatom),ios
       IF(ios /= 0) THEN
          WRITE(6,*) iatnr(jatom),( pos(j,index),j=1,3 ), &
               mult(jatom),isplit(jatom),index 
          WRITE(6,*) 'ERROR IN STRUCT FILE READ'
          STOP
       ENDIF
       IF (mult(jatom) .EQ. 0) THEN
          WRITE (6,6000) jatom, index, mult(jatom)
          STOP
       ENDIF
       DO m=1,mult(jatom)-1                                     
          index=index+1                                            
          READ(20,1031) iatnr(jatom),( pos(j,index),j=1,3)         
       ENDDO
       READ(20,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom), &
            zz(jatom)
       dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)           
       rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jrj(jatom)-1) )           
       READ(20,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
    ENDDO
    ndif=index
    CALL doreallocate(pos, 3, ndif)
    READ(20,1151) iord
    nsym=iord
    ALLOCATE(iz(3,3,nsym+1),tau(3,nsym+1),inum(nsym+1))
    allocate(rotij(1:3,1:3,ndif),tauij(1:3,ndif))
    DO j=1,iord
       READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
    ENDDO
    allocate(tinv(nsym))

1000 FORMAT(A80)                                                       
1010 FORMAT(A4,22X,I4,a20,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(3X,I5,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1040 FORMAT(3X,I5,4X,F10.7,3X,F10.7,3X,F10.7,/,11X,I6,17X,I2)  
1031 FORMAT(3X,I5,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)                                             
1101 FORMAT(3(3I2,F11.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE init_struct

  SUBROUTINE latgen_struct
    USE defs
    IMPLICIT NONE
    INTEGER    :: i,j
    REAL*8     :: sqrt3,rvfac,sinab,sinbc,cosab,cosac,cosbc,wurzel,x1,y1,z1
    REAL*8     :: det

    SQRT3=SQRT(3.D0)
    ALPHA(1)=ALPHA(1)*PI/180.0D0                                             
    ALPHA(2)=ALPHA(2)*PI/180.0D0                                             
    ALPHA(3)=ALPHA(3)*PI/180.0D0                                             
    PIA(1)=2.D0*PI/AA                                                 
    PIA(2)=2.D0*PI/BB                                                 
    PIA(3)=2.D0*PI/CC
    cosab=COS(alpha(3))
    cosac=COS(alpha(2))
    cosbc=COS(alpha(1))
    sinab=SIN(alpha(3))
    sinbc=SIN(alpha(1))

    br1_rec=zero; br2_rec=zero

    SELECT CASE (LATTIC(1:1))
    CASE ('H')       !.....HEXAGONAL LATTICE
       br1_rec(1,1)=2.d0/sqrt3*pia(1)                                        
       br1_rec(1,2)=1.d0/sqrt3*pia(1)                                        
       br1_rec(2,2)=pia(2)                                                   
       br1_rec(3,3)=pia(3)                                                   

       br2_rec(1:3,1:3)=br1_rec(1:3,1:3)
       rvfac=2.d0/SQRT(3.d0)                                             
       ortho=.FALSE.
    CASE ('S','P' )
       wurzel=SQRT(sinbc**2-cosac**2-cosab**2+2*cosbc*cosac*cosab)
       br1_rec(1,1)= sinbc/wurzel*pia(1)
       br1_rec(1,2)= (-cosab+cosbc*cosac)/(sinbc*wurzel)*pia(2)
       br1_rec(1,3)= (cosbc*cosab-cosac)/(sinbc*wurzel)*pia(3)
       br1_rec(2,2)= pia(2)/sinbc
       br1_rec(2,3)= -pia(3)*cosbc/sinbc
       br1_rec(3,3)= pia(3)
       !
       br2_rec(1:3,1:3)=br1_rec(1:3,1:3)
       rvfac= 1.d0/wurzel
       ortho=.TRUE.
       IF(ABS(alpha(1)-pi/2.d0).GT.test) ortho=.FALSE.
       IF(ABS(alpha(2)-pi/2.d0).GT.test) ortho=.FALSE.
       IF(ABS(alpha(3)-pi/2.d0).GT.test) ortho=.FALSE.
    CASE ( 'F')
       BR1_REC(1,1)=PIA(1)                                                   
       BR1_REC(2,2)=PIA(2)                                                   
       BR1_REC(3,3)=PIA(3)                                                   

!     definitions according to column, rows convention for BR2_REC
       BR2_REC(1,1)=-PIA(1)                                                  
       BR2_REC(1,2)= PIA(1)                                                  
       BR2_REC(1,3)= PIA(1)                                                  
       BR2_REC(2,1)= PIA(2)                                                  
       BR2_REC(2,2)=-PIA(2)                                                  
       BR2_REC(2,3)= PIA(2)                                                  
       BR2_REC(3,1)= PIA(3)                                                  
       BR2_REC(3,2)= PIA(3)                                                  
       BR2_REC(3,3)=-PIA(3)                                                  
       !                                                                       
       RVFAC=4.D0                                                        
       ORTHO=.TRUE.                                             
    CASE ( 'B' )
       BR1_REC(1,1)=PIA(1)                                                   
       BR1_REC(2,2)=PIA(2)                                                   
       BR1_REC(3,3)=PIA(3)                                                   
       !                                                                       
       BR2_REC(1,1)= ZERO                                                    
       BR2_REC(1,2)= PIA(1)                                                  
       BR2_REC(1,3)= PIA(1)                                                  
       BR2_REC(2,1)= PIA(2)                                                  
       BR2_REC(2,2)= ZERO                                                    
       BR2_REC(2,3)= PIA(2)                                                  
       BR2_REC(3,1)= PIA(3)                                                  
       BR2_REC(3,2)= PIA(3)                                                  
       BR2_REC(3,3)= ZERO
       !
       RVFAC=2.D0
       ORTHO=.TRUE.                                             
    CASE( 'C' )
       SELECT CASE (lattic(2:3))
       CASE ( 'XY' )
          BR1_REC(1,1)=PIA(1)
          BR1_REC(2,2)=PIA(2)
          BR1_REC(3,3)=PIA(3)

          BR2_REC(1,1)= PIA(1)
          BR2_REC(1,2)= PIA(1)
          BR2_REC(1,3)= ZERO                                                  
          BR2_REC(2,1)=-PIA(2)
          BR2_REC(2,2)= PIA(2)
          BR2_REC(2,3)= ZERO                                                 
          BR2_REC(3,1)= ZERO                                                  
          BR2_REC(3,2)= ZERO                                                 
          BR2_REC(3,3)= PIA(3)

          RVFAC=2.D0                                                        
          ORTHO=.TRUE.                                             
       CASE( 'XZ ' )
          IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
             BR1_REC(1,1)=PIA(1)
             BR1_REC(2,2)=PIA(2)
             BR1_REC(3,3)=PIA(3)

             BR2_REC(1,1)= PIA(1)
             BR2_REC(1,2)= zero                                                   
             BR2_REC(1,3)= PIA(1)
             BR2_REC(2,1)= zero
             BR2_REC(2,2)= PIA(2)
             BR2_REC(2,3)= zero 
             BR2_REC(3,1)=-PIA(3)
             BR2_REC(3,2)= zero   
             BR2_REC(3,3)= PIA(3)
             RVFAC=2.0                                                         
             ORTHO=.TRUE.                                             
          ELSE
             write(6,*) '  gamma not equal 90' !CXZ MONOCLINIC CASE 
             BR1_REC(1,1)= PIA(1)/SINAB 
             BR1_REC(1,2)= -PIA(2)*COSAB/SINAB
             BR1_REC(2,2)= PIA(2)
             BR1_REC(3,3)= PIA(3)

             BR2_REC(1,1)= PIA(1)/SINAB 
             BR2_REC(1,2)= -PIA(2)*COSAB/SINAB
             BR2_REC(1,3)= PIA(1)/SINAB 
             BR2_REC(2,1)= zero 
             BR2_REC(2,2)= PIA(2)
             BR2_REC(2,3)= zero
             BR2_REC(3,1)=-PIA(3)
             BR2_REC(3,2)= zero
             BR2_REC(3,3)= PIA(3)

             RVFAC=2.0/SINAB                                                   
             ORTHO=.FALSE.                                             
          ENDIF
       CASE( 'YZ' )
          BR1_REC(1,1)=PIA(1)
          BR1_REC(2,2)=PIA(2)
          BR1_REC(3,3)=PIA(3)
          BR2_REC(1,1)= PIA(1)
          BR2_REC(1,2)= zero
          BR2_REC(1,3)= zero  
          BR2_REC(2,1)= zero  
          BR2_REC(2,2)= PIA(2)
          BR2_REC(2,3)= PIA(2)
          BR2_REC(3,1)= zero  
          BR2_REC(3,2)=-PIA(3)
          BR2_REC(3,3)= PIA(3)
          RVFAC=2.0d0                                                         
          ORTHO=.TRUE.
       END SELECT
    CASE ( 'R' )
       BR1_REC(1,1)=1.D0/SQRT(3.D0)*PIA(1)
       BR1_REC(1,2)=1.D0/SQRT(3.D0)*PIA(1)
       BR1_REC(1,3)=-2.d0/sqrt(3.d0)*PIA(1)
       BR1_REC(2,1)=-1.0d0*PIA(2)
       BR1_REC(2,2)=1.0d0*PIA(2)
       BR1_REC(2,3)=zero
       BR1_REC(3,1)=1.0d0*PIA(3)
       BR1_REC(3,2)=1.0d0*PIA(3)
       BR1_REC(3,3)=1.0d0*PIA(3)
       
       br2_rec(1:3,1:3)=br1_rec(1:3,1:3)
       
       RVFAC=6.D0/SQRT(3.D0)
       ORTHO=.FALSE.
    CASE DEFAULT
       STOP 'LATGEN - Wrong lattice'
    END SELECT
    vol=aa*bb*cc/rvfac
    CALL invert_struct(br1_rec,br1_dir)
    det=0.d0                                                          
    DO i=1,3                                                      
       det=det+br1_dir(i,1)*br1_rec(i,1)                                       
    ENDDO
    br1_dir(1:3,1:3)=br1_dir(1:3,1:3)*2.d0*PI/det
    CALL invert_struct(br2_rec,br2_dir)
    det=0.d0                                                          
    DO i=1,3                                                      
       det=det+br2_dir(i,1)*br2_rec(i,1)                                       
    ENDDO
    br2_dir(1:3,1:3)=br2_dir(1:3,1:3)*2.d0*PI/det
    write(6,'(a)')   '  convenction (i,j)     (j)'
    write(6,'(a)')   '                      ax bx cx '
    write(6,'(a)')   '                  (i) ay by cy '
    write(6,'(a,/)') '                      az bz cz '
    WRITE(6,*)  '------------BR1_REC-----------'
    WRITE(6,106)((BR1_rec(I,J),I=1,3),J=1,3)                              
    WRITE(6,*)  '------------BR2_REC-----------'
!br_dir are transposed br_rec**(-1)*2*pi
    WRITE(6,106)((BR2_rec(I,J),I=1,3),J=1,3)                              
    WRITE(6,*)  '------------BR1_DIR-----------'
    WRITE(6,106)((br1_dir(I,J),I=1,3),J=1,3)                              
    WRITE(6,*)  '------------BR2_DIR-----------'
    WRITE(6,106)((br2_dir(I,J),I=1,3),J=1,3)                              
106 FORMAT(3(10X,3F10.5,/))                      
  END SUBROUTINE latgen_struct

  SUBROUTINE invert_struct(mat1,inv_mat1)
    REAL*8        :: mat1(3,3),inv_mat1(3,3)
    inv_mat1(1,1)=mat1(2,2)*mat1(3,3)-mat1(3,2)*mat1(2,3)                 
    inv_mat1(2,1)=mat1(3,2)*mat1(1,3)-mat1(1,2)*mat1(3,3)                 
    inv_mat1(3,1)=mat1(1,2)*mat1(2,3)-mat1(2,2)*mat1(1,3)                 
    inv_mat1(1,2)=mat1(2,3)*mat1(3,1)-mat1(3,3)*mat1(2,1)                 
    inv_mat1(2,2)=mat1(3,3)*mat1(1,1)-mat1(1,3)*mat1(3,1)                 
    inv_mat1(3,2)=mat1(1,3)*mat1(2,1)-mat1(2,3)*mat1(1,1)                 
    inv_mat1(1,3)=mat1(2,1)*mat1(3,2)-mat1(3,1)*mat1(2,2)                 
    inv_mat1(2,3)=mat1(3,1)*mat1(1,2)-mat1(1,1)*mat1(3,2)                 
    inv_mat1(3,3)=mat1(1,1)*mat1(2,2)-mat1(2,1)*mat1(1,2)                 
  END SUBROUTINE invert_struct

    subroutine write_struct(itap)
  
      implicit none
  
      integer itap
      integer                :: index,i,j,j1,j2,m,jatom
      
      write (itap,1000) title
      write (itap,1010) lattic,nat,sgroup
      write (itap,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
      index=0
      do jatom=1,nat
         index=index+1
         write(itap,1030) iatnr(jatom),(pos(j,index),j=1,3 ), &
              mult(jatom),isplit(jatom) 
         do m=1,mult(jatom)-1                                     
            index=index+1                                            
            write(itap,1031) iatnr(jatom),(pos(j,index),j=1,3)         
         enddo
         write(itap,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom),zz(jatom)
         write(itap,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
      enddo
     
      write(itap,1151) iord
      DO j=1,iord
         write(itap,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),j,tinv(j)
      ENDDO
      
1000  format(a80)                                                       
1010  format(a4,22x,i4,a20/,&
           'MODE OF CALC=RELA unit=bohr')                                 
1020  format(6f10.6)                                          
1030  FORMAT('AT ',I5,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,&
           '          M',i6,'          ISPLIT=',i2)          
1031  FORMAT(3X,I5,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)                   
1050  format(a10,' NPT=',i5,'  R0=',f10.9,' RMT=',f10.5,'   Z:',f10.5)
1051  format('LOCAL ROT MATRIX:   ',3f10.7,/,20x,3f10.7,/,20x,3f10.7)
1151  format(i4,'      NUMBER OF SYMMETRY OPERATIONS')
1101  FORMAT(3(3I2,F11.8/),2I8)
      
    end subroutine write_struct

END MODULE struct


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 module rotations

   use struct
   use defs

  real*8, allocatable     :: fi(:),teta(:)
  integer, allocatable    :: ncmopt(:)
  integer                 :: nmag
  logical, allocatable    :: matom(:)
  integer, allocatable    :: tmatom(:)
  integer ideb,lmax
  character ncmmode*30
  real*8 cfmix,ncmq(3)
  logical nmagmod,sspiral,fixnmag

  contains
    
    subroutine init_rotations

      integer i,ifixnmag

      allocate (fi(ndif),teta(ndif),matom(ndif),ncmopt(ndif),tmatom(ndif))
      
      
      fi=0.0d0
      teta=0.0d0
      ncmopt=0
      nmag=ndif
      ncmq=0.0d0
      do i=1,ndif
         tmatom(i)=i
      enddo
      lmax=6
      ifixnmag=0
      rewind(4)
      read(4,*,end=4,err=4) ncmmode
      read(4,*,end=4,err=4) (ncmq(i),i=1,3)
      do i=1,ndif
         read(4,*) fi(i),teta(i),ncmopt(i)
      enddo
      read(4,*) cfmix
      read(4,*,err=1,end=1) nmag,(tmatom(i),i=1,nmag)
      read(4,*,err=2,end=2) lmax
      read(4,*,err=3,end=3) ifixnmag
      goto 3
1     nmag=0
2     lmax=6      
3     fixnmag=.true.      
      ifixnmag=1
4     continue

      fixnmag=.true.
      if (ifixnmag.eq.0) fixnmag=.false.

      write(0,*) 'lmax=6'
      write(0,*) 'fixnmag',fixnmag

      fi=fi*pi/180.0d0
      teta=teta*pi/180.0d0

      qq=sqrt(ncmq(1)**2+ncmq(2)**2+ncmq(3)**2)
      sspiral=.false.
      if (qq.gt.1.0d-5) sspiral=.true.

! true for magnetic
      matom=.true.
      do i=1,nmag
         matom(tmatom(i))=.false.         
      enddo      
      
      nmagmod=.true.
      do i=1,ndif
         nmagmod=nmagmod.and.(.not.matom(i))
      enddo

    end subroutine init_rotations

    subroutine write_inncm
      
      integer nmag1

      write(3,'(a)') ncmmode
      write(3,'(3f12.6)') ncmq
      do i=1,ndif
         write(3,'(2f12.6,i5)') mod(fi(i)/(pi/180.0d0),360.0d0),&
              mod(teta(i)/(pi/180.0d0),360.0d0),ncmopt(i)
      enddo
      write(3,'(f12.6)') cfmix
      nmag1=0
      write(3,'(30i5)') nmag1,(tmatom(i),i=1,nmag)
      write(3,'(i5)') lmax

    end subroutine write_inncm

  end module rotations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module Ylm_rot
    integer    :: npgop
    real*8, allocatable     :: pg_sym_oper(:,:,:,:),pgtinv(:,:)
    complex*16, allocatable :: Ylm_rot_mat(:,:,:,:,:)
    complex*16, allocatable :: Ylm_rot_matti(:,:,:,:,:)
    complex*16, allocatable :: Ylm_rot_mat_dmat(:,:,:,:,:)
    complex*16, allocatable :: spin_rot_mat(:,:,:,:,:)
    integer, allocatable    :: npgopat(:)

  contains
    
    subroutine init_Ylm_rot

      use struct, only:nat,nsym
      use rotations, only: lmax

      npgop=nsym


      allocate (pg_sym_oper(nat,npgop,3,3),pgtinv(nat,npgop),npgopat(nat))
      allocate (Ylm_rot_mat(nat,npgop,0:lmax,-lmax:lmax,-lmax:lmax))       
      allocate (Ylm_rot_matti(nat,npgop,0:lmax,-lmax:lmax,-lmax:lmax))       
      allocate (Ylm_rot_mat_dmat(nat,npgop,0:lmax,-lmax:lmax,-lmax:lmax)) 
      allocate (spin_rot_mat(nat,npgop,3,2,2)) 

      Ylm_rot_mat=(0.0d0,0.0d0)
      Ylm_rot_matti=(0.0d0,0.0d0)
      Ylm_rot_mat_dmat=(0.0d0,0.0d0)      
      spin_rot_mat=(0.0d0,0.0d0) 

    end subroutine init_Ylm_rot

  end module Ylm_rot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
