MODULE defs
  REAL*8,PARAMETER       :: CLIGHT= 137.0359895d0
  REAL*8,PARAMETER       :: PI=     3.1415926535897932d0
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

 module struct

  CHARACTER*80             :: title
  CHARACTER*4              :: lattic,irel
  REAL*8                   :: aa,bb,cc,pia(3),alpha(3)
  REAL*8                   :: vol

  INTEGER                  :: nat,ndif
  INTEGER,ALLOCATABLE      :: mult(:),jrj(:),iatnr(:),isplit(:)
  REAL*8,ALLOCATABLE       :: r0(:),dx(:),rmt(:),zz(:)
  CHARACTER*10,ALLOCATABLE :: aname(:)
  REAL*8,POINTER           :: pos(:,:)

  INTEGER                  :: nsym,iord
  INTEGER,POINTER          :: iz(:,:,:),inum(:)
  REAL*8,POINTER           :: tau(:,:)
  LOGICAL                  :: ortho
  REAL*8, allocatable                   :: rotloc(:,:,:)
  REAL*8                   :: br1_rec(3,3),br2_rec(3,3)
  REAL*8                   :: br1_dir(3,3),br2_dir(3,3)

 contains

 subroutine init_struct

    USE reallocate
    USE defs
    IMPLICIT NONE

    integer nato
    INTEGER                :: ios
    REAL*8                 :: test0

    INTEGER                :: index,i,j,j1,j2,m,jatom
    logical lmult

    test0=1.D-5

    rewind(20)

    read (20,1000) title
    lmult=.false.
    if (trim(title).eq.'octavetmp') lmult=.true.
    read (20,1010) lattic,nat,irel
    ALLOCATE(aname(nat),mult(0:nat),jrj(nat),r0(nat),dx(nat),rmt(nat),zz(nat),rotloc(3,3,nat),iatnr(nat),isplit(nat))
    ALLOCATE (pos(3,1000*nat))
    mult(0)=0
    read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    IF(ABS(ALPHA(1)).LT.test0) ALPHA(1)=ninety
    IF(ABS(ALPHA(2)).LT.test0) ALPHA(2)=ninety
    IF(ABS(ALPHA(3)).LT.test0) ALPHA(3)=ninety
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
       IF(ios /= 0) THEN
          WRITE(6,*) iatnr(jatom),( pos(j,index),j=1,3 ), &
               mult(jatom),isplit(jatom) 
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
    ALLOCATE(iz(3,3,nsym),tau(3,nsym),inum(nsym))
    DO j=1,iord
       READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
    ENDDO
  
1000 FORMAT(A80)                                                       
1010 FORMAT(A4,22X,I4,/,13X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,11X,I6,17X,I2)          
1031 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)                                             
1101 FORMAT(3(3I2,F11.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN READING STRUCT : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I4,3X,'INDEX=',I3,3X,'MULT=',I3)

 end subroutine init_struct

 subroutine finit_struct

   deallocate(aname,mult,jrj,r0,dx,rmt,zz,rotloc,iatnr,isplit,pos)

 end subroutine finit_struct
   
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
    CALL inversa(br1_rec,br1_dir)
    br1_dir(1:3,1:3)=br1_dir(1:3,1:3)*2.d0*PI
    CALL inversa(br2_rec,br2_dir)
    DO i=1,3                                                      
       DO j=1,3                                                      
          br2_dir(i,j)=br2_dir(i,j)*2.d0*PI
       ENDDO
    ENDDO

    ALPHA(1)=ALPHA(1)/(PI/180.0D0)                                            
    ALPHA(2)=ALPHA(2)/(PI/180.0D0)                                             
    ALPHA(3)=ALPHA(3)/(PI/180.0D0)

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
    inv_mat1=transpose(inv_mat1)
  END SUBROUTINE invert_struct

   subroutine inversa(a,ainv)
     implicit real*8 (a-h,o-z)
     dimension a(3,3),ainv(3,3)
     det= a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
          +a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3) &
          -a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
     ainv(1,1) =(   a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
     ainv(2,1) =( - a(2,1) * a(3,3) + a(2,3) * a(3,1) ) / det
     ainv(3,1) =(   a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
     ainv(1,2) =( - a(1,2) * a(3,3) + a(1,3) * a(3,2) ) / det
     ainv(2,2) =(   a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
     ainv(3,2) =( - a(1,1) * a(3,2) + a(1,2) * a(3,1) ) / det
     ainv(1,3) =(   a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det
     ainv(2,3) =( - a(1,1) * a(2,3) + a(1,3) * a(2,1) ) / det
     ainv(3,3) =(   a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det
     return
   end subroutine inversa

 end module struct

 module nstruct

   character*80             :: n_title
   character*4              :: n_lattic,n_irel
   real*8                   :: n_aa,n_bb,n_cc,n_alpha(3)
   integer,allocatable      :: n_mult(:),n_jrj(:),n_iatnr(:),n_isplit(:)
   real*8,allocatable       :: n_r0(:),n_dx(:),n_rmt(:),n_zz(:)
   character*10,allocatable :: n_aname(:)
   real*8,pointer           :: n_pos(:,:)

   integer                  :: n_nat,n_ndif    
   integer                  :: n_nsym,n_iord
   integer,pointer          :: n_iz(:,:,:),n_inum(:)
   real*8,pointer           :: n_tau(:,:)
   real*8, allocatable      :: n_rotloc(:,:,:)
   real*8   hex2rho(3,3),rho2hex(3,3)
 
 contains

   subroutine init_nstruct(nat)

     ALLOCATE(n_aname(nat),n_mult(0:nat),n_jrj(nat),&
          n_r0(nat),n_dx(nat),n_rmt(nat),n_zz(nat),n_rotloc(3,3,nat),&
          n_iatnr(nat),n_isplit(nat),n_pos(3,nat))
     ALLOCATE(n_iz(3,3,1),n_tau(3,1),n_inum(1))

   end subroutine init_nstruct

   subroutine write_nstruct
  
     implicit none
     
     integer                :: index,i,j,j1,j2,m,jatom
     integer sgnum
     character sgname*8
     logical lmult

     sgnum=1
     sgname='P1'

     lmult=.false.
     if (trim(n_title).eq.'octavetmp') lmult=.true.

     write (20,1000) n_title
     write (20,1010) n_lattic,n_nat,sgnum,sgname
     write (20,1020) n_aa,n_bb,n_cc,n_alpha(1),n_alpha(2),n_alpha(3)
     index=0
     do jatom=1,n_nat
        index=index+1
        if (lmult) then
           write(20,1040) n_iatnr(jatom),( n_pos(j,index),j=1,3 ), &
                n_mult(jatom),n_isplit(jatom) 
        else
           write(20,1030) n_iatnr(jatom),( n_pos(j,index),j=1,3 ), &
                n_mult(jatom),n_isplit(jatom) 
        endif
        do m=1,n_mult(jatom)-1                                     
           index=index+1                                            
           write(20,1031) n_iatnr(jatom),(n_pos(j,index),j=1,3)         
        enddo
        write(20,1050) n_aname(jatom),n_jrj(jatom),n_r0(jatom),n_rmt(jatom), &
             n_zz(jatom)
        write(20,1051) ((n_rotloc(i,j,jatom),i=1,3),j=1,3)                
     enddo
     write(20,1151) n_iord
     DO j=1,n_iord
        write(20,1101) ( (n_iz(j1,j2,j),j1=1,3),n_tau(j2,j),j2=1,3 ),n_inum(j)
     ENDDO
     
1000 format(a80)                                                       
1010 format(a4,22x,i4,2x,i3,a8/,&
          'MODE OF CALC=RELA unit=bohr')                                 
1020 format(6f10.6)                                          
1030 FORMAT('AT ',I5,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,&
          '          MULT=',i2,'          ISPLIT=',i2)          
1040 FORMAT('AT ',I5,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,&
          '          M',i6,'          ISPLIT=',i2)          
1031 FORMAT(3X,I5,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)                   
1050 format(a10,' NPT=',i5,'  R0=',f10.9,' RMT=',f10.5,'   Z:',f10.5)
1051 format('LOCAL ROT MATRIX:   ',3f10.7,/,20x,3f10.7,/,20x,3f10.7)
1151 format(i4,'      NUMBER OF SYMMETRY OPERATIONS')
1101 FORMAT(3(3I2,F11.8/),I8)
     
   end subroutine write_nstruct

   subroutine def_rho2hex

     use struct, only: inversa

     real*8      hex2ort(3,3),ort2rho(3,3)

     hex2ort(1,1)=0.0d0
     hex2ort(1,2)=1.0d0
     hex2ort(1,3)=0.0d0
     hex2ort(2,1)=sqrt(3.0d0)/2.0d0
     hex2ort(2,2)=-0.5d0
     hex2ort(2,3)=0.0d0
     hex2ort(3,1)=0.0d0
     hex2ort(3,2)=0.0d0
     hex2ort(3,3)=1.0d0
     ort2rho(1,1)=1.0d0/sqrt(3.0d0)
     ort2rho(1,2)=1.0d0/sqrt(3.0d0)
     ort2rho(1,3)=-2.0d0/sqrt(3.0d0)
     ort2rho(2,1)=-1.0d0
     ort2rho(2,2)=1.0d0
     ort2rho(2,3)=0.0d0
     ort2rho(3,1)=1.0d0
     ort2rho(3,2)=1.0d0
     ort2rho(3,3)=1.0d0
     hex2rho=matmul(hex2ort,ort2rho)
     call inversa(hex2rho,rho2hex)

   end subroutine def_rho2hex

 end module nstruct

