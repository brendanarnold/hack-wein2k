MODULE struct

  real*8       b2a,a2b
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


 CONTAINS
  SUBROUTINE init_struct
    USE reallocate
    IMPLICIT NONE

    INTEGER                :: ios
    REAL*8                 :: test,ninety
!loop indexs
    INTEGER                :: index,i,j,j1,j2,m,jatom

    test=1.D-5
    ninety=90.0D0
    b2a=1.0d0!0.521977d0
    a2b=1.0d0/0.521977d0
    read (20,1000) title
    read (20,1010)  lattic,nat,cform,irel
    nato=nat
    REL=.TRUE.                                     
    IF(IREL.EQ.'NREL') REL=.FALSE.                                    
    ALLOCATE(aname(nato),mult(nato),jrj(nato),r0(nato),dx(nato),rmt(nato))
    allocate(zz(nato),rotloc(3,3,nato),iatnr(nato),isplit(nato),v(nato))
    v=0.0d0
    ALLOCATE (pos(3,48*nat))
    read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    IF(ABS(ALPHA(1)).LT.test) ALPHA(1)=ninety
    IF(ABS(ALPHA(2)).LT.test) ALPHA(2)=ninety
    IF(ABS(ALPHA(3)).LT.test) ALPHA(3)=ninety
    INDEX=0
    DO jatom=1,NAT
       INDEX=INDEX+1
       READ(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), &
            mult(jatom),isplit(jatom) 
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
1010 FORMAT(A4,22X,I4,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(3X,I5,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(3X,I5,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)                                             
1101 FORMAT(3(3I2,F11.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE init_struct

END MODULE struct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atom_list

  use struct, only : nat

  real*8, allocatable  :: list_pos(:,:)
  real*8, allocatable  :: list_rpos(:,:)
  integer, allocatable :: list_iat(:)
  integer, allocatable :: list_ind(:)
  real*8               :: br(3,3)
  integer :: nlist,natlist

  contains

  subroutine init_atom_list

    natlist=300*nat

    allocate (list_pos(3,natlist))
    allocate (list_rpos(3,natlist))
    allocate (list_iat(natlist))
    allocate (list_ind(natlist))

    nlist=0
  end subroutine init_atom_list

end module atom_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module povprop

  use struct, only : nat,ndif

  real*8, allocatable  :: rgb(:,:),rgb_tmp(:,:)
  character*5,allocatable :: rgb_name(:)
  integer              :: ncol 
  real*8, allocatable  :: rgba(:,:)
  integer              :: nvdwr
  real*8, allocatable  :: vdwr(:),vdwr_tmp(:)
  character*5,allocatable :: vdwr_name(:)
  real*8, allocatable  :: cyll(:),cylw(:)
  real*8, allocatable  :: coneh(:),conew(:)
  real*8, allocatable  :: boxrgb(:)
  real*8               :: boxw
  real*8, allocatable  :: cargb(:)
  real*8               :: caw,cal
  integer              :: iboxfin,icafin,iatfin,imagfin
  integer              :: ibondcol,ibondfin
  real*8, allocatable  :: rgbbond(:)
  integer, allocatable :: bond(:,:)
  real*8,  allocatable :: bondd(:)
  character*5,allocatable :: bond_name(:,:)   
  real*8, allocatable  :: bond_dis(:)
  integer              :: nbond,maxnb
  real*8               :: lgradbond,bondw
  integer, allocatable :: supercell(:)
  real*8, allocatable  :: qspir(:)
  real*8        :: minmm,maxmm,mincm

contains
  
  subroutine init_povprop
    
    maxnb=200+30*ndif**2
    allocate (rgb(3,nat),rgb_tmp(3,nat),rgb_name(nat))
    allocate (vdwr(nat),vdwr_tmp(nat),vdwr_name(nat))
    allocate (rgba(3,nat))
    allocate (cyll(nat),cylw(nat))
    allocate (coneh(nat),conew(nat))    
    allocate (boxrgb(3),cargb(3),rgbbond(3))    
    allocate (bond(maxnb,2),bondd(maxnb))  
    allocate (bond_name(maxnb,2), bond_dis(maxnb))
    allocate (supercell(3))
    allocate (qspir(3)) 

    bond_name(:,:)='     '
    do i=1,nat
       rgb(1,i)=1.0d0
       rgb(2,i)=0.0d0
       rgb(3,i)=0.0d0
       rgb_name(i)='    ' 
       vdwr_name(i)='    ' 
       vdwr(i)=0.5d0
       cyll(i)=1.0d0
       cylw(i)=0.15d0
       coneh(i)=0.4d0
       conew(i)=0.25d0
    end do
    
    boxrgb(1:3)=0.0d0
    boxw=0.05
    cargb(1:3)=0.0d0
    caxw=0.05
    cal=10.0d0

  end subroutine init_povprop
  
end module povprop
