MODULE param
      INTEGER            LMAX2, LOMAX
      INTEGER            NRAD, nloat, nrf
      integer            iblck
      real*8             clight
      parameter (IBLCK=   64)
      PARAMETER (IBLOCK= 128)
!.....Optimize IBLOCK and IBLCK for your hardware (32-255)
      PARAMETER (LMAX2=   3)                                              
      PARAMETER (LOMAX=   3)
      PARAMETER (NDIM=  2*lmax2+1)
      PARAMETER  (NDIM2= 2*NDIM)                                              
      INTEGER          :: NATO=   0
! for ncom parameter check format 1003 in l2main.frc
      INTEGER          :: NDIF=   0
      INTEGER          :: NKPT= 0                                              
      INTEGER          :: NMAT=0
      PARAMETER (NRAD=  881)                                              
      INTEGER          :: NSYM=   0 
      INTEGER          :: NUME=   0
! for x-dos set lxdos to 3
      parameter (lxdos= 3)
      parameter (nloat= 3)
      parameter (nrf=4)
      PARAMETER (CLIGHT=137.0359895D0)
END MODULE param

MODULE struct
  USE param
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
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)                                             
1101 FORMAT(3(3I2,F10.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE init_struct

END MODULE struct

MODULE case
  COMPLEX*16,ALLOCATABLE     :: cf(:,:,:)
  INTEGER,ALLOCATABLE    :: isum(:,:),iatom(:),lcase(:)
  INTEGER                :: natom

 CONTAINS
  SUBROUTINE init_coef(ndim2)
    READ(5,*)NATOM
    ALLOCATE(CF(NDIM2,NDIM2,natom),ISUM(0:ndim2,natom),iatom(natom),lcase(natom))
    DO I=1,NATOM
    READ(5,*)IATOM(I),LCASE(I)
    write(6,*)'atom number: ',IATOM(I),'   orbit. num. l=',LCASE(I)
    END DO
  END SUBROUTINE init_coef
END MODULE case

MODULE sym2
   INTEGER,ALLOCATABLE    :: idet(:)
   REAL*8,ALLOCATABLE     :: opimat(:,:,:),phase(:),tmat(:,:,:)
 CONTAINS
  SUBROUTINE init_sym2(iord)
   ALLOCATE(idet(iord),opimat(3,3,iord),phase(iord),tmat(3,3,iord))
   idet=1
  END SUBROUTINE init_sym2
END MODULE sym2

MODULE com
   INTEGER,ALLOCATABLE    :: nb(:)
   INTEGER                :: nband,iso
   REAL*8                 :: emin,emax,ef
 CONTAINS
  SUBROUTINE init_com(nkpt)
   ALLOCATE(nb(nkpt))
  END SUBROUTINE init_com
END MODULE com

MODULE abc
   INTEGER,ALLOCATABLE    ::  kx(:),ky(:),kz(:)
!      ,nume,nmat
   REAL*8,ALLOCATABLE     ::  bkx(:),bky(:),bkz(:)
   REAL*8,ALLOCATABLE     ::  e(:),xden(:,:)
   COMPLEX*16,ALLOCATABLE ::  a(:,:,:),alm(:,:,:,:),xqtl(:,:,:)
 CONTAINS
  SUBROUTINE init_abc(nume,nmat,lmax2,ndim2,nrf)
   ALLOCATE(kx(nmat),ky(nmat),kz(nmat),bkx(nmat),bky(nmat),bkz(nmat))
   ALLOCATE(a(nmat,nume,2),alm(2*lmax2+1,nume,nrf,2))
   ALLOCATE(e(nume),xden(ndim2,nume),xqtl(ndim2,ndim2,nume))
   END SUBROUTINE init_abc
END MODULE abc
  

