      SUBROUTINE INIT(FL,NTYPE,NLO1,IZ,IIZ,TAU,IORD,BR1,BR2,DR1,DB1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
! 
      LOGICAL          FL(FLMAX),FL1,FL2,ORTHO
      CHARACTER*4      LATT
      CHARACTER*80     TITLE
      CHARACTER*11     ANA
!      DIMENSION        RO(NATO),DXM(NATO),RMT(NATO), &
!                       ANA(NATO),JRJ(NATO)
      integer, allocatable ::     NATOM(:)    ! nato
      DIMENSION        IZ(3,3,NSYM),TAU(3,NSYM),IIZ(3,3,NSYM)
      DIMENSION        ALPHA(3),POS(3)
      DIMENSION        BR1(3,3),BR2(3,3)
      DIMENSION        DR1(3,3),DR2(3,3)
      DIMENSION        DB1(3,3)
      integer, allocatable :: NLO(:),NLOV(:),NLON(:) ! nato
      DIMENSION        ELO(0:LOMAX,nloat),emist(0:LMAX2)
!******************************************************************************
!
!.....write date
      CALL WRTDATE(5)
      CALL WRTDATE(6)

!.....input from .struct
      READ(20,530) TITLE
      READ(20,541) LATT,NTYPE
      READ(20,530) 
      READ(20,580) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      allocate (natom(ntype),NLO(Ntype),NLOV(Ntype),NLON(Ntype))
      INDEX=0
      DO 10 ITY=1,NTYPE
        INDEX=INDEX+1
        READ(20,543) (POS(I),I=1,3)
        READ(20,544) NATOM(ITY)
        DO 12 INA=2,NATOM(ITY)
          INDEX=INDEX+1
!          IF(INDEX.GT.NDIF) STOP 'ERROR in init: Too many atoms'
          READ(20,543) (POS(I),I=1,3)
 12    CONTINUE 
        READ(20,545) ANA,JRJ,RO,RMT
!        DXM(ITY)=DLOG(RMT(ITY)/RO(ITY))/(JRJ(ITY)-1)
        READ(20,*) 
        READ(20,*) 
        READ(20,*) 
 10   CONTINUE
      FL1=.FALSE.
      FL2=.FALSE.
      READ(20,*) IORD 
      DO 14 I=1,IORD                                                    
        READ(20,512) ((IZ(I1,I2,I),I2=1,3),TAU(I1,I),I1=1,3)  
        READ(20,514) NOIORD
        CALL INVMATI(IZ(1,1,I),IIZ(1,1,I))
!
        STAU=0.D0
        ISIZ=0
        DO 20 I1=1,3
          STAU=STAU+DABS(TAU(I1,I))
          ISIZ=ISIZ+ABS(IZ(I1,I1,I)+1)
          DO 20 I2=1,3
            IF(I2.NE.I1) ISIZ=ISIZ+ABS(IZ(I1,I2,I))
 20     CONTINUE    
        IF(STAU.GT.1E-8) FL1=.TRUE.
        IF(ISIZ.EQ.0)    FL2=.TRUE.
 14   CONTINUE
      IF(FL2) FL(3)=.TRUE.
!
!.....tranformation matrices
!     BR1(3,3)      -  from lapw1 k-list coord. into Cartesian coord (recipr. space)
!     BR2(3,3)      -  from reciprocal coord.   into Cartesian coord (recipr. space)
!     DR1,DR2       -  inv(BR1), inv(BR2)
!     DB1(3,3)      -  DR2*BR1; from lapw1 k-list coord. to reciprocal coord.
!
      CALL LATGEN2( LATT,ALPHA,BR1,BR2,ORTHO)
      CALL INVMAT(BR1,DR1)
      CALL INVMAT(BR2,DR2)
      CALL MATMM2(DR2,BR1,DB1)
!     
!
!.....read header in vector file
      NNLO=0
      DO 210 ITY=1,NTYPE
        READ(10) EMIST
        READ(10) ELO
        IF(FL(2)) READ(9) EMIST
        IF(FL(2)) READ(9) ELO
        NLOV(ITY)=NNLO
        NLO(ITY)=0
        DO 211 L = 0,LOMAX
!              write(*,*) 'lo', l,elo(l,1),elo(l,2)
          IF (ELO(L,1).LT.(995.D+0)) THEN
            NNLO=NNLO+((2*L+1))*NATOM(ITY)
          ENDIF
          IF (ELO(L,2).LT.(995.D+0)) THEN
            NNLO=NNLO+((2*L+1))*NATOM(ITY)
            NLO(ITY)=NLO(ITY)+((2*L+1))*NATOM(ITY)
          ENDIF
 211     CONTINUE
 210  CONTINUE
      DO 212 ITY=1,NTYPE
         NLON(ITY)=NNLO-NLO(ITY)-NLOV(ITY)
 212  CONTINUE
      NLO1=NLO(1)+NLON(1)+NLOV(1)
!      write(*,*) 'nlo', nnlo
      nlo1=nnlo
!
!.....output
      WRITE(6,529) TITLE,LATT
      INDEX=0
!      DO 310 ITY=1,NTYPE
!        INDEX=INDEX+1
!        WRITE(6,563) ITY,ANA(ITY),(POS(I,INDEX),I=1,3)
!        DO 312 INA=2,NATOM(ITY)
!          INDEX=INDEX+1
!          WRITE(6,564) (POS(I,INDEX),I=1,3)
!  312   CONTINUE 
!  310 CONTINUE 
!
!......check if case.in1 is consistent with inversion symmetry
       IF((FL(1).AND.FL(3)).OR.((.NOT.FL(1)).AND.(.NOT.FL(3))))THEN
        WRITE(*,*)
        WRITE(6,*)
        IF(FL(3)) THEN
           WRITE(*,*) 'WARNING: ',&
           'case.in1c exists, but inversion symmetry in case.struct'
           WRITE(6,*) 'WARNING: ',&
           'case.in1c exists, but inversion symmetry in case.struct'
         ELSE
           WRITE(*,*) 'WARNING: ',&
           'case.in1 exists, but no inversion symmetry in case.struct'   
           WRITE(6,*) 'WARNING: ',&
           'case.in1 exists, but no inversion symmetry in case.struct'   
         ENDIF
           WRITE(6,*) &
           'Although irrep does not use input from case.in1(c), the program'
           WRITE(6,*) &
           'uses case.in1(c) to determine if complex eigenvectors or not.'
          WRITE(*,*)
          WRITE(6,*)
      ENDIF
!
      IF(FL(2)) FL(1)=.TRUE.
!
!.....write out flags
      IF(     FL1)   WRITE(6,'(A23)',advance='NO') ' Non-symmorphic crystal'
      IF(.NOT.FL1)   WRITE(6,'(A19)',advance='NO') ' Symmorphic crystal'
      IF(     FL(3)) WRITE(6,'(A24)')   ' with inversion symmetry'
      IF(.NOT.FL(3)) WRITE(6,'(A27)')   ' without inversion symmetry'
      IF(     FL(1)) WRITE(6,'(A23)')   ' Complex eigenfunctions'
      IF(.NOT.FL(1)) WRITE(6,'(A20)')   ' Real eigenfunctions'
      IF(     FL(2)) WRITE(6,'(A45)')    &
                   ' Spin-orbit eigenfunctions (->time inversion)'
      IF(.NOT.FL(2)) WRITE(6,'(A29)')    &
                   ' No spin-orbit eigenfunctions'
      IF(     FL(4)) WRITE(6,'(A18)')    &
                   ' Spin-polarization'
      IF(.NOT.FL(4)) WRITE(6,'(A21)')    &
                   ' No spin-polarization'
!
      WRITE(6,590)
      WRITE(6,595) ((BR1(I1,I2),I2=1,3),I1=1,3)
      WRITE(6,591)
      WRITE(6,595) ((DB1(I1,I2),I2=1,3),I1=1,3)
!  
      RETURN
 512  FORMAT(2(3I2,F10.5,/), (3I2,F10.5))
 514  FORMAT(I8)
 529  FORMAT(1X,A40,/,1X,A4,' lattice')
 530  FORMAT(A80)
 541  FORMAT(A4,23X,I3)
 543  FORMAT(12X,F10.7,3X,F10.7,3X,F10.7)
 544  FORMAT(15X,I2)
 545  FORMAT(A10,5X,I5,5X,F10.9,5X,F10.8)
 563  FORMAT(' Atom=',I3,' name=',A10,' position:',3(F12.7))
 564  FORMAT(35X,3(F12.7))
 580  FORMAT(6F10.7,10X,F9.6)
 590  FORMAT(//,' Transformations:',/, &
             ' From k-list in lapw1 to Cartesian coord. system (BR1)')
 591  FORMAT(' From k-list in lapw1 to prim. reciprocal space  (DB1)')
 595  FORMAT(3(3F7.4,/))
!
      END  
