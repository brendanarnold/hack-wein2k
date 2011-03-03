      PROGRAM MIXER                                                     
!                                                                       
!                     F I L E     M I X E R                             
!                     ACTUAL VERSION                                    
!C     THIS ALGORITHM LINKS THE OUTPUT OF LAPW2  AND  CORE    BY         
!       1. CALCULATING NEW TOTAL CHARGE DENSITY AS SUM  OF CLMVAL,      
!          CLMSC AND CLMCORE                                            
!       2. MIXING THE NEW AND OLD DENSITY ACCORDING TO THE PRATT-BROYD      
!          SCHEME                                                       
!       3. INTEGRATING CHARGES IN THE MUFFIN-TIN SPHERES                
!       4. SUMS AND PRINTS THE TOTAL ENERGY                             
!                                                                       
!      PAY ATTENTION TO THE FOLLOWING ITEMS :                           
!          THE MIXING SCHEME IS SET BY A SWITCH:'PRATT' OR 'BROYD'      
!       1. THE MIXING FACTOR WEIGHTS THE NEW CHARGE DENSITY             
!       2. THE LM and K-LIST is updated from clmval
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!
      CHARACTER*19   nkktext
      CHARACTER*4    IREL
      CHARACTER*5    MIX,NORM                                               
      CHARACTER*11   STATUS,FORM                                        
      CHARACTER*20   TITLE                                              
      CHARACTER*60   COMM1,COMM2,COMM3,COMM4                            
      CHARACTER*79   MARGN                                             
      CHARACTER*80   FNAME                                              
      CHARACTER*80   fndmat1,fndmat2,fndmat3,fninc
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN
      LOGICAL        CORY1,CORY2,SCY1,SCY2,VALY1,VALY2,EQIVAL           
      LOGICAL        OLD1,OLD2,cut,itertwice,inversion                                  
      LOGICAL,allocatable ::       FORCE(:)    ! NATO+1)
!      COMPLEX*16        ROKNEW,ROKOLD,ROKSC,ROKVL,ROKMIX                
      complex*16,allocatable :: cfft(:,:,:),ust(:,:,:)
!     >>> Start of additional code for lattice changes
      logical panic
      complex*16,allocatable :: rpanic(:,:) 
      integer,allocatable :: kpanic(:,:) 
!     >>> End of additional code for lattice changes
!                                                                       
      COMMON /SYMETR/  TAU(3,NSYM),IZ(3,3,NSYM),             &
                       INUM(NSYM),IORD                                  
      COMMON /GENER/   BR1(3,3),BR2(3,3),avec(3,3)                    
      dimension RHO(NRAD,2),R(NRAD)
      real*8,allocatable :: CLMOLD(:,:,:,:),CLMNEW(:,:,:,:)
           
      integer,allocatable ::   KVECVL(:,:),KVECSC(:,:),KVCOLD(:,:)            
      complex*16,allocatable :: ROKNEW(:,:),ROKOLD(:,:),ROKSC(:,:),ROKVL(:,:),ROKMIX(:,:)
      real*8,allocatable :: PRATT0(:),RNOT(:),DX(:),ZZ(:)
      integer,allocatable :: JRI(:),ISPLIT(:),LM(:,:,:),LMMAX(:)
      INTEGER, ALLOCATABLE    :: jatom1(:,:,:),ll(:,:,:)
      REAL*8, ALLOCATABLE     :: alx(:,:,:),aly(:,:,:),alz(:,:,:)
      COMPLEX*16, ALLOCATABLE :: dmat(:,:,:,:,:)
      REAL*8, ALLOCATABLE     :: alx_old(:,:,:),aly_old(:,:,:),alz_old(:,:,:)
      COMPLEX*16, ALLOCATABLE :: dmat_old(:,:,:,:,:)
!
!        LATTIC  - lattice type
!        NAME(i) - name of atom i in the unit-cell (e.g. 'Titanium')
!
      CHARACTER*4        LATTIC
      CHARACTER*10,allocatable ::       NAME(:)
      COMMON  /CHAR/     LATTIC
      SAVE    /CHAR/
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        IATNR(i)    - atom index of (inequivalent) atom i
!                      also indicates cubic and non-cubic symmetries
!                      IATNR(i) .GT. 0 ... cubic symmetry
!                      IATNR(i) .LT. 0 ... non-cubic symmetry
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - volume of the direct lattice unit-cell
!
      DOUBLE PRECISION   VI,rotmat(3,3),rotinv(3,3)

      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3)
      INTEGER,allocatable ::  IATNR(:), MULT(:)
      real*8,allocatable :: POS(:,:)
      real*8,allocatable ::  RMT(:), V(:),ROTLOC(:,:,:)
      COMMON  /STRUK/    ALAT, ALPHA, PIA, VI
      SAVE    /STRUK/
      dimension qpw(2),brat(3,3)
      real*8,allocatable ::  FORSUM(:,:),qel(:,:)    
      CHARACTER*10,allocatable ::   ANAME(:)
      real*8,allocatable ::  RHO0(:,:,:),hyperf(:,:,:)
      real*8,allocatable :: fint(:,:),fglob(:,:)              
!---------------------------------------------------------------------  
!                                                                       
      fndmat1='case.dmat_dummy'
      fndmat2='case.dmatup_dummy'
      fndmat3='case.dmatdn_dummy'
      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in MIXER')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
 783  CONTINUE
         READ (1,*,END=784,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         IF(iunit==71) fndmat1=fname
         IF(iunit==72) fndmat2=fname
         IF(iunit==73) fndmat3=fname
         IF(iunit==7) fninc=fname
      GOTO 783
 784  CONTINUE
      CLOSE (1)
!
      cut=.false.
      itertwice=.false.
      VALY1  =.FALSE.                                                   
      SCY1   =.FALSE.                                                   
      CORY1 =.FALSE.                                                    
      VALY2  =.FALSE.                                                   
      SCY2   =.FALSE.                                                   
      CORY2 =.FALSE.                                                    
      OLD1  =.FALSE.                                                    
      OLD2  =.FALSE.
!     Inversion flag
      inversion = .FALSE.                                                    
!                                                                       
!.....READ TAPE20=POTE                                                  
!                                                                       
      READ(20,1000) TITLE                                               
      READ (20,5010) LATTIC, NAT, IREL
      nato=nat
      allocate ( FORCE(NATO+1),NAME(NATO),aname(nato))
      allocate ( CLMOLD(NRAD,NCOM,NATO,2),CLMNEW(NRAD,NCOM,NATO,2))
      allocate ( PRATT0(NATO),RNOT(NATO),DX(NATO),ZZ(NATO))
      allocate ( JRI(NATO),ISPLIT(NATO),LM(2,NCOM,NATO),LMMAX(NATO))
      allocate ( IATNR(NATO), MULT(NATO))
      allocate ( POS(3,NatO*48))
      allocate ( RMT(NATO), V(NATO),ROTLOC(3,3,NATO))
      allocate ( FORSUM(NATO,4),qel(natO,2))    
      allocate ( RHO0(NATO,4,2),hyperf(natO,4,2))
      allocate ( fint(natO,3),fglob(natO,3))              
      ALLOCATE(jatom1(nat,4,3),ll(nat,4,3))
      ALLOCATE(alx(nat,4,3),aly(nat,4,3),alz(nat,4,3))
      ALLOCATE(dmat(nat,4,-3:3,-3:3,3))
      ALLOCATE(alx_old(nat,4,3),aly_old(nat,4,3),alz_old(nat,4,3))
      ALLOCATE(dmat_old(nat,4,-3:3,-3:3,3))
      jatom1=0

      rho0(1:nato,1:4,1:2)=0.0d0
      hyperf(1:nato,1:4,1:2)=0.0d0
      CLMOLD(1:nrad,1:ncom,1:nato,1:2)=0.0d0
      CLMNEW(1:nrad,1:ncom,1:nato,1:2)=0.0d0 
      PI=ACOS(-1.0d0)                                                         
!
!        read in lattice constants
!
      READ (20,5020) ALAT(1), ALAT(2), ALAT(3), &
                     ALPHA(1), ALPHA(2), ALPHA(3)
      IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
      IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
      IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
      call latgen(nat)
!      if(cform.eq.'NEW ') then
        assign 2021 to iform1
        assign 2071 to iform2
        assign 2076 to iform3
!      else
!        assign 2020 to iform1
!        assign 2070 to iform2
!        assign 2075 to iform3
!      end if                                    
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                 
 2070 FORMAT(3X,3I5,2E15.7)                                             
 2071 FORMAT(3X,3I5,2E19.12)                                             
 2075 FORMAT(3X,3I5,4E15.7)                                             
 2076 FORMAT(3X,3I5,2E19.12,2e11.3)                                             
      INDEX=0                                                           
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
!                                                                       
      DO 5 JATOM = 1,NAT                                                
         INDEX=INDEX+1                                                  
!
         READ(20,1040) IATNR(JATOM),POS(1,INDEX),POS(2,INDEX),             &
                       POS(3,INDEX),MULT(JATOM),ISPLIT(JATOM)            
         IF (MULT(JATOM) .EQ. 0) THEN
!
!           illegal number of equivalent atoms
!
            GOTO 930
         ENDIF
            MULTND=MULT(JATOM)-1                                        
            DO 10 M=1,MULTND                                            
            INDEX=INDEX+1                                               
            READ(20,1041) IATNR(JATOM),POS(1,INDEX),POS(2,INDEX),          &
                          POS(3,INDEX)                                   
 10         CONTINUE                                                    
         READ(20,1050) ANAME(JATOM),JRI(JATOM),RNOT(JATOM),RMT(jatom),         &
                       ZZ(JATOM)                                        
         DX(JATOM)=LOG(RMT(jatom)/RNOT(JATOM)) / (JRI(JATOM)-1)               
         READ (20,5060) ((ROTLOC(I1,I2,JATOM),I2=1,3),I1=1,3)
         V(JATOM)=4.D0*PI/3.D0*RMT(JATOM)**3
 5    CONTINUE                                                          
         ndif=index
!                                                                       
!     READ SYMMETRY-OPERATIONS FROM TAPE20=POTE                         
      call nn(nat,mult,rmt,pos)
      READ(20,1055) IORD                                                
      DO 12 J=1,IORD                                                    
      READ(20,1056) ( (IZ(J1,J2,J),J1=1,3),TAU(J2,J), J2=1,3 ),INUM(J)
!     Check if inversion is present
      if((IZ(1,1,J).eq.-1).and.(IZ(2,2,J).EQ.-1).and.(IZ(3,3,J).EQ.-1))then
!     Diagonal elements OK, check others
         if((IZ(2,1,J).EQ.0).and.(IZ(3,1,J).EQ.0).and.(IZ(1,2,J).EQ.0).and. &
          (IZ(3,2,J).EQ.0).and.(IZ(1,3,J).EQ.0).and.(IZ(2,3,J).EQ.0))then
!     Off-diagonals OK, insist on zero translations (may not be always correct)
          if((abs(tau(1,J)).lt.1d-5).and.(abs(tau(2,J)).lt.1d-5).and.(abs(tau(3,J)).lt.1d-5))then
!     We have inversion
               inversion=.true.
          endif
         endif
      endif
 12   CONTINUE
!     ALL OF POTE IS READ                                               
!                                                                       
!.....LOOK, WHICH CHARGE-DENSITY TAPES EXIST                            
!.....READ PRATT-FACTORS; INPUT RHO=RHONEW*PRATT + (1-PRATT)*RHOOLD     
!                                                                       
      READ(17,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) VALY1 =.TRUE.                                      
      READ(18,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) SCY1  =.TRUE.                                      
      READ(19,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) CORY1=.TRUE.                                       
      READ(10,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) OLD1=.TRUE.                                        
      READ(47,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) VALY2 =.TRUE.                                      
      READ(48,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) SCY2  =.TRUE.                                      
      READ(49,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) CORY2=.TRUE.                                       
      IF(VALY2) THEN                                                    
        READ(60,2032,IOSTAT=ITAPE)                                      
        IF(ITAPE.EQ.0) OLD2 =.TRUE.                                     
      ENDIF                                                             
      JSPIN=1                                                           
      IF(VALY2) JSPIN=2                                                 
!                                                                       
!
!     read scf file and sums up total energy and forces
!
      CALL SCFANA(NAT,JSPIN,FORSUM,ESUM,force,iscf,cut,itertwice,mult,fninc)
!
!              
      WRITE(6,1070) ISCF
      NORM='YES'
      supn=0.0d0
      READ(5,702,err=699)MIX,supn,norm
!      READ(5,702)MIX,NORM
 699  continue
      IF (MIX.EQ.'INTER'.OR.MIX.EQ.'PRAT0') THEN
         READ(5,*) (pratt0(i),i=1,nat),pratt
      ELSE 
         READ(5,*)QMX
         nbroyd=20
!        These seem reasonable defaults to use
!        If the charge is oscillating, it may be the Plane waves (not the
!        charge itself; more tests need to be done
         scl1=1.D0
         scl2=1.D0
         nmax=nbroyd+10
         read(5,*,end=698) scl1,scl2
         scl1=scl1/ndif
         read(5,*,end=698) nbroyd,nmax
 698     continue
      ENDIF                                              
      QMX1=1.d0                                                           
      IF(MIX.EQ.'PRATT') QMX1=QMX                                       
!     PRATT'S SCHEME:                                                   
      IF(.NOT.OLD1) THEN                                                
         QMX1=1.d0                                                           
         MIX='PRATT'                                                       
         QMX=1.d0                                                            
         ENDFILE 31                                                        
      ENDIF                                                             
!      WRITE(6,1061) MIX, QMX                                            
!      WRITE(21,1061) MIX, QMX                                           
      IF (MIX.EQ.'PRATT'.OR.MIX.EQ.'BROYD'.or.mix.eq.'BROYG') THEN
         DO 13 JATOM =1,NAT                                              
           PRATT0(JATOM)=QMX                                              
 13      CONTINUE                                                        
         PRATT=QMX1                                            
      ENDIF

      ifilenum=70
      CALL read_denmat(ifilenum,nat,ndm_1,jatom1,ll,alx,aly,alz,dmat)
      ifilenum=73
      CALL read_denmat(ifilenum,nat,ndm_2,jatom1,ll,alx_old,aly_old,alz_old,dmat_old)
      ndm=ndm_2
      IF(ndm_1 /= ndm_2) then
        PRATT=min(0.1d0,qmx)
        DO jatom=1,nat
                pratt0(jatom)=pratt
        enddo
        MIX='PRATT'
      endif
!                                                                       
!.....MIX DENSITY IN SPHERES                                            
!                                                                       
      DO 15 JATOM=1,NAT   
         do 16 jrj=1,jri(jatom)
 16      r(jrj)=rnot(jatom)*exp(dx(jatom)*(jrj-1))
!
         IF(VALY1) THEN                                                 
            READ(17,1980)                                               
            READ(17,2000) LMMAX(JATOM)                                  
            DO 20 LM1=1,LMMAX(JATOM)                                    
               READ(17,2010) LM(1,LM1,JATOM),LM(2,LM1,JATOM)               
               READ(17,iform1) (CLMNEW(J,LM1,JATOM,1),J=1,JRI(JATOM))
       IF(LM1.EQ.1) THEN                                                
       RHO0(JATOM,1,1)=CLMNEW(1,LM1,JATOM,1)/RNOT(JATOM)**2/SQRT(4.*PI) 
       if(jspin.eq.2) call hyper(clmnew(1,1,jatom,1),r, &
       hyperf(jatom,1,1),zz(jatom),dx(jatom),sqrt(4.d0*pi))
       ENDIF                                                            
 20         READ(17,2031)                                               
            READ(17,2033)                                               
         ENDIF                                                          
         IF(VALY2) THEN                                                 
            READ(47,1980)                                               
            READ(47,2000) LMMAX(JATOM)                                  
            DO 520 LM1=1,LMMAX(JATOM)                                   
               READ(47,2031)                                               
               READ(47,iform1) ( CLMNEW(J,LM1,JATOM,2), J=1,JRI(JATOM))      
       IF(LM1.EQ.1) THEN                                                
       RHO0(JATOM,1,2)=CLMNEW(1,LM1,JATOM,2)/RNOT(JATOM)**2/SQRT(4.*PI) 
       if(jspin.eq.2) call hyper(clmnew(1,1,jatom,2),r, &
       hyperf(jatom,1,2),zz(jatom),dx(jatom),sqrt(4.d0*pi))
       ENDIF                                                            
 520        READ(47,2031)                                               
            READ(47,2033)                                               
         ENDIF                                                          
         IF(SCY1.and.MIX.NE.'INTER')  THEN 
            READ(18,1980)                                               
            READ(18,2000) LMMAX(JATOM)                                  
            DO 25 LM1=1,LMMAX(JATOM)                                    
               READ(18,2031)                                               
               READ(18,iform1) ( RHO(J,1), J=1,JRI(JATOM))                   
       IF(LM1.EQ.1) then
            RHO0(JATOM,2,1)=RHO(1,1)/RNOT(JATOM)**2/SQRT(4.*PI)  
            if(jspin.eq.2) call hyper(rho(1,1),r, &
            hyperf(jatom,2,1),zz(jatom),dx(jatom),sqrt(4.d0*pi))
       END IF
            DO 30 J=1,JRI(JATOM)                                        
            CLMNEW(J,LM1,JATOM,1)=CLMNEW(J,LM1,JATOM,1) + RHO(J,1)      
!      IF(LM1.EQ.1)CLMNEW(J,LM1,JATOM)=CLMNEW(J,LM1,JATOM)*SQRT(4.*PI)  
 30         CONTINUE                                                    
 25         READ(18,2031)                                               
            READ(18,2033)                                               
         ENDIF                                                          
         IF(SCY2.and.MIX.NE.'INTER')  THEN   
            READ(48,1980)                                               
            READ(48,2000) LMMAX(JATOM)                                  
            DO 525 LM1=1,LMMAX(JATOM)                                   
            READ(48,2031)                                               
            READ(48,iform1) ( RHO(J,2), J=1,JRI(JATOM))                   
       IF(LM1.EQ.1)then 
            RHO0(JATOM,2,2)=RHO(1,2)/RNOT(JATOM)**2/SQRT(4.d0*PI)  
            if(jspin.eq.2) call hyper(rho(1,2),r, &
            hyperf(jatom,2,2),zz(jatom),dx(jatom),sqrt(4.d0*pi))
       END IF
            DO 530 J=1,JRI(JATOM)                                       
            CLMNEW(J,LM1,JATOM,2)=CLMNEW(J,LM1,JATOM,2) + RHO(J,2)      
!      IF(LM1.EQ.1)CLMNEW(J,LM1,JATOM)=CLMNEW(J,LM1,JATOM)*SQRT(4.*PI)  
 530        CONTINUE                                                    
 525        READ(48,2031)                                               
            READ(48,2033)                                               
         ENDIF                                                          
!                                                                       
         IF(MIX.NE.'INTER') THEN 
            DO 32 J=1,JRI(JATOM)                                           
            DO 32 ISPIN=1,JSPIN                                            
 32      CLMNEW(J,1,JATOM,ISPIN)=CLMNEW(J,1,JATOM,ISPIN)*SQRT(4.d0*PI)    
         ENDIF                                                          
         IF(CORY1.and.MIX.NE.'INTER') THEN 
            READ(19,1980)                                               
            READ(19,2000) LM1                                           
            READ(19,2031)                                               
            READ(19,iform1) ( RHO(J,1), J=1,JRI(JATOM) )                  
            RHO0(JATOM,3,1)=RHO(1,1)/RNOT(JATOM)**2/4./PI/JSPIN         
            if(jspin.eq.2) call hyper(rho(1,1),r, &
            hyperf(jatom,3,1),zz(jatom),dx(jatom),0.5d0)
            DO 35 J=1,JRI(JATOM)                                        
 35           CLMNEW(J,1,JATOM,1)=CLMNEW(J,1,JATOM,1) + RHO(J,1)/JSPIN  
            READ(19,2031)                                               
            READ(19,2033)                                               
!           Look at core charge
            DO J=1,JRI(JATOM)
              R(J)=RNOT(JATOM)*EXP( DX(JATOM)*DBLE(J-1) )
            enddo
            CALL CHARGE(R,DX(JATOM),RHO(1,1),1,JRI(JATOM), QCOR)
            write(6,*)':CINT Core integral, Spin Up atom ',jatom, QCOR/JSPIN
!
         ENDIF                                                          
         IF(CORY2.and.MIX.NE.'INTER') THEN
            READ(49,1980)                                               
            READ(49,2000) LM1                                           
            READ(49,2031)                                               
            READ(49,iform1) ( RHO(J,2), J=1,JRI(JATOM) )                  
            RHO0(JATOM,3,2)=RHO(1,2)/RNOT(JATOM)**2/4./PI/JSPIN         
            if(jspin.eq.2) call hyper(rho(1,2),r, &
            hyperf(jatom,3,2),zz(jatom),dx(jatom),0.5d0)
            DO 535 J=1,JRI(JATOM)                                       
 535           CLMNEW(J,1,JATOM,2)=CLMNEW(J,1,JATOM,2) + RHO(J,2)/JSPIN 
            READ(49,2031)                                               
            READ(49,2033)                                               
!           Look at core charge
            DO J=1,JRI(JATOM)
              R(J)=RNOT(JATOM)*EXP( DX(JATOM)*DBLE(J-1) )
            enddo
            CALL CHARGE(R,DX(JATOM),RHO(1,2),1,JRI(JATOM), QCOR)
            write(6,*)':CINT Core integral, Spin Dn atom ',jatom, QCOR/JSPIN
         ENDIF                                                          
!                                                                       
            DO 501 ISPIN=1,JSPIN                                        
 501   RHO0(JATOM,4,ISPIN)=CLMNEW(1,1,JATOM,ISPIN)/RNOT(JATOM)**2/4./PI 
       IF(OLD1) THEN                                                    
       READ(10,1980)                                                    
       READ(10,2000) LMMAX1
       if(lmmax1.ne.lmmax(jatom)) then
         write(*,*) 'INFO: LM-LIST in CLMSUM changed in MIXER'
         ENDFILE 31
         rewind 31
       endif
!       LMMAX(JATOM)=LMMAX1                     
         DO 40 LM1=1,LMMAX1                                             
            READ(10,2031)                                               
            READ(10,iform1) (CLMOLD(J,LM1,JATOM,1),J=1,JRI(JATOM) )       
 40      READ(10,2031)                                                  
         READ(10,2033)                                                  
       ENDIF                                                            
!                                                                       
       IF(OLD2) THEN                                                    
       IF(JSPIN.EQ.2) THEN                                              
       READ(60,1980)                                                    
       READ(60,2000) LMMAX1                                             
         DO 540 LM1=1,LMMAX1                                            
            READ(60,2031)                                               
            READ(60,iform1) ( CLMOLD(J,LM1,JATOM,2), J=1,JRI(JATOM) )     
540      READ(60,2031)                                                  
         READ(60,2033)                                                  
       END IF                                                           
       END IF                                                           
!                                                                       
  15   CONTINUE 
       DO 3837 ISPIN=1,JSPIN                                            
       IF(ISPIN.EQ.1) MARGN='UP'                                        
       IF(ISPIN.EQ.2) MARGN='DN'                                        
       IF(JSPIN.EQ.1) MARGN='TO'                                        
       WRITE(21,3839)                                                        
 3837  WRITE(21,3838) (MARGN,JATOM,jatom, &
      (RHO0(JATOM,J,ISPIN),J=1,4),JATOM=1,NAT)   
       IF(JSPIN.EQ.2) THEN                                              
       DO 3836 JATOM=1,NAT                                              
       DO 3836 J=1,4                                                    
 3836  RHO0(JATOM,J,1)=RHO0(JATOM,J,1)+RHO0(JATOM,J,2)                  
       WRITE(21,3839)
       MARGN='TO'                                        
       WRITE(21,3838) (MARGN,JATOM,jatom, &
      (RHO0(JATOM,J,1),J=1,4),JATOM=1,NAT)       
 3839  FORMAT(/,7X,'DENSITY AT NUCLEUS'/                                  &
       ,8X,'JATOM',4X,'VALENCE',7X,'SEMI-CORE',12X,'CORE',11X,'TOTAL'/)   
 3838  FORMAT(':R',a2,i3.3,':',1X,I3,4F16.6)   
! Hyperfine fileds with Thomson radius
       do 3833 jatom=1,nat
       hyperf(jatom,4,1)=hyperf(jatom,1,1)+hyperf(jatom,2,1)+ &
       hyperf(jatom,3,1)
       hyperf(jatom,4,2)=hyperf(jatom,1,2)+hyperf(jatom,2,2)+ &
       hyperf(jatom,3,2)
       write(21,3834) jatom
       write(21,3835) (hyperf(jatom,j,1),j=1,4)
       write(21,3835) (hyperf(jatom,j,2),j=1,4)
      write(21,3832)jatom, &
      ((hyperf(jatom,j,1)-hyperf(jatom,j,2))*524.3,j=1,4)
 3833  continue
 3834  format(//,7x,'SPINDENSITIES AT THE NUCLEUS (THOMSON) FOR ', &
       'ATOM',i3,/,10X,'VALENCE',7X,'SEMI-CORE',12X,'CORE',11X,'TOTAL'/)
 3835  format(5x,4f16.6)
 3832  format(':HFF',i3.3,':',1x,'HFF:',f8.3,3f16.3,' (KGAUSS)')
       END IF                                                           
!                                                                  
!.....READ DENSITY IN INTERSTITIAL                                       
!                                                                       
      IF(VALY1)  THEN                                                   
         READ(17,1980)                                                  
!        READ(17,2060) NKKVL
  read(17,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkkvl
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkkvl
 6768 continue
         nwav=nkkvl  
      allocate ( KVECVL(3,NWAV),KVECSC(3,NWAV),KVCOLD(3,NWAV))            
      allocate ( ROKNEW(NWAV,2),ROKOLD(NWAV,2),ROKSC(NWAV,2),ROKVL(NWAV,2),ROKMIX(NWAV,2))
      DO 2 J    =1,NWAV                                                 
      DO 2 ISPIN=1,2                                                    
 2       ROKNEW(J,ISPIN)=(0.0d0,0.0d0)                                      
!         IF(NKKVL.gt.NWAV) goto 900
         DO 140 J=1,NKKVL                                               
 140     READ(17,iform2)  (KVECVL(JX,J),JX=1,3),ROKVL(J,1)                
      ENDIF                                                             
      IF(VALY2)  THEN                                                   
         READ(47,1980)                                                  
!         READ(47,2060) NKKVL                                            
  read(47,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6769) nkkvl
  goto 6770
 6769 read(nkktext,'(13x,i6)') nkkvl
 6770 continue
         DO 5140 J=1,NKKVL                                              
 5140     READ(47,iform2)  (KVECVL(JX,J),JX=1,3),ROKVL(J,2)               
      ENDIF                                                             
      IF(SCY1.and.MIX.NE.'INTER')   THEN   
         READ(18,1980)                                                  
!         READ(18,2060) NKKSC                                            
  read(18,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6771) nkksc
  goto 6772
 6771 read(nkktext,'(13x,i6)') nkksc
 6772 continue
         DO 145 J=1,NKKSC                                               
 145     READ(18,iform2) (KVECSC(JX,J),JX=1,3),ROKSC(J,1)                 
      ENDIF                                                             
      IF(SCY2.and.MIX.NE.'INTER')   THEN  
         READ(48,1980)                                                  
!         READ(48,2060) NKKSC                                            
  read(48,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6773) nkksc
  goto 6774
 6773 read(nkktext,'(13x,i6)') nkksc
 6774 continue
         DO 5145 J=1,NKKSC                                              
5145     READ(48,iform2) (KVECSC(JX,J),JX=1,3),ROKSC(J,2)                 
      ENDIF                                                             
!                                                                       
      IF(OLD1)THEN                                                      
         READ(10,1980)                                                     
!         READ(10,2060) NKKOLD 
  read(10,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6775) nkkold
  goto 6776
 6775 read(nkktext,'(13x,i6)') nkkold
 6776 continue
         if(nkkold.ne.nwav) then
         nkkold=min(nkkold,nwav)          
         write(*,*) 'INFO: K-LIST in CLMSUM changed in MIXER',jvl
         ENDFILE 31
         rewind 31
         endif 
         DO 150 J=1,NKKOLD                                              
 150     READ(10,iform2) (KVCOLD(JX,J),JX=1,3),ROKOLD(J,1)                
      ELSE                                                              
         NKKOLD=NKKVL                                                      
         DO 151 J=1,NKKVL                                                  
            ROKOLD(J,1)=0.d0                                              
            DO 151 JX=1,3                                                     
  151    KVCOLD(JX,J)=KVECVL(JX,J)                                         
      ENDIF                                                             
!                                                                       
      IF(OLD2) THEN                                                     
         IF(JSPIN.EQ.2) THEN                                               
            READ(60,1980)                                                     
!            READ(60,2060) NKKOLD                                              
  read(60,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6777) nkkold
  goto 6778
 6777 read(nkktext,'(13x,i6)') nkkold
 6778 continue
         if(nkkold.ne.nwav)  nkkold=min(nkkold,nwav)          
            DO 5150 J=1,NKKOLD                                             
 5150       READ(60,iform2) (KVCOLD(JX,J),JX=1,3),ROKOLD(J,2)                
         ELSE                                                              
            DO 152 J=1,NKKVL                                                  
               ROKOLD(J,2)=0.d0  
               DO 152 JX=1,3    
  152       KVCOLD(JX,J)=KVECVL(JX,J) 
         ENDIF                                                             
      ENDIF                                                            
!                                                                       
!                                                                       
!.....START SEARCH-LOOP OVER ALL VECTORS OF OLD PLANE WAVE DENSITY      
!                                                                       
!     >>> New version for lattice changes
      panic=.false.
      DO 88 JVL=1,NKKVL
      ROKNEW(jvl,1)=ROKVL(JVL,1)
      IF(VALY2) ROKNEW(jvl,2)=ROKVL(JVL,2)
      IF(SCY1.and.MIX.NE.'INTER')  &
         ROKNEW(jvl,1)=ROKNEW(Jvl,1) + ROKsc(JVL,1)
      IF(SCY2.and.MIX.NE.'INTER') &
         ROKNEW(jvl,2)=ROKNEW(Jvl,2) + ROKsc(JVL,2)
!     Same value ?
      IF(EQIVAL(KVECVL(1,JVL),KVCOLD(1,JVL))) GOTO 88
!     Set panic flag
      if(.not.panic)then
         panic=.true.
         allocate(kpanic(3,nkkold),rpanic(nkkold,2))
         kpanic=kvcold
         rpanic=rokold
!     Eliminate broyden files
        endfile 31
        rewind 31
        endfile 32
        rewind 32
        write(*,*)'Note: k-list has changed'
        write(22,*)':WARNING: K-list has changed'
      endif
!........different K-vector list in clmsum and clmval
!         write(*,*) 'INFO: K-LIST in CLMSUM changed in MIXER',jvl
 188     KVCOLD(1,JVL)=KVECVL(1,JVL)
         KVCOLD(2,JVL)=KVECVL(2,JVL)
         KVCOLD(3,JVL)=KVECVL(3,JVL)
!        Search for equivalent h,k,l
         DO 89 J=Max(JVL-300,1),NKKOLD
            IF((KVECVL(1,JVL) .eq. Kpanic(1,J)) .and. &
               (KVECVL(2,JVL) .eq. Kpanic(2,J)) .and. &
               (KVECVL(3,JVL) .eq. Kpanic(3,J)) )  then
                ROKold(jvl,1)=Rpanic(J,1)
                if(old2) ROKold(jvl,2)=Rpanic(J,2)
                write(6,*)'Found changed K-vector at ',J,jvl
                goto 88
            end if
 89      CONTINUE
!        No h,k,l found
         write(6,*)'Did not find the k vector, set to zero !',jvl
         rokold(jvl,1)=0.d0
         if(old2) rokold(jvl,2)=0.d0
  88     CONTINUE
      if(panic)deallocate(kpanic,rpanic)
!     >>> End of new code for lattice changes                                            
! begin old code
!       DO 88 JVL=1,NKKVL                                              
!      ROKNEW(jvl,1)=ROKVL(JVL,1)                                  
!      IF(VALY2) ROKNEW(jvl,2)=ROKVL(JVL,2)                        
!      IF(SCY1.and.MIX.NE.'INTER')  &
!         ROKNEW(jvl,1)=ROKNEW(Jvl,1) + ROKsc(JVL,1)      
!      IF(SCY2.and.MIX.NE.'INTER') &
!         ROKNEW(jvl,2)=ROKNEW(Jvl,2) + ROKsc(JVL,2)
!      if(jvl.gt.nkkold) goto 188      
!      IF(EQIVAL(KVECVL(1,JVL),KVCOLD(1,JVL))) GOTO 88                  
!!........different K-vector list in clmsum and clmval
!         write(*,*) 'INFO: K-LIST in CLMSUM changed in MIXER',jvl
!         ENDFILE 31
!         rewind 31
! 188     KVCOLD(1,jvl)=KVECVL(1,JVL)                                 
!         KVCOLD(2,jvl)=KVECVL(2,JVL)                                 
!         KVCOLD(3,jvl)=KVECVL(3,JVL)                                 
!         DO 89 J=jvl,NKKOLD                                               
!            IF(EQIVAL(KVECVL(1,JVL),KVCOLD(1,J))) then
!            ROKold(jvl,1)=ROKold(J,1)                                  
!            if(old2) ROKold(jvl,2)=ROKold(J,2)                                 
!            goto 88
!            end if         
! 89      CONTINUE 
!         rokold(jvl,1)=0.d0  
!            if(old2) rokold(jvl,2)=0.d0  
!  88     CONTINUE                                                       
! end old code
      nkkold=nkkvl                                       
      NKKNEW =NKKOLD                                                    
!..................................................................
!.....integrate charges, check normalization of charge density
!.....calc iff parameters
!      call iffpar(ifft1,ifft2,ifft3,iff1,iff2,iff3,nkknew)
      call iffpar(iff1,iff2,iff3,nkknew,kvcold)
      allocate ( cfft(IFF1,IFF2,IFF3), &
            UST((IFF1+1),(IFF2+1),(IFF3+1)))
!.....calc stepfunction
      CALL REAN0a(NKKnew,Kvcold,NAT,IFF1,IFF2,IFF3,UST,mult,rmt,v,pos,inversion)
!
      write(21,*) '    '
      write(6,*) '      CHARGES OF NEW CHARGE DENSITY'
      write(21,*)'      CHARGES OF NEW CHARGE DENSITY'
      margn='N'
      call integr(nkknew,kvcold,roknew,nat,mult,jspin,jri,rnot, &
        dx,r,clmnew,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
      margn='01'
      call normal(roknew,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn &
      ,norm,supn)
!
      if(old1) then
         write(21,*) '    '
         write(6,*) '      CHARGES OF OLD CHARGE DENSITY'
         write(21,*)'      CHARGES OF OLD CHARGE DENSITY'
         margn='O'
         call integr(nkknew,kvcold,rokold,nat,mult,jspin,jri,rnot, &
           dx,r,clmold,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
      margn='02'
         call normal(rokold,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn &
         ,norm,supn)
!........convergence-test
         write(21,*) '    '
         write(21,*)'      CONVERGENCE TEST'
         call distan(nkknew,kvcold,roknew,rokold,rokmix, &
           nat,jspin,jri,rnot,dx,mult,r,rho, &
           clmnew,clmold,iff1,iff2,iff3,qpw,qel,cfft,ust,qtot)
      end if
!...................................................................
!     MIX DENSITIES
!
      IF(MIX.EQ.'INTER'.OR.MIX.EQ.'PRAT0') MIX='PRATT'
      IF(MIX.EQ.'PRATT'.AND.OLD1) THEN
         DO jatom=1,nat
            DO lm1=1,lmmax(jatom)
               DO J=1,JRI(JATOM)                                        
                  CLMNEW(J,LM1,JATOM,1)=CLMNEW(J,LM1,JATOM,1)*PRATT0(JATOM) &
                       + CLMOLD(J,LM1,JATOM,1) * (1.0d0-PRATT0(JATOM))     
               ENDDO
            ENDDO
         ENDDO
         IF(JSPIN.EQ.2.AND.OLD2) THEN      
            DO jatom=1,nat
               DO lm1=1,lmmax(jatom)
                  DO J=1,JRI(JATOM)                                       
                     CLMNEW(J,LM1,JATOM,2)=CLMNEW(J,LM1,JATOM,2)*PRATT0(JATOM) &
                          + CLMOLD(J,LM1,JATOM,2) * (1.0d0-PRATT0(JATOM))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         ! GM: MIX density matrices
         DO idm=1,ndm
            DO  JATOM=1,NAT
               IF(jatom1(jatom,1,idm).EQ.0) EXIT
               DO iorb=1,4
                  IF(jatom1(jatom,iorb,idm).EQ.0) EXIT      
                  alx(jatom,iorb,idm)=alx(jatom,iorb,idm)*pratt0(jatom)+ &
                       alx_old(jatom,iorb,idm)*(1.d0-pratt0(jatom))
                  aly(jatom,iorb,idm)=aly(jatom,iorb,idm)*pratt0(jatom)+ &
                       aly_old(jatom,iorb,idm)*(1.d0-pratt0(jatom))
                  alz(jatom,iorb,idm)=alz(jatom,iorb,idm)*pratt0(jatom)+ &
                       alz_old(jatom,iorb,idm)*(1.d0-pratt0(jatom))
                  l=ll(jatom,iorb,idm)
                  DO m=-l,l
                     DO mp=-l,l
                        dmat(jatom,iorb,m,mp,idm)=dmat(jatom,iorb,m,mp,idm)*pratt0(jatom)+ &
                             dmat_old(jatom,iorb,m,mp,idm)*(1.d0-pratt0(jatom))
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         ! GM.
      ENDIF
      IF(MIX.NE.'PRATT') pratt=1.0d0
      DO J=1,NKKNEW                                                 
         DO ISPIN=1,JSPIN                                              
            ROKMIX(J,ISPIN)=ROKNEW(J,ISPIN)*PRATT+ROKOLD(J,ISPIN)*      &
                 (1.d0-PRATT)                                                
         ENDDO
      ENDDO

      IF(MIX.NE.'PRATT')  then
! ****************************************************
!     Added, April 2006 by LDM
!     Check maximum qel, not just qtot
      qbig=-1.D0
      qrms=0.0D0
      irms=0.D0
      DO JATOM=1,NAT
         DO ISPIN=1,JSPIN
            qbig=max(qbig,qel(JATOM,ISPIN))
            qrms=qrms+qel(JATOM,ISPIN)*qel(JATOM,ISPIN)*MULT(JATOM)
            irms=irms+MULT(JATOM)
         ENDDO
      ENDDO
!     Something here as a trap, not sure exactly what....
!      qbig=qbig/sqrt(dble(NAT))
!     This might also be good, maybe better
      qrms=sqrt(0.5d0*qrms/irms)
      qbig=qbig/3.D0
      write(6,666)'Big check ',qbig,qrms,qtot
      write(21,666)':BIG check (qbig,qrms,qtot) ',qbig,qrms,qtot
666   format(a,3D12.3)
      qtot=max(qtot,qbig,qrms)
!     End of added trap
! ****************************************************
!       Smoother version, somewhat similar, maybe should go down faster
!       Adder term RFACT, was 4, try 3.5, also tail
           RFACT=3.5D0
           RFACT2=1.5d3
           red3=(1.D0-exp(-RFACT2*QTOT*QTOT))
           red2=(1.d0-exp(-RFACT*QTOT))*red3 
           reduction=1.0D0+0.5D0*red2
!          reduction=1.d0
!          if(QTOT.gt.0.01d0) reduction=1.1d0
!          if(QTOT.gt.0.05d0) reduction=1.2d0
!          if(QTOT.gt.0.10d0) reduction=1.3d0
!          if(QTOT.gt.0.30d0) reduction=1.4d0
!          if(QTOT.gt.0.6d0) reduction=1.5d0
!          if(QTOT.gt.1.d0)  reduction=1.6d0
          
          red1=reduction
          if(nat.gt.20)then
                reduction=reduction-(reduction-1)*red2/1.175d0
          else if(nat.ge.15)then
                reduction=reduction-(reduction-1)*red2/1.25d0
          else if(nat.ge.10)then
                reduction=reduction-(reduction-1)*red2/1.5d0
          endif
             do i=1,nat
! Provide a dependence of the mixing factor on the type (3d,4d,4f,5d,5f) and
! number of atoms in the unit cell
            if(zz(i).gt.22.and.zz(i).lt.29) qmx=qmx/reduction/1.05 !1.3d0
            if(zz(i).gt.42.and.zz(i).lt.47) qmx=qmx/reduction !1.3d0
            if(zz(i).gt.56.and.zz(i).lt.72) qmx=qmx/reduction/1.05 !1.3d0
            if(zz(i).gt.76.and.zz(i).lt.79) qmx=qmx/reduction !1.3d0
            if(zz(i).gt.89.and.zz(i).lt.104) qmx=qmx/reduction !1.3d0
             enddo
!          if(QTOT.gt.0.1d0) qmx=qmx/2.d0                                      
!          if(QTOT.gt.0.6d0) qmx=qmx/2.d0                                      
!          if(QTOT.gt.1.d0) qmx=qmx/2.d0                                      
          if(qtot.gt.1d-3)qmx=qmx/(1.d0+qtot*5.d0)
          if(qmx.lt.0.0025d0) qmx=0.0025d0
!        if(MIX.eq.'BROYT')then
!          CALL QMIXT(CLMNEW,CLMOLD,ROKMIX,ROKOLD,   &
!                     JSPIN,NKKNEW,NAT,LMMAX,JRI,QMX, &
!                     ndm,jatom1,ll,alx,aly,alz,dmat, &
!                     alx_old,aly_old,alz_old,dmat_old,mult,inversion, &
!                     scl1,scl2,nmax,LM)
!        else
          write(21,4141) reduction,red1,qmx
4141      format(':REDuction and QMX before broyd:',3f10.4)        
          CALL QMIX5(CLMNEW,CLMOLD,ROKMIX,ROKOLD,   &
                     JSPIN,NKKNEW,NAT,LMMAX,JRI,QMX, &
                     ndm,jatom1,ll,alx,aly,alz,dmat, &
                     alx_old,aly_old,alz_old,dmat_old,mult,inversion, &
                     scl1,scl2,nbroyd,nmax,LM)
!        endif
      ENDIF                                   
      WRITE(6,1061) MIX, QMX                                            
      WRITE(21,1061) MIX, QMX                                           
!                                                                       
      write(21,*) '    '
      write(6,*) '      CHARGES OF MIXED CHARGE DENSITY'
      write(21,*)'      CHARGES OF MIXED CHARGE DENSITY'
         margn='C'
      call integr(nkknew,kvcold,rokmix,nat,mult,jspin,jri,rnot, &
        dx,r,clmnew,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
!
      margn='03'
      call normal(rokmix,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn &
      ,norm,supn)
!.....................................................................
!.....WRITE NEW INPUT DENSITY TO TAPE51,TAPE 52, TAPE11                 
!                                                                       
!     REWIND 11                                                         
      WRITE(51,1970) ISCF                                               
      WRITE(51,78) TITLE,LATTIC,ALAT(1),ALAT(2),ALAT(3),JRI(1)
      WRITE(51,77)                                
      DO 560 JATOM=1,NAT                                                
         WRITE(51,1990) JATOM                                           
         WRITE(51,2001) LMMAX(JATOM)                                    
         DO 565 LM1=1,LMMAX(JATOM)                                      
            WRITE(51,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
            WRITE(51,iform1) ( CLMNEW(J,LM1,JATOM,1), J=1,JRI(JATOM) )    
565         WRITE(51,2031)                                              
560   WRITE(51,2033)                                                    
!                                                                       
      WRITE(51,2051)                                                    
      WRITE(51,1980)                                                    
      WRITE(51,2061) NKKNEW                                             
      DO 5155 J=1,NKKNEW                                                
5155  WRITE(51,iform3) (KVCOLD(JX,J),JX=1,3),ROKMIX(J,1),ROKNEW(J,1)      
!....................................................................
              IF(JSPIN.EQ.2) THEN                                           
!
      WRITE(52,1970) ISCF                                               
      WRITE(52,78) TITLE,LATTIC,ALAT(1),ALAT(2),ALAT(3),JRI(1)
      WRITE(52,77)                               
      DO 660 JATOM=1,NAT                                                
         WRITE(52,1990) JATOM                                           
         WRITE(52,2001) LMMAX(JATOM)                                    
         DO 665 LM1=1,LMMAX(JATOM)                                      
            WRITE(52,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
            WRITE(52,iform1) ( CLMNEW(J,LM1,JATOM,2), J=1,JRI(JATOM) )    
665         WRITE(52,2031)                                              
660   WRITE(52,2033)                                                    
!                                                                       
      WRITE(52,2051)                                                    
      WRITE(52,1980)                                                    
      WRITE(52,2061) NKKNEW                                             
      DO 6155 J=1,NKKNEW                                                
6155  WRITE(52,iform3) (KVCOLD(JX,J),JX=1,3),ROKMIX(J,2),ROKNEW(J,2)      
!                                                                       
      WRITE(11,1970) ISCF                                               
      WRITE(11,78) TITLE,LATTIC,ALAT(1),ALAT(2),ALAT(3),JRI(1)
      WRITE(11,77)                               
      DO 60 JATOM=1,NAT                                                 
         WRITE(11,1990) JATOM                                           
         WRITE(11,2001) LMMAX(JATOM)                                    
         DO 65 LM1=1,LMMAX(JATOM)                                       
            WRITE(11,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
      DO 63 J=1,JRI(JATOM)                                              
 63   CLMNEW(J,LM1,JATOM,1)=CLMNEW(J,LM1,JATOM,1)+CLMNEW(J,LM1,JATOM,2) 
            WRITE(11,iform1) ( CLMNEW(J,LM1,JATOM,1), J=1,JRI(JATOM) )    
 65         WRITE(11,2031)                                              
 60   WRITE(11,2033)                                                    
!                                                                       
      WRITE(11,2051)                                                    
      WRITE(11,1980)                                                    
      WRITE(11,2061) NKKNEW                                             
      DO 155 J=1,NKKNEW                                                 
      ROKMIX(J,1)=ROKMIX(J,1)+ROKMIX(J,2)                               
      ROKNEW(J,1)=ROKNEW(J,1)+ROKNEW(J,2)                               
 155  WRITE(11,iform3) (KVCOLD(JX,J),JX=1,3),ROKMIX(J,1),ROKNEW(J,1)      
      END IF                                                            

      close(71)
      close(72)
      close(73)
      open(71,file=fndmat1,form='formatted')
      open(72,file=fndmat2,form='formatted')
      open(73,file=fndmat3,form='formatted')
      ifilenum=70
      CALL write_denmat(ifilenum,nat,ndm,jatom1,ll,alx,aly,alz,dmat)

! ..................................................................
!                                                                       
      if(itertwice) then
      cut=.true.
      write(21,888) iscf
 888  format(/,':WARN : Iteration ',i3,' occurs more than once in this scf file!!!', &
      /,7x,' ETOT is not calculated correctly !!!!!', &
      /,7x,' You should save your calculations using save_lapw and continue')
      endif

      DO 5719 JATOM=1,NAT
      if(force(jatom)) then
 !check consistency of forces
         ixfor=0
         iyfor=0
         izfor=0
         do lm1=1,lmmax(jatom)
          if(LM(1,LM1,JATOM).eq.1.and.LM(2,LM1,JATOM).eq.1) ixfor=1
          if(LM(1,LM1,JATOM).eq.-1.and.LM(2,LM1,JATOM).eq.1) iyfor=1
          if(LM(1,LM1,JATOM).eq.1.and.LM(2,LM1,JATOM).eq.0) izfor=1
         enddo
         if(ixfor.eq.0.and.abs(forsum(jatom,2)).gt.0.0001d0) then
           write(21,716) jatom,forsum(jatom,2)
           cut=.true.
           forsum(jatom,2)=0.d0
         endif
         if(iyfor.eq.0.and.abs(forsum(jatom,3)).gt.0.0001d0) then
           write(21,717) jatom,forsum(jatom,3)
           cut=.true.
           forsum(jatom,3)=0.d0
         endif
         if(izfor.eq.0.and.abs(forsum(jatom,4)).gt.0.0001d0) then
           write(21,718) jatom,forsum(jatom,4)
           cut=.true.
           forsum(jatom,4)=0.d0
         endif
       endif
 5719 continue
 716     format(':WARN  : X-FORCE for atom ',i3,' is not zero as required by symmetry:',f10.4)   
 717     format(':WARN  : Y-FORCE for atom ',i3,' is not zero as required by symmetry:',f10.4)   
 718     format(':WARN  : Z-FORCE for atom ',i3,' is not zero as required by symmetry:',f10.4)   

      if(cut) then
         WRITE(21,2073) ESUM                                               
         WRITE(6,2073) ESUM                                                
      else
         WRITE(21,2074) ESUM                                               
         WRITE(6,2074) ESUM                                                
      endif
      if(force(NATO+1)) WRITE(21,720)
      if(force(NATO+1)) WRITE(6,720)

      DO 5720 JATOM=1,NAT

      if(force(jatom)) then
         WRITE(21,710) JATOM,JATOM,(FORSUM(JATOM,j),j=1,4)
         WRITE(6,710) JATOM,JATOM,(FORSUM(JATOM,j),j=1,4)
!ad
      do i=1,3
      do j=1,3
      rotmat(i,j)=rotloc(i,j,jatom)
      enddo
      enddo
!
      call invmat(rotmat,rotinv)
!ad
      do i=1,3
      fglob(jatom,i)=0.0
      do j=1,3
      fglob(jatom,i)=fglob(jatom,i)+rotinv(j,i)*forsum(jatom,j+1)
      enddo
      enddo
!ad
      call forcint(pi,br2,brat,alat,alpha,lattic)
!ad
      do i=1,3
      fint(jatom,i)=0.0
      do j=1,3
      fint(jatom,i)=fint(jatom,i)+brat(j,i)*fglob(jatom,j)
      enddo
      enddo
 
      endif
 5720 continue

      if(force(NATO+1))       write(6,721)
      if(force(NATO+1))       write(21,721)

      DO 5721 JATOM=1,NAT
      if(force(jatom)) then
      WRITE(21,711) JATOM,JATOM,(FGLOB(JATOM,j),j=1,3)
      WRITE(6,711) JATOM,JATOM,(FGLOB(JATOM,j),j=1,3)
      endif
 5721 CONTINUE

      if(force(NATO+1))       write(6,722)
      if(force(NATO+1))       write(21,722) 
  
      DO 5722 JATOM=1,NAT
      if(force(jatom)) then
      WRITE(21,712) JATOM,JATOM,(Fint(JATOM,j),j=1,3)
      WRITE(6,712) JATOM,JATOM,(Fint(JATOM,j),j=1,3)
      endif
 5722 CONTINUE

      CALL ERRCLR(ERRFN)
      STOP ' MIXER END'                                                  
!
!        error handling
!
  910 INFO = 1
!
!        'mixer.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('MIXER',ERRMSG)
      GOTO 999
!
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('MIXER',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('MIXER',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('MIXER',ERRMSG)
      GOTO 999
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      WRITE (ERRMSG,9050) JATOM, MULT(JATOM), INDEX
      CALL OUTERR('MIXER',ERRMSG)
      CALL OUTERR('MIXER','MULT .EQ. 0')
      GOTO 999
  940 INFO = 4
      CALL OUTERR('MIXER','Maximum number of K-points exceeded.')
      GOTO 999
!                                                                       
  960 INFO = 7
!
!        Error reading file 'mixer.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('MIXER',ERRMSG)
      GOTO 999
  999 STOP 'MIXER - Error'
!                                                                       
 77   FORMAT(1X,9(F9.7))                                                
 78   FORMAT(1X,A20,A4,3F10.6,I5)                                       
 700  FORMAT(I3,A76)                                                    
 701  FORMAT(I3,A4)                                                     
 702  FORMAT(A5,f8.2,a5)                                                
! 702  FORMAT(2A5,13X,A60)                                                
 703  FORMAT(A8,F8.5,A60)                                               
 704  FORMAT('CONTINUE',F8.5,A60)                                       
 705  FORMAT(F10.3,6X,A60)                                              
 710  FORMAT(':FOR',i3.3,':',1x,i3,'.ATOM',4F15.3)
 711  FORMAT(':FCA',i3.3,':',1x,i3,'.ATOM',15x,3F15.3)
 712  FORMAT(':FGL',i3.3,':',1x,i3,'.ATOM',12x,3F16.9)
 720  FORMAT (7x,'TOTAL FORCE IN mRy/a.u. = |F|',5x,'Fx',13x, &
              'Fy',13x,'Fz')
 721  FORMAT  &
      (/,7x,'TOTAL FORCE WITH RESPECT TO GLOBAL CARTESIAN COORDINATES:')
 722  FORMAT  &
      (/,7x,'TOTAL FORCE WITH RESPECT TO THE GLOBAL COORDINATE SYSTEM:')
 743  FORMAT(3X,A76)                                                    
 1000 FORMAT(A20)                                                       
 1040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1041 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1055 FORMAT(I4)                                                        
 1056 FORMAT(3(3I2,F10.5,/),I8)                                          
 1061 FORMAT(7X,A5,' MIXING SCHEME WITH',F6.3)                          
 1070 FORMAT(I3,'. ITERATION',//)                                       
 1970 FORMAT(3X,'             TOTAL CHARGE DENSITY GENERATED ',        &
             'BY',I3,'. ITERATION ')                 
 1980 FORMAT(3X)                                                        
 1990 FORMAT(3X,'ATOMNUMBER =',I3,5X,10A4)                             
 2000 FORMAT(15X,I3//)                                                  
 2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
 2010 FORMAT(15X,I3,5X,I2/)                                             
 2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                         
 2031 FORMAT(/)                                                         
 2032 FORMAT(//)                                                        
 2033 FORMAT(///)                                                       
 2051 FORMAT(3X,'MIXED TOTAL CHARGE DENSITY IN INTERSTITIAL',            &
            3X,'NEW TOTAL CHARGE DENSITY')                              
 2060 FORMAT(/,13X,I6)                                                  
 2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')                                     
 2073 FORMAT(//,':ENE  :',1X,'*WARNING** TOTAL ENERGY IN Ry =',F20.6/) 
 2074 FORMAT(//,':ENE  :',1X,'********** TOTAL ENERGY IN Ry =',F20.6/) 
 5010 FORMAT(A4,23X,I3,/,13X,A4)
 5060 FORMAT(20X,3F10.7)
 5020 FORMAT(6F10.7)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
 9050 FORMAT('MULT(',I3,'=',I3,', INDEX=',I3)
      END                                                               
