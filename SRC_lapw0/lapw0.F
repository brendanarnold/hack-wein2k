PROGRAM LAPW0
!
!              P O T E N T I A L   G E N E R A T O R
!              ACTUAL VERSION
!
!     LAPW0 CALCULATES THE TOTAL POTENTIAL BY SOLVING THE POISSON-
!     EQUATION  AND USING THE LSD, GGA OR MGGA APPROXIMATIONS
!     FULL NON-MUFFIN-TIN TOTAL POTENTIAL COMPLEX VERSION
!
  USE parallel
  USE defs
  USE ams
  USE struct
  USE efeld
  USE symetr
  USE densit
  USE work
  USE vresp1
  IMPLICIT real*8 (a-h,o-z)
  INCLUDE 'param.inc'
#ifdef Parallel
  include 'mpif.h'
#endif
!
!     THE PARAMETERS MUST SATISFY THE FOLLOWING EQUATIONS :
!         NCOM :     NCOM+3 >= SUM(L,M=0) + 2*SUM(L,M>0)
!                    ... DIMENSION OF YKA, QQ,..
!
  LOGICAL          LR2V, force,force_buf
  LOGICAL,ALLOCATABLE :: logefg(:)
  CHARACTER*19     nkktext
  CHARACTER*11     STATUS,FORM
  CHARACTER*4      SWITCH,R2V
  CHARACTER*5      switc1,switch1,switch2             !EECE switch1 added
  CHARACTER*67     ERRMSG
  CHARACTER*77     RMARGI
  CHARACTER*80     fname,fnavn(99)
  CHARACTER*80     DEFFN, ERRFN
  integer            :: I,J,k,JJ,J1,J2,iunit,irecl,il,lmma,nuc
  integer            :: JSPIN,IXCPOT,ISCF,JATOM,IND,LMMULT1
  integer            :: INDEX,JROW,JCOL,MU,JX,LM1,L,M,NKK,NKK1,LLMM
  integer            :: LL,LMX,NCON,NCON1,IBES,JRJ,LMOD,IRAD,IK,ISPIN
  integer            :: JRE,ITAPS,ITAPN,INFO
  integer            :: IFF1,IFF2,IFF3,II1,II2,II3,IFF1T,IFF2T,IFF3T
  integer            :: iffx2,iffy2,iffz2
  REAL*8             :: MP,bes(0:20),av(3),dfak(42),pi,sqfp,sqrt2,pi4
  REAL*8             :: SUPN,DELCC,QIN,ARG,RK,RK1,d1,rnuc
  REAL*8             :: T947,TE947,TS947,TIMTOT,TMULTI,PMULTI,TCOULI
  REAL*8             :: PCOULI,TXCSPH,PXCSPH,TXCPW,PXCPW,TENERG
  REAL*8             :: PENERG,TREANA,PREANA,TIMEY,TIME0,TIME1,TIME2
  REAL*8             :: TIME3,TIME4,TIME5,TIME6,TCOUL,PCOUL,time4a,time4b
  REAL*8             :: RJ1,RJ2,ACORR,VOLINT,VLIN1,VLIN2,SURFIN
  REAL*8             :: ZSHIFT,VINF,VOUTF,EKVOL1,UST0
  REAL*8             :: ETOTAL
  REAL*8             :: VC              !EECE
  REAL*8             :: c_kub(0:10,0:10)
  REAL*8             ::  vfield(ncom+3,nrad)
  COMMON /norm_kub/ c_kub
!bext
  real*8             :: bextene
  integer            :: imfield,nmod,nsp,natorb
!bext
  REAL*8,ALLOCATABLE :: sigma(:,:)
!  REAL*8,allocatable :: vxc(:,:,:,:),exc(:,:,:,:)
!  real*8,allocatable :: vclmsp(:,:,:,:)
  COMPLEX*16       IMAG,IMAG1,QOUT
  COMPLEX*16       RHOKIN,RHOKT
  COMPLEX*16       CP,CQ(ncom+3),CVOUT(ncom+3)
  COMPLEX*16,ALLOCATABLE   ::  IMAGA(:,:)
  COMPLEX*16,ALLOCATABLE   ::  cfft(:,:,:)
  COMPLEX*16,ALLOCATABLE   ::  fft(:,:,:,:)
  COMPLEX*16,ALLOCATABLE   ::  ust(:,:,:,:)
  
  COMPLEX*16,ALLOCATABLE   ::  potk(:),buf1(:),vkcoul(:),vkexc(:,:), &
       vkvxc(:,:),cvalue(:),taup(:)
  REAL*8,ALLOCATABLE       ::  absk(:),rhopw0(:),rhopw0r(:)
  REAL*8,ALLOCATABLE       ::  fhf(:,:)
  REAL*8,ALLOCATABLE       ::  vtot(:,:,:)

  COMPLEX*16,ALLOCATABLE :: fc(:,:)
  INTEGER,ALLOCATABLE    :: lmmtmx(:),lmmult(:,:,:),lmmax1(:)
  INTEGER,ALLOCATABLE    :: lmmax(:),lm(:,:,:)
  REAL*8,ALLOCATABLE     :: qel(:),vcoul0(:),rldm1(:),rldm2(:)
  COMPLEX*16,ALLOCATABLE :: qq(:,:),qqnew(:,:)

  TYPE efg_type
     REAL*8 :: v20,v21,v22
     REAL*8 :: v21m,v22m
     REAL*8 :: v20srf,v21srf,v22srf
     REAL*8 :: v21msr,v22msr,vzz
     REAL*8 :: winkel,eta
     REAL*8 :: qmat(3,3)
     REAL*8 :: tens(3)
     REAL*8 :: eivec(3,3)
  END TYPE efg_type

  TYPE (efg_type),ALLOCATABLE,TARGET :: efgb(:)
!  TYPE (efg_type), POINTER    :: efgbpointer

#ifdef Parallel
  INTEGER               :: efgb_blocklength(16),efgb_displace(16)
  INTEGER               :: efgb_address(17),efgb_typelist(16)
  INTEGER               :: typempi_efgb
#endif
! EECE additional variables 
        real*8 espxc(2),vcoef
        real*8 ,allocatable:: dexc(:,:),dvxc(:,:)  ! energy tests
        real*8 ,allocatable:: elnum(:),fk(:),elnumsp(:,:),fxck(:,:),fck(:)
        integer ifil
! end EECE

  DATA         IMAG    /(0.D0,1.D0)/
!-----------------------------------------------------------------------
!

  CALL init_ams
  ALLOCATE(yka(ncom+3),r(nrad),value(nrad),rhopw0(nrad+100),rhopw0r(nrad+100))
  ALLOCATE(RLDM1(nrad),rldm2(nrad))
  force_buf = .FALSE.
  dfak(1)=1.d0
  DO i=2,41
     dfak(i)= dfak(i-1)*(2*i-1)
  ENDDO
  pi=4.0d+0*ATAN(1.0d+0)
  sqfp=SQRT(4.d0*pi)
  sqrt2=SQRT(2.D0)
  pi4=4.D0*pi
  
  CALL GTFNAM(DEFFN,ERRFN)
  CALL ERRFLG(ERRFN,'Error in LAPW0')
  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
8000 CONTINUE
  READ (1,*,END=8001,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
#ifdef Parallel
  IF(iunit.EQ.6) THEN 
     il=LEN_TRIM(fname)
     i=myid
     k=1000
     DO j=1,3
        k=k/10
        fname(il+j:il+j)=CHAR(ICHAR("0")+i/k)
        i=i-i/k*k
     ENDDO
  ENDIF
#endif
  fnavn(IUNIT)=FNAME
  OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
  GOTO 8000
8001 CONTINUE
  CLOSE (1)
  CLOSE (10)
  CALL CPUTIM(TIME0)
!      CALL XUFLOW(0)
!
!
!.....READ SWITCH AND IXCPOT
!     SWITCH :  = MULT, COUL, EXCH, POT , TOT
!     IXCPOT :  = 1 :  JMW SPIN-POLARIZED
!               = 2 :  HL-NEW
!               = 3 :  HL-OLD (FREEMAN)
!               = 4 :  VBH
  supn=ZERO
  read(4,'(a5,f8.2)',end=3029,err=3029) switc1,supn 
3029 READ(5,995) SWITCH,IXCPOT
  IF (SWITCH.NE.'MULT'.AND.SWITCH.NE.'COUL'.AND.SWITCH.NE.'EXCH' &
       .AND.SWITCH.NE.'POT ') SWITCH='TOT '
  LR2V=.FALSE.
  switch1=''
  switch2=''
  READ(5,997,END=998) R2V,switch1,switch2           !EECE switch1 added
 998  continue
  faciff=1.d0
  iffx2=0
  iffy2=0
  iffz2=0
  if (switch2 .EQ. ' IFFT') then
    read(5,*,err=1999) iffx2,iffy2,iffz2,faciff
    goto 1998
!    read(5,*) faciff
  endif
 1999  faciff=1.d0
       iffx2=0
       iffy2=0
       iffz2=0
 1998  continue

  if (switch1 .EQ. ' HYBR') then

     if (IXCPOT .EQ. 5) then
        IF(myid.EQ.0)  write(21,*) '      HYBRID LDA-FOCK-ALPHA calculation'
        write(6,*)  '      HYBRID LDA-FOCK-ALPHA calculation'
        IXCPOT=6
     elseif (IXCPOT .EQ. 13) then
        IF(myid.EQ.0)  write(21,*) '      HYBRID PBE-FOCK-ALPHA calculation'
        write(6,*)  '      HYBRID PBE-FOCK-ALPHA calculation'
        IXCPOT=33
     elseif (IXCPOT .EQ. 18) then
        IF(myid.EQ.0)  write(21,*) '      HYBRID B3PW91 calculation'
        write(6,*)  '      HYBRID B3PW91 calculation'
        IXCPOT=6
     endif

  elseif (switch1 .EQ. ' EECE') then

     if (IXCPOT .EQ. 5) then
        IF(myid.EQ.0)  write(21,*) '      LDA-EECE calculation'
        write(6,*)  '      LDA-EECE calculation'
     elseif (IXCPOT .EQ. 13) then
        IF(myid.EQ.0)  write(21,*) '      PBE-EECE calculation'
        write(6,*)  '      PBE-EECE calculation'
     endif

  endif

  if(switch1 .EQ. ' HYBR') switch1=' EECE'

  IF(R2V.EQ.'R2V ')LR2V=.TRUE.
!bext
       imfield=0
       read(7,*,end=5556,err=5556)nmod,nsp,natorb,Bextene
       if(nmod.eq.3)then
       imfield=1
       endif
 5556 continue                        
!bext
!.....READ TAPE21=SCFDATA UNTIL BOTTOM AND INCREASE CYCLENUMBER BY ONE
!
  READ(8,2044)  ISCF
  ISCF=ISCF+1
  IF(ISCF.GT.999) ISCF=1
  if(switch1.ne.' EECE') then
      IF(myid.EQ.0)  WRITE(21,705) ISCF,ISCF
  endif
!bext
       if(imfield.ne.0) then
       IF(myid.EQ.0)  write(21,5557) Bextene
       write(6,5557) Bextene
5557  format('        Interaction with Bext included, muB*Bext=',f12.7,' Ry')
       endif 
!bext
!     Initialize E-field
      call EPOT0  (myid)

!  iefeld=0
!  refeld=0.d0
!!!new efeld by L.Marks
!        read(5,5551,end=5555,err=5555)DEFFN
!5551 format(a)
!        read(deffn,*,end=5552,err=5552) iefeld,refeld,wefeld
!        goto 5553
!5552    wefeld=0.3
!        read(deffn,*,end=5555,err=5555) iefeld,refeld
!5553    continue
!        if(mod(abs(iefeld),1000).gt.2000) stop 'iefeld too large'
!!  read(5,*) (jefeld(i),i=1,abs(iefeld))
!       IF(myid.EQ.0) write(21,5554)refeld, iefeld,wefeld
!       write(6,5554) refeld, iefeld,wefeld
!5554    format('       E-field of ',F10.5,' Ry added for', i5, 'coeff, cut ',F10.5)
! Establish values ramp if needed
!  if(abs(iefeld) .gt. 0)then !999)then
!        iefeld1=mod(abs(iefeld),1000)
!        iefmode=abs(iefeld)/1000
!        call makeback(fefeldpw,fefeldrmt,wefeld,iefeld1,refeld,iefmode)
!  endif
!!!  read(5,*,end=5555,err=5555) iefeld,refeld
!!!  if(mod(abs(iefeld),1000).gt.100) stop 'iefeld too large'
!!!!  read(5,*) (jefeld(i),i=1,abs(iefeld))
!!!       IF(myid.EQ.0) write(21,*) '       E-field of ',refeld,' Ry added for', iefeld, 'coeff'
!!!       write(6,*) '        E-field of ',refeld,' Ry added for', iefeld, 'coeff'
5555 CONTINUE
!
  CALL readstruct
  call haveinversion

  ALLOCATE(kkk(3,nsym))
  ALLOCATE(fc(ncom+3,nat),lmmtmx(nat),lmmult(2,ncom+3,nat))
  ALLOCATE(qel(nat),qq(ncom+3,nat),qqnew(ncom+3,nat))
  ALLOCATE(lmmax1(nat),lmmax(nat),lm(2,ncom+3,nat))
  ALLOCATE(vsph(nat))
  ALLOCATE(v(1:nrad,1:ncom,1:nat))
  ALLOCATE(vcoul2(1:nrad,1:nat))
  ALLOCATE(taup(nsym))
  ALLOCATE(clm(1:nrad,1:ncom,1:nat),clmsp(1:nrad,1:ncom,1:nat,1:2))
  ALLOCATE(logefg(1:nat))
!EECE start
  allocate(dexc(ncom+3,2),dvxc(ncom+3,2))
  allocate(elnum(nat),fk(ncom+3),elnumsp(nat,2),fxck(ncom+3,2),fck(ncom+3))
!EECE stop

  logefg=.false.
  v(1:nrad,1:ncom,1:nat)=zero

!  etotcorr=0.d0
  DO jatom=1,nat
     vsph(jatom)=pi4/3.D0*rmt(jatom)**3
!     etotcorr=etotcorr+supn*zz(jatom)*mult(jatom)
  ENDDO

  CALL LATGEN(NAT)
  WRITE(6,800)
  WRITE(6,805) TITLE
  WRITE(6,810) LATTIC
  WRITE(6,819) IXCPOT
  if(supn.ne.0.d0) write(6,*) ' CHARGED CELL with ',supn
  
  if(lattic(1:1).eq.'H'.or.lattic(1:1).eq.'R') then
     WRITE(6,820) AA,BB,CC
  else
     WRITE(6,820) AA,BB,CC,ALPHA
  endif
  WRITE(6,840) NAT
  WRITE(6,850) IREL
  WRITE(6,860) ISCF
  IF(myid.EQ.0) then  
   WRITE(21,840) NAT
   WRITE(21,805) TITLE
   WRITE(21,810) LATTIC
   WRITE(21,819) IXCPOT
 822 format( ':WARN  :        CHARGED CELL with ',f7.3,/,8x, &
    'an energy correction like C Q**2/(L eps) is not included (PRB51,4014; PRB73,35215)')
   if(supn.gt.1.d-6) write(21,822) supn
   if(supn.lt.-1.d-6) write(21,822) supn
   if(lattic(1:1).eq.'H'.or.lattic(1:1).eq.'R') then
     WRITE(21,820) AA,BB,CC
   else
     WRITE(21,820) AA,BB,CC,ALPHA
   endif
   WRITE(21,821) VOL
   WRITE(21,850) IREL
  endif
  INDEX=0
  DO JATOM=1,NAT
     WRITE(6,1060) ANAME(JATOM)
     DO JROW=1,3
        WRITE(6,1070) ( ROTLOC(JCOL,JROW,JATOM),JCOL=1,3 )
     ENDDO
     DO MU=1,MULT(JATOM)
        INDEX=INDEX+1
        WRITE(6,1080) MU,( POS(JX,INDEX),JX=1,3 )
        DO JROW=1,3
           WRITE(6,1070) (ROTIJ(JCOL,JROW,INDEX),JCOL=1,3)
        ENDDO
     ENDDO
  ENDDO
!
!.....read total charge density CLM from TAPE7=CLMSUM
!     norm of CLM(0,0)=CLM(0,0)*R*R*SQRT(4.*PI), else CLM(0,0)*R*R
!
  DO JATOM=1,NAT
     READ(8,2000) LMMAX(JATOM)
     LMMAX1(JATOM)=LMMAX(JATOM)
     IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
     DO LM1=1,LMMAX(JATOM)
!            WRITE(*,*) JATOM,LMMAX(JATOM),LM1
        READ(8,2010) L,M
        if(l.gt.lmax2) goto 971
        LM(1,LM1,JATOM)=L
        LM(2,LM1,JATOM)=M
        IF(M.NE.0) THEN
           LMMAX1(JATOM)=LMMAX1(JATOM)+1
           IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
           LM(1,LMMAX1(JATOM),JATOM)=L
           LM(2,LMMAX1(JATOM),JATOM)=-M
        END IF
        READ(8,2021) ( CLM(J,LM1,JATOM), J=1,JRI(JATOM))
        READ(8,2031)
     ENDDO
     READ(8,1980)
     READ(8,2030)
  ENDDO

  READ(12,2030,END=3031)
  READ(13,2030,END=3031)
  JSPIN=2
  GOTO 3032
3031 JSPIN=1
3032 CONTINUE

  IF(jspin.EQ.1) THEN
     WRITE(6,*) '      NON-SPINPOLARIZED CALCULATION'
     IF(myid.EQ.0)  WRITE(21,*) '      NON-SPINPOLARIZED CALCULATION'
  ELSE
     WRITE(6,*) '      SPINPOLARIZED CALCULATION'
     IF(myid.EQ.0)  WRITE(21,*) '      SPINPOLARIZED CALCULATION'
  END IF

  DO jatom=1,nat
     IF(jspin.EQ.2) THEN                                               
        READ(12,2000) lmma                                             
        lmmax(jatom)=lmma                  !eece
        LMMAX1(JATOM)=LMMAX(JATOM)         !eece
        READ(13,2000) lmma                                             
        DO lm1=1 ,lmma                                             
           READ(12,2010) l,m                                              
           if(l.gt.lmax2) goto 971         !eece   use lm-list from unit 12
           LM(1,LM1,JATOM)=L
           LM(2,LM1,JATOM)=M
           IF(M.NE.0) THEN
             LMMAX1(JATOM)=LMMAX1(JATOM)+1
             IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
             LM(1,LMMAX1(JATOM),JATOM)=L
             LM(2,LMMAX1(JATOM),JATOM)=-M
           END IF                                !eece end
           READ(12,2021) ( clmsp(j,lm1,jatom,1), j=1,jri(jatom))          
           READ(13,2010) l,m                                              
           READ(13,2021) ( clmsp(j,lm1,jatom,2), j=1,jri(jatom))          
           READ(13,2031)                                                  
           READ(12,2031)                                                  
           if(switch1.eq.' EECE')then     ! EECE start
            do j=1,jri(jatom)
             if (l.eq.0)then
              clmsp(j,lm1,jatom,1)=clmsp(j,lm1,jatom,1)*sqfp
              clmsp(j,lm1,jatom,2)=clmsp(j,lm1,jatom,2)*sqfp
             endif
             clm(j,lm1,jatom)=clmsp(j,lm1,jatom,1)+clmsp(j,lm1,jatom,2)
            enddo
           endif                          !EECE end
        ENDDO
        READ(12,1980)                                                  
        READ(12,2030)                                                  
        READ(13,1980)                                                  
        READ(13,2030)                                                  
     ELSE                                                              
           DO  lm1=1,lmmax(jatom)                                      
        DO  j=1,jri(jatom)                                          
              clmsp(j,lm1,jatom,1)=clm(j,lm1,jatom)*0.5D0                     
              clmsp(j,lm1,jatom,2)=clm(j,lm1,jatom)*0.5D0                     
           ENDDO
        ENDDO
     END IF
  ENDDO

  IF(ixcpot.eq.12.or.ixcpot.eq.45.or.ixcpot.eq.40.or.ixcpot.eq.41.or. &
       ixcpot.eq.42.or.ixcpot.eq.43.or.ixcpot.eq.44.or.ixcpot.eq.27)  THEN
     allocate(vclmsp(1:nrad,1:ncom,1:nat,1:2))
     DO jatom=1,nat
        IF(JSPIN.EQ.2) THEN                                               
           if(jatom.eq.1) then
              READ(29,2030) 
              READ(30,2030)
           endif
           READ(29,2000) LMMA                                             
           READ(30,2000) LMMA                                             
           DO LM1=1 ,LMMA                                             
              READ(29,2010) L,M
              READ(29,2021) ( VCLMSP(J,LM1,JATOM,1), J=1,JRI(JATOM))
              READ(30,2010) L,M
              READ(30,2021) ( VCLMSP(J,LM1,JATOM,2), J=1,JRI(JATOM))
              READ(30,2031)
              READ(29,2031)
           ENDDO
           READ(29,1980)                                                  
           READ(29,2030)                                                  
           READ(30,1980)                                                  
           READ(30,2030)                                                  
        ELSE                                                              
           if(jatom.eq.1) then
              READ(28,2030) 
           END IF
           READ(28,2000) LMMA                                             
           DO LM1=1 ,LMMA                                             
              READ(28,2010) L,M
              READ(28,2021) ( VCLMSP(J,LM1,JATOM,1), J=1,JRI(JATOM))
              READ(28,2031)
           ENDDO
           READ(28,1980)                                                  
           READ(28,2030)                                                  
              DO LM1=1,LMMAX(JATOM)
           DO J=1,JRI(JATOM)                                          
                 VCLMSP(J,LM1,JATOM,1)=VCLMSP(J,LM1,JATOM,1)*0.5D0 
                 VCLMSP(J,LM1,JATOM,2)=VCLMSP(J,LM1,JATOM,1)
              ENDDO
           ENDDO
        END IF
     ENDDO
  END IF
     
#ifdef Parallel
  CALL begend(nat)
#else
  nstart(0)=1
  nstop(0)=nat
  coul=.TRUE.
  xc=.TRUE.
  icoul=0
#endif

  ALLOCATE(efgb(1:nat))
  DO jatom=1,nat
     efgb(jatom)%v20=zero ;efgb(jatom)%v21=zero; efgb(jatom)%v22=zero
      efgb(jatom)%v21m=zero ; efgb(jatom)%v22m=zero
      efgb(jatom)%v20srf=zero ; efgb(jatom)%v21srf=zero ; efgb(jatom)%v22srf=zero
      efgb(jatom)%v21msr=zero ; efgb(jatom)%v22msr=zero ; efgb(jatom)%vzz=zero
      efgb(jatom)%winkel=zero ; efgb(jatom)%eta=zero
      efgb(jatom)%qmat(1:3,1:3)=zero
      efgb(jatom)%tens(3)=zero
      efgb(jatom)%eivec(1:3,1:3)=zero
   ENDDO

  ALLOCATE(sigma(2,1:nat))
  sigma=zero

  ALLOCATE(exc(1:nrad,1:ncom,1:nat,1:2),vxc(1:nrad,1:ncom,1:nat,1:2))
!  allocate(vclmsp(1:nrad,1:ncom,1:nat,1:2))
  vxc=zero
  exc=zero
!  vclmsp=zero

!
!.....read plane wave density RHOK(K) from CLMSUM and average YLM
!
! READ(8,2060) NKK
  read(8,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkk
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkk
 6768 continue
  allocate (vkcoul(1:nkk),vkexc(1:nkk,1:2), &
       vkvxc(1:nkk,1:2),cvalue(nkk),absk(nkk))
  allocate (kzz(3,nkk),inst(nkk))
  allocate (rhok(1:nkk),rhoksp(1:nkk,1:2))
  VKCOUL(1:nkk)=ZEROC
  VKEXC(1:nkk,1:2)=ZEROC
  VKVXC(1:nkk,1:2)=ZEROC
  
  READ(8,2071)  ( KZZ(JX,1), JX=1,3 ), RHOKIN
  RHOK(1)=RHOKIN
  ABSK(1)=0.0
  INST(1)=1
  !
  IFF1=0
  IFF2=0
  IFF3=0
  DO 40 J=2,NKK
     READ(8,2071)  ( KZZ(JX,J), JX=1,3 ), RHOKIN
     RHOK(J)=RHOKIN
9    AV(1)=BR1(1,1)*KZZ(1,J)+BR1(1,2)*KZZ(2,J)+BR1(1,3)*KZZ(3,J)
     AV(2)=BR1(2,1)*KZZ(1,J)+BR1(2,2)*KZZ(2,J)+BR1(2,3)*KZZ(3,J)
     AV(3)=BR1(3,1)*KZZ(1,J)+BR1(3,2)*KZZ(2,J)+BR1(3,3)*KZZ(3,J)
     ABSK(J)=SQRT(AV(1)*AV(1)+AV(2)*AV(2)+AV(3)*AV(3))
     CALL STERN(J,IND,KKK,TAUP)
     INST(J)=IND
     DO 42 JJ=1,IND
!
!     transfer karth into internal coord.
!     find max k-vector
!
        IF(ORTHO) THEN
           II1=NINT((AVEC(1,1)*KKK(1,JJ)+AVEC(1,2)*KKK(2,JJ)+ &
                AVEC(1,3)*KKK(3,JJ))/AA)
           II2=NINT((AVEC(2,1)*KKK(1,JJ)+AVEC(2,2)*KKK(2,JJ)+ &
                AVEC(2,3)*KKK(3,JJ))/BB)
           II3=NINT((AVEC(3,1)*KKK(1,JJ)+AVEC(3,2)*KKK(2,JJ)+ &
                AVEC(3,3)*KKK(3,JJ))/CC)
        ELSE
           II1=KKK(1,JJ)
           II2=KKK(2,JJ)
           II3=KKK(3,JJ)
        END IF
        IF(II1.GT.IFF1) IFF1=II1
        IF(II2.GT.IFF2) IFF2=II2
        IF(II3.GT.IFF3) IFF3=II3
        IF(KKK(1,JJ).GT.IFF1) IFF1=KKK(1,JJ)
        IF(KKK(2,JJ).GT.IFF2) IFF2=KKK(2,JJ)
        IF(KKK(3,JJ).GT.IFF3) IFF3=KKK(3,JJ)
42   CONTINUE
40 CONTINUE
!
!     Set E-field values if needed
      CALL EPOTINIT(NKK,myid)
!
!         jj=0
!!!!      if(abs(iefeld).gt.999) then
!      if(abs(iefeld).gt.0) then
!! select all (0,0,kz)
!         DO J=2,NKK
!!           if(KZZ(1,J).eq.0.and.KZZ(2,J).eq.0.and.KZZ(3,J).ne.0) then
!           if(KZZ(1,J).eq.0.and.KZZ(2,J).eq.0.and.KZZ(3,J).ne.(KZZ(3,J)/2)*2) then
!            jj=jj+1
!            if(jj.gt.mod(abs(iefeld),1000)) goto 43
!            jefeld(jj)=j
!           endif
!         enddo
!!!      else if(iefeld.ne.0) then
!!!! select even (0,0,kz)
!!!         jj=0
!!!         DO J=2,NKK
!!!           if(KZZ(1,J).eq.0.and.KZZ(2,J).eq.0.and.KZZ(3,J).ne.(KZZ(3,J)/2)*2) then
!!!            jj=jj+1
!!!            if(jj.gt.abs(iefeld)) goto 43
!!!             jefeld(jj)=j
!!!           endif
!!!         enddo
!      endif
! 43 continue
!    kefeld=min(jj,mod(abs(iefeld),1000))
!
!     set IFF values for Fast Fourier calculation
!     check dimensions
!
  IFF1=(IFF1+1)*2
  IFF2=(IFF2+1)*2
  IFF3=(IFF3+1)*2
   if(iffx2.gt.0) iff1=iffx2
   if(iffy2.gt.0) iff2=iffy2
   if(iffz2.gt.0) iff3=iffz2
  IFF1=IFF1 * faciff
  IFF2=IFF2 * faciff
  IFF3=IFF3 * faciff
  IFF1T=IFF1 
  IFF2T=IFF2
  IFF3T=IFF3
  CALL IFFLIM(IFF1T,IFF1)
  CALL IFFLIM(IFF2T,IFF2)
  CALL IFFLIM(IFF3T,IFF3)

  allocate (cfft(iff1,iff2,iff3),fft(iff1,iff2,iff3,2))
  allocate (ust(iff1,iff2,iff3,2))
  write(6,*) 'iff1-3,faciff',iff1,iff2,iff3,faciff
  IF(myid.EQ.0)  write(21,'(":IFFT  : FFT-parameters:",3i5," Factor:",f5.2)') iff1,iff2,iff3,faciff

!Check for task parallisme
  IF(.NOT.coul) GOTO 777

!
!.....calculate multipolmoments of true charge density inside spheres from CLM
!
!      DO 20 JATOM=1,NAT
  DO jatom=nstart(myid),nstop(myid)
     DO J=1,JRI(JATOM)
        R(J)=R0(JATOM)*EXP(DX(JATOM)*dble(J-1))
     ENDDO
     DO LM1=1,LMMAX(JATOM)
        L=IABS( LM(1,LM1,JATOM) )
        DO J=1,JRI(JATOM)
           VALUE(J)=CLM(J,LM1,JATOM)*R(J)**L
        ENDDO
        CALL CHARGE(R,DX(JATOM),VALUE,1,JRI(JATOM),QIN)
        QQ(LM1,JATOM)=-QIN
        IF(L.EQ.0) THEN
           QEL(JATOM) =QIN
           QQ(1,JATOM)=QQ(1,JATOM)/SQFP + ZZ(JATOM)/SQFP
!....for charged cells
           delcc=supn/vol
        ENDIF
     ENDDO
  ENDDO
!
!.....write out multipolmoments
!
  WRITE(6,870)

!      DO 35 JATOM=1,NAT
  DO JATOM=NSTART(myid),NSTOP(myid)
     WRITE(6,2040) JATOM,ANAME(JATOM),ZZ(JATOM),LMMAX(JATOM), &
          ( POS(J,JATOM), J=1,3 )
     IF(IATNR(JATOM).GT.0) THEN
        LMMAX1(JATOM)=LMMAX(JATOM)
     ELSE
        IND=LMMAX(JATOM)
        DO 36 LM1=1,LMMAX(JATOM)
           IF(LM(2,LM1,JATOM).NE.0) THEN
              IND=IND+1
              QQ(IND,JATOM)=QQ(LM1,JATOM)
           END IF
36       CONTINUE
     END IF
     CALL MULTFC(FC,JATOM,LMMAX1(JATOM),LM(1,1,JATOM),IATNR(JATOM))
     CALL MULTSU(QQ,fc,JATOM,LMMAX1(JATOM),IATNR(JATOM),lm)
     LMMULT1=0
     DO 38 LM1=1,LMMAX1(JATOM)
!
! combines qq of l,m and -l,m and defines new indices lmmult
!
        LMMULT1=LMMULT1+1
        QQ(LMMULT1,JATOM)=QQ(LM1,JATOM)
        LMMULT(1,LMMULT1,JATOM)=LM(1,LM1,JATOM)
        LMMULT(2,LMMULT1,JATOM)=LM(2,LM1,JATOM)
        IF (LMMULT1.EQ.1) GOTO 38
        IF ((IABS(LMMULT(1,LMMULT1,JATOM)) .EQ. &
             IABS(LMMULT(1,LMMULT1-1,JATOM))) .AND. &
             (LMMULT(2,LMMULT1,JATOM) .EQ. &
             LMMULT(2,LMMULT1-1,JATOM))) THEN
           QQ(LMMULT1-1,JATOM)=QQ(LMMULT1-1,JATOM)+QQ(LMMULT1,JATOM)
           LMMULT1=LMMULT1-1
        END IF
38   CONTINUE
     LMMTMX(JATOM)=LMMULT1

     DO LMMULT1=1,LMMTMX(JATOM)
        WRITE(6,2050) LMMULT(1,LMMULT1,JATOM), &
             LMMULT(2,LMMULT1,JATOM),QQ(LMMULT1,JATOM)
     ENDDO

     IF(IATNR(JATOM).GT.0) THEN
        DO LM1=1,LMMAX(JATOM)
           IF(lm(2,lm1,jatom).NE.0) THEN
              LMMAX1(JATOM)=LMMAX1(JATOM)+1
              IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
              LM(1,LMMAX1(JATOM),JATOM)=lm(1,lm1,jatom)
              LMMULT(1,LMMAX1(JATOM),JATOM)=lm(1,lm1,jatom)
              LM(2,LMMAX1(JATOM),JATOM)=-lm(2,lm1,jatom)
              LMMULT(2,LMMAX1(JATOM),JATOM)=-lm(2,lm1,jatom)
           END IF
        ENDDO
     ENDIF
  ENDDO

!
!.....calculate multipolmoments of plane waves in spheres
!
  INDEX=1
  do jatom=1,NSTART(myid)-1
     INDEX=INDEX+MULT(JATOM)
  enddo
  DO JATOM=NSTART(myid),NSTOP(myid)
     WRITE(6,2040) JATOM,ANAME(JATOM),ZZ(JATOM),LMMAX(JATOM), &
          ( POS(J,INDEX),J=1,3 )
     LLMM=LMMTMX(JATOM)
     DO LM1=1,LLMM
        CQ(LM1)=(0.d0,0.d0)
     ENDDO
     DO J=2,NKK
        ARG=RMT(JATOM)*ABSK(J)
        CALL SPHBES(lmax2+1,ARG,BES)
        CALL YKAV(J,JATOM,YKA,LMMTMX,LMMULT,LMMAX(JATOM))
!           IF(J.LT.10) WRITE(6,3848) JATOM,J,LM1,(YKA(J1),J1=1,9)
! 3848       FORMAT(3I3,9(2F5.2,3X))
        RHOKT=RHOK(J)/ABSK(J)
        DO LM1=1,LLMM
           LL=IABS( lmmult(1,LM1,JATOM) )
           CQ(lm1)=CQ(lm1)+RHOKT*YKA(LM1)*BES(LL+1)
        ENDDO
     ENDDO
     DO LM1=1,LLMM
        LL=IABS( LMMULT(1,LM1,JATOM) )
        QOUT=PI4*RMT(JATOM)**(LL+2)
        IF(LL.EQ.0) CQ(LM1)=CQ(LM1)+RHOK(1)*RMT(JATOM)/3./SQFP
        CQ(LM1)=CQ(LM1)*QOUT*IMAG**LL
        QQNEW(LM1,JATOM)=QQ(LM1,JATOM) + CQ(LM1)
        WRITE(6,996)LMMULT(1,LM1,JATOM),LMMULT(2,LM1,JATOM),CQ(LM1)
     ENDDO
     INDEX=INDEX+MULT(JATOM)
  ENDDO

!
!..... calculate pseudo multipolmoments and define QQs for M<0
!
!krimmel
!      WRITE(57,*) 9
!krimmel
  INDEX=1
  do jatom=1,NSTART(myid)-1
     INDEX=INDEX+MULT(JATOM)
  enddo
!      DO 49 JATOM=1,NAT
  DO JATOM=NSTART(myid),NSTOP(myid)
     LLMM=LMMTMX(JATOM)
     IF(IATNR(JATOM).GT.0) THEN
        DO LM1=1,LMMAX(JATOM)
           IF(lm(2,lm1,jatom).NE.0) THEN
              LLMM=LLMM+1
              IF(llmm.GT.NCOM+3) GOTO 970
              QQNEW(llmm,JATOM)=CONJG( QQNEW(lm1,JATOM) )
           END IF
        ENDDO
     END IF
     WRITE(6,2040) JATOM,ANAME(JATOM),ZZ(JATOM),LMMAX(JATOM), &
          (POS(J,INDEX),J=1,3)
!krimmel
!         WRITE(57,*) LLMM
!krimmel
     DO LMX=1,LLMM
        WRITE(6,994) ( LMMULT(J,LMX,JATOM),J=1,2 ),QQNEW(LMX,JATOM)
     ENDDO
     INDEX=INDEX+MULT(JATOM)
  ENDDO
  CALL CPUTIM(TIME1)


  IF(SWITCH.EQ.'MULT') STOP 'MULT OK'
  IF(SWITCH.EQ.'EXCH') GOTO 445
  RK=0.0
  DO JATOM=1,NAT
     RK1=RMT(JATOM)*ABSK(NKK)
     RK=MAX(RK,RK1)
  ENDDO
!      NCON=RK-5.
!      IF(RK.GT.17.) NCON=RK-6.
!      IF(NCON.LT.8) NCON=8
!      IF(NCON.GT.19) NCON=19
  NCON=9
 1212 WRITE(6,880)
  WRITE(6,1520) NCON,RK

! Buffer for parallel
  allocate(buf1(nkk),potk(nkk))
  buf1(1:nkk)=zeroc
  potk(1:nkk)=zeroc


!
!.....calculate plane wave Coulomb potentiaL
!
! Find some parameters for the next loop
!        DO JATOM=1,NAT
        DO JATOM=NSTART(myid),NSTOP(myid)
        IB_MAX=-999
        L_MAX=-999
        L_MIN=999
        LLMM=LMMTMX(JATOM)
        IF(IATNR(JATOM).GT.0) LLMM=LMMAX1(JATOM)
        DO LM1=1,LLMM
           LL=IABS( LMMULT(1,LM1,JATOM) )
           NCON1=NCON-LL
           IBES=LL+NCON1+1
           IF(IBES.GT.20) GOTO 901
           IF(IBES.GT.IB_MAX)IB_MAX=IBES
           IF(LL.GT.L_MAX)L_MAX=LL
           IF(LL.LT.L_MIN)L_MIN=LL
        ENDDO
        ENDDO
        L_MIN=MIN(0,L_MIN)
        L_MAX=max(L_MAX+2,20)
!        write(*,*)'Lrange ',L_MAX,L_MIN
!        write(*,*)'IB Range ',IB_MAX
!        write(*,*)'Jatom range ',NAT
        allocate(IMAGA(L_MIN:L_MAX,NAT))
!        DO JATOM=1,NAT
        DO JATOM=NSTART(myid),NSTOP(myid)
        DO LL=L_MIN,L_MAX
                NCON1=NCON-LL
                IBES=LL+NCON1+1
                MP=DFAK(IBES+1)/DFAK(LL+1)
                IMAGA(LL,JATOM)=((-IMAG/RMT(JATOM))**LL)*MP*MULT(JATOM)
        ENDDO
        ENDDO

  DO J=2,NKK
!         DO 48 JATOM=1,NAT
     DO JATOM=NSTART(myid),NSTOP(myid)
#ifdef Parallel
        ARG=RMT(JATOM)*ABSK(J)
!        CALL SPHBES(20,ARG,BES)
        CALL SPHBES(IB_MAX,ARG,BES)
#else
!       Do not recalculate bessel functions unless we have to
        JJ=max(1,JATOM-1)
        if( (JATOM.EQ.1) .OR. (ABS(RMT(JJ)-RMT(JATOM)).GT.1D-6))then
                ARG=RMT(JATOM)*ABSK(J)
!               CALL SPHBES(20,ARG,BES)
                CALL SPHBES(IB_MAX,ARG,BES)
        endif
#endif

        LLMM=LMMTMX(JATOM)
        IF(IATNR(JATOM).GT.0) LLMM=LMMAX1(JATOM)
        CALL YKAV(J,JATOM,YKA,LMMTMX,LMMULT,LLMM)
        DO LM1=1,LLMM
           LL=IABS( LMMULT(1,LM1,JATOM) )
           NCON1=NCON-LL
           IBES=LL+NCON1+1
!           IF(IBES.GT.20) GOTO 901
!           CP=(-IMAG)**LL/(RMT(JATOM)**LL)/(ARG**(NCON1+1))* &
!                BES(IBES)*QQNEW(LM1,JATOM)*CONJG(YKA(LM1))
!           MP=DFAK(IBES+1)/DFAK(LL+1)
!           CP=CP*MP
!           buf1(J)=buf1(J)+CP*MULT(JATOM)
           CP=  IMAGA(LL,JATOM)/(ARG**(NCON1+1))* &
                BES(IBES)*QQNEW(LM1,JATOM)*CONJG(YKA(LM1))
            buf1(J)=buf1(J)+CP
        ENDDO
     ENDDO
  ENDDO
#ifdef Parallel
  call MPI_ALLREDUCE(buf1,potk,nkk,MPI_DOUBLE_COMPLEX,MPI_SUM &
       ,comm_coul,ierr)
!      call MPI_ALLREDUCE(buf1,potk,nkk,MPI_DOUBLE_COMPLEX,MPI_SUM &
!           ,MPI_COMM_WORLD,ierr)
#else
  DO j=2, nkk
     potk(j)=buf1(j)
  ENDDO
#endif
  TMP=4.D0*PI/VOL
  TMP1=-8.D0*PI
  DO j =2, nkk
!     POTK(J)=4.D0*PI*POTK(J)/VOL*INST(J)
!     POTK(J)=POTK(J)-RHOK(J)
!     POTK(J)=-8.D0*PI*POTK(J)/ABSK(J)/ABSK(J)
     POTK(J)=(TMP*POTK(J)*INST(J)-RHOK(J))*TMP1/(ABSK(J)*ABSK(J))
     VKCOUL(J)=POTK(J)*0.5D0
  ENDDO
  deallocate(buf1)
  CALL CPUTIM(TIMEy)
  t947=0.d0
!
!.....calculate Coulomb potential inside muffin-tin spheres
!
  allocate(fhf(0:3,1:nat),vcoul0(1:nat))
  fhf = zero ; vcoul0 = zero

  NKK1=NKK-1
  INDEX=1
  do jatom=1,NSTART(myid)-1
     INDEX=INDEX+MULT(JATOM)
  enddo
!      DO 50 JATOM=1,NAT
  DO 50 JATOM=NSTART(myid),NSTOP(myid)
!         INDEX=INDEX + MULT(JATOM)
!         IF(LR2V) WRITE(11,*) 'COUL. POT',LMMAX(JATOM)
     LLMM=LMMAX(JATOM)
     JRJ=JRI(JATOM)
     DO 54 J=1,JRJ
        R(J)=R0(JATOM)*EXP(dble(J-1)*DX(JATOM))
        Rhopw0r(J)=R(J)
54   CONTINUE
     DO J=jrj+1,JRJ+100
        Rhopw0r(J)=R0(JATOM)*EXP(dble(J-1)*DX(JATOM))
     enddo
     CALL CPUTIM(ts947)
!
!.......Coulomb pot at R-Mt
!
     rnuc=2.2677D-05*(atom_mass(int(zz(jatom)))**(1.d0/3.d0))
     write(6,*)'amass, r0:',atom_mass(int(zz(jatom))),rnuc
     DO J=1,JRJ
        if (r(j).ge.rnuc) exit
     ENDDO
     nuc=J
     rnuc=r(j)
     write(6,*)'nuc:',nuc
     rhopw0=0.d0
        IF(SWITCH.EQ.'COUL')  call rhopw(rhok(1),1.d0/sqrt(4.d0*pi),ABSK(1),rhopw0r,jrj,rhopw0)
     DO LM1=1,LLMM
        CVOUT(LM1)=(0.d0,0.d0)
     ENDDO
     DO 55 J=2,NKK
        ARG=RMT(JATOM)*ABSK(J)
        CALL SPHBES(lmax2+1,ARG,BES)
        CALL YKAV(J,JATOM,YKA,LMMTMX,LMMULT,LLMM)
        LMMULT1=0
        DO LM1=1,LLMM
           L=IABS( LM(1,LM1,JATOM) )
           M=LM(2,LM1,JATOM)
!     selects indices of yka
           LMMULT1=LMMULT1+1
           IF (LMMULT1.NE.1) THEN
              IF ((L.EQ.IABS(LMMULT(1,LMMULT1-1,JATOM))).AND. &
                   (M.EQ.LMMULT(2,LMMULT1-1,JATOM))) THEN
                 LMMULT1=LMMULT1-1
              END IF
           END IF
           CVOUT(LM1)=CVOUT(LM1) + POTK(J)*BES(L)*YKA(LMMULT1)
        ENDDO
        IF(SWITCH.EQ.'COUL')  call rhopw(rhok(j),yka,ABSK(J),rhopw0r,jrj,rhopw0)
 55  CONTINUE
!
!        IF(SWITCH.EQ.'COUL') write(3,'(2f10.6)') ( rhopw0r(j),sqfp*rhopw0(j),j=1,jrj+100)
        IF(SWITCH.EQ.'COUL') then
          write(3,'("ATOM",i4,/,(2f15.8))') jatom,( rhopw0r(j),sqfp*rhopw0(j),j=jrj-1,jrj+1)
        endif
     DO LM1=1,LLMM
        M=LM(2,LM1,JATOM)
        IF(M.NE.0) THEN
           IMAG1=(1.d0,0.d0)
           IF(LM(1,LM1,JATOM).LT.0) IMAG1=-IMAG
           IF(MOD(M,2).EQ.1) IMAG1=-IMAG1
           CVOUT(lm1)=CVOUT(lm1)/IMAG1*SQRT2
        ENDIF
     ENDDO
     CALL CPUTIM(te947)
     t947=t947+(te947-ts947)
     if(switch1.eq.' EECE')then         !EECE Exact exchange
! EECE electron numbers
         do ispin=1,2
           DO J=1,JRJ
              VALUE(J)=CLMSP(J,1,JATOM,ISPIN)
           ENDDO
           CALL CHARG2 (R,DX(JATOM),VALUE,1,jrj,VLIN1)
           ACORR=VALUE(1)*R(1)*0.5D0
           elnumsp(jatom,ispin)=vlin1+acorr
         enddo
         elnum(jatom)=elnumsp(jatom,1)+elnumsp(jatom,2)
         write(6,553) jatom,(elnumsp(jatom,ispin),ispin=1,2),elnum(jatom)
553   format(i4,3f12.6,' atom, in sphere electron numbers: up, down, total')
! EECE: calculation of Slater integrals from spherical density
      if(elnum(jatom).gt.0.00001)then ! only atoms selected for EECE
       WRITE(50,*) JATOM,LMMAX(JATOM),' Atom, number of LM components'
       WRITE(51,*) JATOM,LMMAX(JATOM),' Atom, number of LM components'
       do L=0,6,2
        do irad=1,jrj
         do j=1,irad-1
          vaLue(j)=cLm(j,1,jatom)*r(j)**L/r(irad)**(L+1)
         enddo
! r'>r
         do j=irad,jrj
          vaLue(j)=cLm(j,1,jatom)*r(irad)**L/r(j)**(L+1)
         enddo
           CALL CHARG2 (R,DX(JATOM),VALUE,1,IRAD,VLIN1)
           CALL CHARG3 (R,DX(JATOM),VALUE,IRAD,JRJ,VLIN2)
           ACORR=VALUE(1)*R(1)*0.5D0
 if(irad.gt.nrad)write(99,*)' irad,nrad',irad,nrad
 if(l   .lt.1   )write(99,*)'    l,ncom',l,ncom
 if(l   .gt.ncom)write(99,*)'    l,ncom',l,ncom
 if(jatom.gt.nat)write(99,*)'jatom,nat ',jatom,nat
           v(irad,L+1,jatom)=VLIN1+VLIN2+ACORR
        enddo
           DO J=1,JRJ
              VALUE(J)=v(J,L+1,JATOM)*cLm(j,1,jatom)/eLnum(jatom)**2
           ENDDO
           CALL CHARG2 (R,DX(JATOM),VALUE,1,jrj,VLIN1)
           ACORR=VALUE(1)*R(1)*0.5D0
           fk(L+1)=(vLin1+acorr)
        write(6,551)jatom,L,fk(L+1)*2.d0
        IF(myid.EQ.0)  write(21,551)jatom,L,fk(L+1)*2.d0
        write(50,555)jatom,L,fk(L+1)*2.d0
        write(51,555)jatom,L,fk(L+1)*2.d0
551   format(2i4,f12.6,'  atom, L, SLater integraL F^2L, Ry')
555   format(2i4,e20.12,'  atom, L, SLater integraL F^2L, Ry')
      enddo ! end of SLater parameters caLcuLation
     endif
!
     DO LM1=1,LLMM
             CALL potfac (vinf,voutf,Lm(1,Lm1,jatom), &
                            Lm(2,Lm1,jatom),iatnr(jatom))
        L=IABS( LM(1,LM1,JATOM) )
        M=LM(2,LM1,JATOM)
        LMOD=0
        IF(L.EQ.0) LMOD=1
        DO IRAD=1,JRJ
! EECE caLcuLation of Vc in the sphere using 'in sphere rho' onLy
!  r'<r
     do j=1,irad-1
      vaLue(j)=cLm(j,Lm1,jatom)*r(j)**L/r(irad)**(L+1)
     enddo
! r'>r
     do j=irad,jrj
      value(j)=clm(j,lm1,jatom)*r(irad)**l/r(j)**(l+1)
     enddo
           CALL CHARG2 (R,DX(JATOM),VALUE,1,IRAD,VLIN1)
           CALL CHARG3 (R,DX(JATOM),VALUE,IRAD,JRJ,VLIN2)
           ACORR=VALUE(1)*R(1)*0.5D0
           VC=(VLIN1+VLIN2+ACORR)
! EECE: the coefficient is the same as for VOLINT
! for 00 comp. vcoef=2*sqrt(4*pi), coef 2 because v is in Ry
! sqrt(4*pi) because integral over angles (4pi) is multiplied
! by norm of clm_00 (1/sqrt(4pi))
           vcoef =8.0d0*PI/(2*L+1)*VINF
           v(irad,lm1,jatom) =VC*vcoef
     if(jatom.eq.1)then
!!     write(60+lm1,677)r(irad),clm(irad,lm1,jatom),v(irad,lm1,jatom)
677  format(6e13.5)
     endif
    enddo   !instead of 53      CONTINUE
   ENDDO    !end of LM1 loop
   else     !EECE end, original LAPW0 part follows
!
!     v-coul inside spheres
!

     DO LM1=1,LLMM
             CALL potfac (vinf,voutf,lm(1,lm1,jatom), &
                            lm(2,lm1,jatom),iatnr(jatom))
!        CALL POTFAC (VINF,VOUTF,LM1,IATNR(JATOM))
        L=IABS( LM(1,LM1,JATOM) )
        M=LM(2,LM1,JATOM)
        LMOD=0
        IF(L.EQ.0) LMOD=1
        DO J=1,JRJ
                RLDM1(J)=R(J)**L
        ENDDO
        DO J=1,JRJ
                RLDM2(J)=( 1.0D0-(R(J)/RMT(JATOM))**(2*L+1))/R(J)**(L+1)
        ENDDO
        DO 53 IRAD=1,JRJ
!           j1=irad
!           rj1=(1.d0-(r(j1)/RMT(jatom))**(2*l+1))/r(j1)**(l+1)
           rj1=rldm2(irad)
           DO J=1,irad-1
!              J2=J
!              VALUE(J)=CLM(J,LM1,JATOM)*R(J2)**L * rj1
              VALUE(J)=(CLM(J,LM1,JATOM)-delcc*r(j)*r(j)*pi4*LMOD)*RLDM1(J) * rj1
           ENDDO
!           j2=irad
!           rj2=r(j2)**l
           rj2=RLDM1(IRAD)
           DO J=irad,JRJ
!              J1=j
!              VALUE(J)=CLM(J,LM1,JATOM)*RJ2 /R(J1)**(L+1)* &
!                   ( 1.0D0-(R(J1)/RMT(JATOM))**(2*L+1) )
               VALUE(J)=(CLM(J,LM1,JATOM)-delcc*r(j)*r(j)*pi4*LMOD)*RJ2*RLDM2(J)
           ENDDO
!
!     integration in two parts, simpson + 3/8 formula at start(end)
!     for even points
!
           CALL CHARG2 (R,DX(JATOM),VALUE,1,IRAD,VLIN1)
           CALL CHARG3 (R,DX(JATOM),VALUE,IRAD,JRJ,VLIN2)
           ACORR=VALUE(1)*R(1)*0.5D0
           VOLINT=VLIN1+VLIN2+ACORR
           VOLINT =VOLINT*8.0*PI/dble(2*L+1)*VINF
!E-field
!           if(abs(iefeld) .gt. 0)then
!                if(irad .eq. 1.and.lm1.eq.1)call  epot1c(LM,JATOM,LMMAX,vfield)
!                VOLINT=VOLINT+vfield(LM1,IRAD)
!           endif
!
           SURFIN=CVOUT(LM1)*(-4.D0*PI)*(R(IRAD)/RMT(JATOM))**L* &
                IMAG**L*VOUTF
           ZSHIFT=(-2.D0*ZZ(JATOM)/R(IRAD)+ &
                2.D0*ZZ(JATOM)/RMT(JATOM))*SQFP
!               if(irad.eq.jri(jatom)) print*,'lm1,volint',lm1,surfin,VOLINT,vfield(LM1,IRAD)
           V(IRAD,LM1,JATOM)=(VOLINT-SURFIN+LMOD*ZSHIFT)
!           IF(L.EQ.2 .AND. M.EQ.0) write(95,*) r(irad),-V(irad,LM1,JATOM)/r(irad)**2
           IF(IRAD.EQ.nuc) THEN
!
!....  EFG:
!
              IF(L.EQ.2 .AND. M.EQ.0) THEN
                 efgb(jatom)%V20SRF=SURFIN
                 efgb(jatom)%V20=-V(nuc,LM1,JATOM)
              ELSE IF(LM(1,LM1,JATOM).EQ. 2.AND. M.EQ.2) THEN
                 efgb(jatom)%V22SRF=SURFIN
                 efgb(jatom)%V22=-V(nuc,LM1,JATOM)
              ELSE IF(LM(1,LM1,JATOM).EQ.-2.AND. M.EQ.2) THEN
                 efgb(jatom)%V22MSR=SURFIN
                 efgb(jatom)%V22M=-V(nuc,LM1,JATOM)
              ELSE IF(LM(1,LM1,JATOM).EQ.2.AND. M.EQ.1) THEN
                 efgb(jatom)%V21SRF=SURFIN
                 efgb(jatom)%V21 =-V(nuc,LM1,JATOM)
              ELSE IF(LM(1,LM1,JATOM).EQ.-2.AND. M.EQ.1) THEN
                 efgb(jatom)%V21MSR=SURFIN
                 efgb(jatom)%V21M=-V(nuc,LM1,JATOM)
              ENDIF
           ENDIF
           IF(LM1.EQ.1) VCOUL2(IRAD,JATOM)=(VOLINT-SURFIN)/2.D0 + &
                LMOD*(ZZ(JATOM)/RMT(JATOM))*SQFP
53      CONTINUE
!
!fb
!     calculation of ik. component of Hellman-Feynman force
!     from:  lim(r-->0) V_1m / r * grad( r * y_lmp)
!     
!     component of force: 
!     ll,m --> ik
!     +1,1     1
!     -1,1     2
!     +1,0     3
!
        IF (l.EQ.1) THEN
           force_buf=.TRUE.
           ll=lm(1,lm1,jatom)
           ik=MOD(ll+m+1,3)+1
           fhf(ik,jatom)=v(1,lm1,jatom)/r(1) &
                *zz(jatom)*SQRT(3.d0/4.d0/pi)
!!!! correction due to efield
!!!           if(ik.eq.3.and.pos(3,index).ne.0.d0) then
!!!             if(pos(3,index).lt.0.499d0) then
!!!                  fhf(ik,jatom)=fhf(ik,jatom)-2.d0*refeld/cc*2.d0           
!!!             endif
!!!             if(pos(3,index).gt.0.501d0) then
!!!                  fhf(ik,jatom)=fhf(ik,jatom)+2.d0*refeld/cc*2.d0           
!!!             endif
!!!           endif
        ENDIF
!fe
     ENDDO
     vcoul0(jatom)=v(jri(jatom),1,jatom)/sqfp
     WRITE(6,1335)  jatom,aname(jatom),vcoul0(jatom)

!        electric fieldgradient:
!
     IF(lmmax(jatom).GT.1.AND.iatnr(jatom).LT.0) THEN
        CALL EFG(efgb(jatom)%v20,efgb(jatom)%v21,efgb(jatom)%v22 &
             ,efgb(jatom)%v21m,efgb(jatom)%v22m &
             ,efgb(jatom)%v20srf,efgb(jatom)%v21srf,efgb(jatom)%v22srf &
             ,efgb(jatom)%v21msr,efgb(jatom)%v22msr,efgb(jatom)%vzz &
             ,efgb(jatom)%qmat,efgb(jatom)%tens,efgb(jatom)%eivec &
             ,efgb(jatom)%winkel,efgb(jatom)%eta &
             ,rnuc,jatom )
        logefg(jatom)=.TRUE.
     ENDIF
!
    endif    !EECE the original part remained
!   INDEX=INDEX+MULT(JATOM)
50 CONTINUE
!
! Include Efield for pw
!     if(abs(iefeld) .gt. 0) call EPOT3(POTK,nkk,kzz)
!fe
!
! Serial output of coulomb potential
     IF (LR2V) THEN
#ifdef Parallel            
        IF(myid.NE.0) THEN
           CALL MPI_RECV(i,1,MPI_INTEGER,myid-1,MPI_ANY_TAG,comm_coul &
                ,statusmpi,ierr)
        ENDIF
#endif
        
        IF(myid.EQ.0) THEN
           OPEN(10,FILE=fnavn(10))
           WRITE(10,*) '    COULOMB POTENTIAL'
           WRITE(10,*) '    NORM OF VC(L,M)=VC(L,M)*r**2'
           WRITE(10,*)
        else
           OPEN(10,FILE=fnavn(10),POSITION='append')
        ENDIF
        
        DO jatom=nstart(myid),nstop(myid)
           WRITE(10,1990) jatom
           WRITE(10,2001) lmmax(jatom)
           DO  j=1,jri(jatom)
              r(j)=r0(jatom)*EXP(dx(jatom)*(j-1))
           ENDDO
           DO lm1=1,lmmax(jatom)
              WRITE(10,2011) lm(1,lm1,jatom),lm(2,lm1,jatom)
              WRITE(10,2021) (v(j,lm1,jatom)*r(j)*r(j), &
                   j=1,jri(jatom))
              WRITE(10,2031)
           ENDDO
           WRITE(10,2030)
        ENDDO
        IF(myid.EQ.icoul) THEN
!        IF(myid.EQ.0) THEN
           WRITE(10,*) '    COULOMB POTENTIAL IN INTERSTITIAL'
           WRITE(10,*)
           WRITE(10,2061) nkk
           if(inversion)then
           DO j=1,nkk
              WRITE(10,2071) (kzz(jx,j),jx=1,3), dble(potk(j)),0.D0
           ENDDO
           else
           DO j=1,nkk
              WRITE(10,2071) (kzz(jx,j),jx=1,3), potk(j)
           ENDDO
           endif
        ENDIF
        CLOSE(10)
#ifdef Parallel
        IF(myid.LT.icoul) THEN
           CALL MPI_SEND(i,1,MPI_INTEGER,myid+1, 0,comm_coul,ierr)
        ENDIF
#endif
     ENDIF

#ifdef Parallel
     ! Reduce efgb for later output
     efgb_typelist(1:16) = MPI_DOUBLE_PRECISION
     efgb_blocklength(1:13) = 1
     efgb_blocklength(14) = 9
     efgb_blocklength(15) = 3
     efgb_blocklength(16) = 9

     CALL MPI_ADDRESS(efgb(1)%v20,   efgb_address(2),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v21,   efgb_address(3),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v22,   efgb_address(4),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v21m,  efgb_address(5),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v22m,  efgb_address(6),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v20srf,efgb_address(7),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v21srf,efgb_address(8),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v22srf,efgb_address(9),  ierr)
     CALL MPI_ADDRESS(efgb(1)%v21msr,efgb_address(10), ierr)
     CALL MPI_ADDRESS(efgb(1)%v22msr,efgb_address(11), ierr)
     CALL MPI_ADDRESS(efgb(1)%vzz,   efgb_address(12), ierr)
     CALL MPI_ADDRESS(efgb(1)%winkel,efgb_address(13), ierr)
     CALL MPI_ADDRESS(efgb(1)%eta,   efgb_address(14), ierr)
     CALL MPI_ADDRESS(efgb(1)%qmat,  efgb_address(15), ierr)
     CALL MPI_ADDRESS(efgb(1)%tens,  efgb_address(16), ierr)
     CALL MPI_ADDRESS(efgb(1)%eivec, efgb_address(17), ierr)
     
     DO i=1,16
        efgb_displace(i)=efgb_address(i+1)-efgb_address(2)
     ENDDO

     CALL MPI_TYPE_STRUCT(16, efgb_blocklength, efgb_displace, efgb_typelist &
            , typempi_efgb, ierr)

     CALL MPI_TYPE_COMMIT(typempi_efgb, ierr)

     IF(myid.NE.0) THEN
        DO j=nstart(myid),nstop(myid)
           CALL MPI_SEND(efgb(j)%v20, 1, typempi_efgb, 0, 0, comm_coul, ierr)
           CALL MPI_SEND(vcoul0(j), 1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
           CALL MPI_SEND(fhf(0,j), 4, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
           CALL MPI_SEND(logefg(j), 1, MPI_LOGICAL, 0, 0, comm_coul, ierr) 
        ENDDO
     ELSE
        DO i=1,icoul
           DO j=nstart(i),nstop(i)
              CALL MPI_RECV(efgb(j)%v20, 1, typempi_efgb, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
              CALL MPI_RECV(vcoul0(j), 1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
              CALL MPI_RECV(fhf(0,j), 4, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
              CALL MPI_RECV(logefg(j), 1, MPI_LOGICAL, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
           ENDDO
        ENDDO
     ENDIF
#endif 

  IF(SWITCH.EQ.'COUL') THEN
     CALL POISSN (LM,LMMAX)
     STOP ' COUL OK'
  ENDIF
445 CONTINUE
  WRITE(6,2031)
  WRITE(6,1540)
!      TIME2=SECOND()
  CALL CPUTIM(TIME2)
!
!.....calculate exchange potential for MT-spheres and for PW
!
!      WRITE(4) VCOUL2

!Task paralism
  IF(.NOT.xc) GOTO 888

777 CONTINUE
  IF(.NOT.xc) goto 889

!  WRITE(6,*) 'hi Im ',myid,' and Im doing excange'
  XCSPH:  DO JATOM=NSTART(myid),NSTOP(myid)
     JRJ=JRI(JATOM)
     DO J=1,JRJ
        R(J)=R0(JATOM)*EXP(dble(J-1)*DX(JATOM))
     ENDDO
    if(switch1.eq.' EECE')then             !EECE
      if(elnum(jatom).gt.0.00001)then ! only atoms selected for EECE
       CALL XCPOT1(LM,JATOM,LMMAX,LR2V,IXCPOT,JSPIN,sigma(1,jatom),imfield,bextene) !EECE
      endif                                !EECE
    else                                   !EECE
       CALL XCPOT1(LM,JATOM,LMMAX,LR2V,IXCPOT,JSPIN,sigma(1,jatom),imfield,bextene)
    endif                                  !EECE
  ENDDO XCSPH
#ifdef Parallel
  call mpi_reduce(force_buf,force,1,MPI_LOGICAL,MPI_LOR,0,comm_coul,ierr)
  IF(myid.NE.0) THEN
     DO j=nstart(myid),nstop(myid)
        CALL MPI_SEND(sigma(1,j), 2, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
     ENDDO
  ELSE
     DO i=1,icoul
        DO j=nstart(i),nstop(i)
           CALL MPI_RECV(sigma(1,j), 2, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
        ENDDO
     ENDDO
  ENDIF
#else
  force=force_buf
#endif

888 CONTINUE

!!$#ifdef Parallel
!!$!Communicate for taskparallism
!!$  if (icoul.ne.ixc) then
!!$  call MPI_BARRIER(MPI_COMM_WORLD)
!!$  if(myid.ge.nat) then
!!$     call MPI_SEND(vxc(1,1,myid-nat+1,1), nrad*ncom, MPI_DOUBLE_PRECISION,myid-nat,tag,MPI_COMM_WORLD,ierr)
!!$     call MPI_SEND(vxc(1,1,myid-nat+1,2), nrad*ncom, MPI_DOUBLE_PRECISION,myid-nat,tag,MPI_COMM_WORLD,ierr)
!!$     CALL MPI_SEND(exc(1,1,myid-nat+1,1), nrad*ncom, MPI_DOUBLE_PRECISION,myid-nat,tag,MPI_COMM_WORLD,ierr)
!!$     CALL MPI_SEND(exc(1,1,myid-nat+1,2), nrad*ncom, MPI_DOUBLE_PRECISION,myid-nat,tag,MPI_COMM_WORLD,ierr)
!!$!for "this silly functional" transfer also 
!!$!     if(ixcpot.EQ.40) then
!!$!        CALL MPI_SEND(vclmsp(1,1,myid-nat+1,1), nrad*ncom, MPI_DOUBLE_PRECISION,myid-nat,tag,MPI_COMM_WORLD,ierr)
!!$!        CALL MPI_SEND(vclmsp(1,1,myid-nat+1,2), nrad*ncom, MPI_DOUBLE_PRECISION,myid-nat,tag,MPI_COMM_WORLD,ierr) 
!!$!     endif
!!$     goto 889
!!$  else
!!$     call MPI_RECV(vxc(1,1,myid+1,1), nrad*ncom, MPI_DOUBLE_PRECISION,myid+nat,tag,MPI_COMM_WORLD,STATUSMPI,ierr)
!!$     call MPI_RECV(vxc(1,1,myid+1,2), nrad*ncom, MPI_DOUBLE_PRECISION,myid+nat,tag,MPI_COMM_WORLD,STATUSMPI,ierr)
!!$     call MPI_RECV(exc(1,1,myid+1,1), nrad*ncom, MPI_DOUBLE_PRECISION,myid+nat,tag,MPI_COMM_WORLD,STATUSMPI,ierr)
!!$     call MPI_RECV(exc(1,1,myid+1,2), nrad*ncom, MPI_DOUBLE_PRECISION,myid+nat,tag,MPI_COMM_WORLD,STATUSMPI,ierr)
!!$!     if(ixcpot.EQ.40) then
!!$!        CALL MPI_RECV(vclmsp(1,1,myid+1,1), nrad*ncom, MPI_DOUBLE_PRECISION,myid+nat,tag,MPI_COMM_WORLD,ierr)
!!$!        CALL MPI_RECV(vclmsp(1,1,myid+1,2), nrad*ncom, MPI_DOUBLE_PRECISION,myid+nat,tag,MPI_COMM_WORLD,ierr) 
!!$!     endif
!!$  endif
!!$  endif
!!$#endif


  CALL CPUTIM(TIME3)
  
  if(switch1.eq.' EECE') goto 3456  !EECE jump over interstitial

!        array UST is used only as working space in XCPOT3

  CALL XCPOT3 (NKK,VKEXC,VKVXC,LR2V,IXCPOT,JSPIN, &
       IFF1,IFF2,IFF3,CFFT,UST,FFT,POTK,imfield,bextene)
  IF(SWITCH.EQ.'EXCH') STOP 'EXCH OK'
  CALL CPUTIM(TIME4)
!     
!     shift to mt-zero
!
  CALL REAN0(NKK,KZZ,IFF1,IFF2,IFF3,UST)
  ISPIN=1
  CALL CPUTIM(TIME4a)
  CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT,UST,UST0)
  EKVOL1=CVALUE(1)*UST0
  CALL CPUTIM(TIME4b)
  IF(SWITCH.EQ.'TOT ') &
        CALL ENERGY(VCOUL0,QEL,LM,LMMAX,ekvol1,ixcpot,etotal)
  CALL CPUTIM(TIME5)
!
!.....transform plane wave total potential to zero inside spheres
!     (stepfunction)
!     CVALUE(K) is now the reanalyzed  PW  potential  coefficient
!
  CALL REAN4 (NKK,IFF1,IFF2,IFF3,UST,FFT,VKeXC)
!  
! sum up the potentials and do the output
!
  SCFWRITE:  IF(myid.EQ.0) THEN
!     OPEN(21,file=fnavn(21),POSITION='append')
     IF(myid.EQ.0) then
     DO jatom=1,nat
        WRITE(21,1335) jatom,aname(jatom),vcoul0(jatom)
        IF(logefg(jatom)) THEN
           WRITE(21,1340)  jatom,efgb(jatom)%vzz
           WRITE(21,1341)  efgb(jatom)%V20,efgb(jatom)%V20SRF
           WRITE(21,1342)  efgb(jatom)%V22,efgb(jatom)%V22SRF
           WRITE(21,1343)  efgb(jatom)%V22M,efgb(jatom)%V22MSR
           WRITE(21,1344)  efgb(jatom)%V21,efgb(jatom)%V21SRF
           WRITE(21,1345)  efgb(jatom)%V21M,efgb(jatom)%V21MSR
           WRITE(21,1350)  efgb(jatom)%qmat(1,1),efgb(jatom)%qmat(1,2), &
                efgb(jatom)%qmat(1,3),efgb(jatom)%TENS(1),0.,0.
           WRITE(21,1350)   efgb(jatom)%qmat(2,1),efgb(jatom)%qmat(2,2), &
                efgb(jatom)%qmat(2,3),0.,efgb(jatom)%TENS(2),0.
           WRITE(21,1350)   efgb(jatom)%qmat(3,1),efgb(jatom)%qmat(3,2), &
                efgb(jatom)%qmat(3,3),0.,0.,efgb(jatom)%TENS(3)
           WRITE(21,1355)   ((efgb(jatom)%eivec(i,j),j=1,3),i=1,3), &
                jatom,efgb(jatom)%winkel
           WRITE(21,1360)  jatom,efgb(jatom)%eta
        ENDIF
     ENDDO
     IF (force) THEN
!       Efeld term
        if(abs(iefeld) .gt. 0) call EFORCE(FHF)
        DO jatom=1,nat
           WRITE(21,710)
           fhf(0,jatom)=SQRT (fhf(1,jatom)**2 +fhf(2,jatom)**2 &
                +fhf(3,jatom)**2)
           WRITE(21,711)jatom,jatom,(fhf(ik,jatom)*1.d3,ik=0,3)
        ENDDO
     ENDIF
     DO jatom=1,nat
        WRITE(21,699) JATOM,JATOM,SIGMA(1,jatom)
        WRITE(21,699) JATOM,JATOM,SIGMA(2,jatom)
     ENDDO
     
!     write(6,*) 'total energy correction due to charge cell (q*Z/2)',etotcorr/2.d0
!     write(21,*) '          total energy correction due to charge cell (q*Z/2)',etotcorr/2.d0
     WRITE(6,1121) ETOTAL      
     WRITE(21,1121) ETOTAL      
! Add efeld terms for convenience here to potk (not used in fact!)
      if(abs(iefeld).gt.0)call epot3(potk,nkk)
! calculate potential at z=0 for workfunction (requires origin in interstitial)
      call workf1(kzz,potk,vkvxc,nkk) 
     endif
  ENDIF SCFWRITE

 3456 continue     !EECE

  allocate(VTOT(NRAD,NCOM,nstart(myid):nstop(myid)))
  SPIN:  DO ISPIN=1,JSPIN
     if(switch1.eq.' EECE') goto 3457          !EECE
!     
     CVALUE(1:NKK)=(0.d0,0.d0)
     TMP=1.D0/VOL
     DO JRE=1,NKK
        CVALUE(JRE)=VKeXC(JRE,ISPIN)*tmp
     ENDDO
            
     DO JATOM=NSTART(myid),NSTOP(myid)
        DO J=1,JRI(JATOM)
           R(J)=R0(JATOM)*EXP(DX(JATOM)*dble(J-1))
        ENDDO
        if(abs(iefeld) .gt. 0)call epot1c(LM,JATOM,LMMAX,vfield)
        DO LM1=1,LMMAX(JATOM)
           DO J=1,JRI(JATOM)
!     E-field, added completely after the energy terms
           if(abs(iefeld) .gt. 0)then
              VTOT(J,LM1,JATOM) = V(J,LM1,JATOM) + &
                   VXC(J,LM1,JATOM,ISPIN) +vfield(LM1,J)
           else
              VTOT(J,LM1,JATOM) = V(J,LM1,JATOM) + &
                   VXC(J,LM1,JATOM,ISPIN)
           endif
           ENDDO
        ENDDO
     ENDDO

      close(9)
     close(10)
     close(11)
     close(16)
     close(17)
     close(18)
     close(19)
     
#ifdef Parallel            
     IF(myid.NE.0) THEN
        CALL MPI_RECV(i,1,MPI_INTEGER,myid-1, MPI_ANY_TAG,comm_coul &
             ,statusmpi,ierr)
     ENDIF
#endif


     
     IF (LR2V) THEN
        IF(myid.EQ.0) THEN
           if(ispin.eq.1) then
             OPEN(9,FILE=fnavn(9))
             OPEN(11,file=fnavn(11))
           else
             OPEN(9,FILE=fnavn(9)(1:LEN_TRIM(fnavn(9)))//'dn')
             OPEN(11,file=fnavn(11)(1:LEN_TRIM(fnavn(11)))//'dn')
           endif
           WRITE(9,*) '    TOTAL POTENTIAL'
           WRITE(9,*) '    NORM OF VC(L,M)=VC(L,M)*r**2'
           WRITE(9,*)
           WRITE(11,1970) ISCF
           WRITE(11,78)
           WRITE(11,77)
        else
           if(ispin.eq.1) then
           OPEN(9,FILE=fnavn(9),POSITION='append')
           OPEN(11,file=fnavn(11),POSITION='append')
           else
             OPEN(9,FILE=fnavn(9)(1:LEN_TRIM(fnavn(9)))//'dn',POSITION='append')
             OPEN(11,file=fnavn(11)(1:LEN_TRIM(fnavn(11)))//'dn',POSITION='append')
           endif
        ENDIF
        
        DO jatom=nstart(myid),nstop(myid)
           WRITE(9,1990) JATOM
           WRITE(9,2001) LMMAX(JATOM)
           WRITE(11,1990) JATOM                                        
           WRITE(11,2001) LMMAX(JATOM)                                 
           DO  J=1,JRI(JATOM)
              R(J)=R0(JATOM)*EXP(DX(JATOM)*dble(J-1))
           ENDDO
           DO LM1=1,LMMAX(JATOM)
              WRITE(9,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
              WRITE(9,2021) (VTOT(J,LM1,JATOM)*R(J)*R(J), &
                   J=1,JRI(JATOM))
              WRITE(9,2031)
              WRITE(11,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
              WRITE(11,2021) (VXC(J,LM1,JATOM,ispin)*r(j)*r(j) &
                   ,J=1,JRI(JATOM))
              WRITE(11,2031)                                           
           ENDDO
           WRITE(9,2030)
           WRITE(11,2030)
        ENDDO
        IF(myid.EQ.icoul) THEN
!        IF(myid.EQ.0) THEN
           WRITE(9,*) '    TOTAL POTENTIAL IN INTERSTITIAL'
           WRITE(9,*)
           WRITE(9,2061) nkk
           WRITE(11,*) '    EXCHANGE-CORRELATION POTENTIAL IN INTERSTITIAL'
           WRITE(11,*)
           WRITE(11,2061) nkk
           if(inversion)then
           DO j=1,nkk
              WRITE(9,2071) (kzz(jx,j),jx=1,3), dble(potk(j)+vkvxc(j,ispin)),0.D0
              WRITE(11,2071) (kzz(jx,j),jx=1,3), dble(vkvxc(j,ispin)),0.D0
           ENDDO
           else
           DO j=1,nkk
              WRITE(9,2071) (kzz(jx,j),jx=1,3), potk(j)+vkvxc(j,ispin)
              WRITE(11,2071) (kzz(jx,j),jx=1,3), vkvxc(j,ispin)
           ENDDO
           endif
        ENDIF
      END IF
      close(9)
      close(11)
!
!.....write out spherical part of total potential to TAPE 16,17 = VSP ,
!     PW and nonspherical part of total potential to TAPE 18,19 = VNS
!
 3457 continue
     if(switch1.ne.' EECE')then          !EECE
      ITAPS=15+ISPIN
      ITAPN=17+ISPIN

      IF(myid.EQ.0) THEN
         open(itaps,file=fnavn(itaps))
         open(itapn,file=fnavn(itapn))
         WRITE(itaps,787) iscf
         WRITE(itaps,*) '    NORM:  V*R'
         WRITE(itaps,*)
         WRITE(itapn,*) '    TOTAL NONSPHERICAL POTENTIAL OF MT-SPHERES'
         WRITE(itapn,*) '    NORM  OF VLM(R) = VLM(R)'
         WRITE(itapn,*)
      else
         open(itaps,file=fnavn(itaps),POSITION='append')
         open(itapn,file=fnavn(itapn),POSITION='append')
      ENDIF
     endif          !EECE
      !
      DO jatom=nstart(myid),nstop(myid)
       if(switch1.ne.' EECE')then          !EECE
         WRITE(itaps,1990) jatom
         WRITE(itapn,1990) jatom
         WRITE(itaps,2001) 1
         WRITE(itapn,2001) lmmax(jatom)-1
       endif                              !EECE
!
      if(switch1.eq.' EECE')then     !EECE potential output
          DO J=1,JRI(JATOM)
           R(J)=R0(JATOM)*EXP(DX(JATOM)*(J-1))
            do k=1,3
             ifil=69+k+3*(ispin-1)
!!             write(ifil,678)r(j),v(j,k,jatom),vxc(j,k,jatom,ispin),vtot(j,k,jatom)
678  format(4f12.6)
            enddo
          ENDDO
     if(elnum(jatom).gt.0.00001d0)then    ! EECE start

       IF (IATNR(JATOM) .LT. 0) THEN     !non-cubic case

          espxc(ispin)=0.
          DO LM1=1,LMMAX(JATOM)
!EECE espxc=Exc-Vxc integrated over the sphere
           do j=1,jri(jatom)
            value(j)=clmsp(j,lm1,jatom,ispin)*  &
              (exc(j,lm1,jatom,ispin)-vxc(j,lm1,jatom,ispin))
            if(lm1.eq.1)value(j)=value(j)/sqfp
           enddo
            CALL CHARG2 (R,DX(JATOM),VALUE,1,jrj,VLIN1)
           ACORR=VALUE(1)*R(1)*0.5d0
           espxc(ispin)= espxc(ispin)+vlin1+acorr
! tests of energy
           do j=1,jri(jatom)
            value(j)=clmsp(j,lm1,jatom,ispin)*  &
              exc(j,lm1,jatom,ispin)
            if(lm1.eq.1)value(j)=value(j)/sqfp
           enddo
            CALL CHARG2 (R,DX(JATOM),VALUE,1,jrj,VLIN1)
           ACORR=VALUE(1)*R(1)*0.5d0
           dexc(lm1,ispin)= vlin1+acorr
           do j=1,jri(jatom)
            value(j)=clmsp(j,lm1,jatom,ispin)*  &
              vxc(j,lm1,jatom,ispin)
            if(lm1.eq.1)value(j)=value(j)/sqfp
           enddo
            CALL CHARG2 (R,DX(JATOM),VALUE,1,jrj,VLIN1)
           ACORR=VALUE(1)*R(1)*0.5d0
           dvxc(lm1,ispin)= vlin1+acorr
! tests of energy end
!  EECE calculation of <g|vxc|g> radial part
           do j=1,jri(jatom)
            value(j)=vxc(j,lm1,jatom,ispin)*clmsp(j,1,jatom,ispin)
           enddo
            CALL CHARG2 (R,DX(JATOM),VALUE,1,jrj,VLIN1)
           ACORR=VALUE(1)*R(1)*0.5D0
           fxck(lm1,ispin)=(vlin1+acorr)/elnumsp(jatom,ispin)
           write(6,652)jatom,lm1,(lm(i,lm1,jatom),i=1,2),   &
           fxck(lm1,ispin),espxc(ispin)
           write(49+ispin,653)jatom,lm1,(lm(i,lm1,jatom),i=1,2), &
           fxck(lm1,ispin),dexc(lm1,ispin),dvxc(lm1,ispin)
652  format(4i4,2f12.6,' atom, lm1, l, m, fxc, Exc-Vxc')
653  format(4i4,3e20.12)
654  format(i4,e20.12)
655  format(i4,f12.6,' atom, XC energy')
         ENDDO
           write(6,655)jatom,espxc(ispin)
           write(49+ispin,654)jatom,espxc(ispin)

       ELSEIF (IATNR(JATOM) .GT. 0) THEN     !cubic case

         ESPXC(ISPIN) = 0D0
         DEXC(:,ISPIN) = 0D0
         DVXC(:,ISPIN) = 0D0
         FXCK(:,ISPIN) = 0D0

         DO 41 LM1=1, LMMAX(JATOM)

           L = IABS(LM(1,LM1,JATOM))
           M = LM(2,LM1,JATOM)

           IF ((L .EQ. 0) .OR. (L .EQ. 3)) THEN

             DO J=1, JRI(JATOM)
               VALUE(J) = CLMSP(J,LM1,JATOM,ISPIN)* &
                 (EXC(J,LM1,JATOM,ISPIN) - VXC(J,LM1,JATOM,ISPIN))
               IF (L .EQ. 0) VALUE(J) = VALUE(J)/SQFP
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             ESPXC(ISPIN) = ESPXC(ISPIN) + VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = CLMSP(J,LM1,JATOM,ISPIN)*EXC(J,LM1,JATOM,ISPIN)
               IF (L .EQ. 0) VALUE(J) = VALUE(J)/SQFP
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             DEXC(LM1,ISPIN) = VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = CLMSP(J,LM1,JATOM,ISPIN)*VXC(J,LM1,JATOM,ISPIN)
               IF (L .EQ. 0) VALUE(J) = VALUE(J)/SQFP
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             DVXC(LM1,ISPIN) = VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = VXC(J,LM1,JATOM,ISPIN)*CLMSP(J,1,JATOM,ISPIN)
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             FXCK(LM1,ISPIN) = (VLIN1 + ACORR)/ELNUMSP(JATOM,ISPIN)

             WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), ESPXC(ISPIN)
             WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), DEXC(LM1,ISPIN), DVXC(LM1,ISPIN)

           ELSEIF ((L .NE. 9) .AND. (M .EQ. 4)) THEN

             WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), ESPXC(ISPIN)
             WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), DEXC(LM1,ISPIN), DVXC(LM1,ISPIN)

             GOTO 41

           ELSEIF ((M .EQ. 6) .OR. (M .EQ. 8) .OR. (M .EQ. 10)) THEN

             WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), ESPXC(ISPIN)
             WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), DEXC(LM1,ISPIN), DVXC(LM1,ISPIN)

             GOTO 41

           ELSEIF ((L .EQ. 4) .OR. (L .EQ. 6) .OR. (L .EQ. 7) .OR. (L .EQ. 9)) THEN

             DO J=1, JRI(JATOM)
               VALUE(J) = C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN)* &
                 (EXC(J,LM1,JATOM,ISPIN) - VXC(J,LM1,JATOM,ISPIN)) + &
                 C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN)* &
                 (EXC(J,LM1+1,JATOM,ISPIN) - VXC(J,LM1+1,JATOM,ISPIN))
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             ESPXC(ISPIN) = ESPXC(ISPIN) + VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN)*EXC(J,LM1,JATOM,ISPIN) + &
                 C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN)*EXC(J,LM1+1,JATOM,ISPIN)
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             DEXC(LM1,ISPIN) = VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN)*VXC(J,LM1,JATOM,ISPIN) + &
                 C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN)*VXC(J,LM1+1,JATOM,ISPIN)
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             DVXC(LM1,ISPIN) = VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = (C_KUB(L,M)*VXC(J,LM1,JATOM,ISPIN) + &
                 C_KUB(L,M+4)*VXC(J,LM1+1,JATOM,ISPIN))*CLMSP(J,1,JATOM,ISPIN)
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR=VALUE(1)*R(1)*0.5D0
             FXCK(LM1,ISPIN) = (VLIN1 + ACORR)/ELNUMSP(JATOM,ISPIN)

             WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), ESPXC(ISPIN)
             WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), DEXC(LM1,ISPIN), DVXC(LM1,ISPIN)

           ELSEIF ((L .EQ. 8) .OR. (L .EQ. 10)) THEN

             DO J=1, JRI(JATOM)
               VALUE(J) = C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN)* &
                 (EXC(J,LM1,JATOM,ISPIN) - VXC(J,LM1,JATOM,ISPIN)) + &
                 C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN)* &
                 (EXC(J,LM1+1,JATOM,ISPIN) - VXC(J,LM1+1,JATOM,ISPIN)) + &
                 C_KUB(L,M+8)*CLMSP(J,LM1+2,JATOM,ISPIN)* &
                 (EXC(J,LM1+2,JATOM,ISPIN) - VXC(J,LM1+2,JATOM,ISPIN))
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             ESPXC(ISPIN) = ESPXC(ISPIN) + VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN)*EXC(J,LM1,JATOM,ISPIN) + &
                 C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN)*EXC(J,LM1+1,JATOM,ISPIN) + &
                 C_KUB(L,M+8)*CLMSP(J,LM1+2,JATOM,ISPIN)*EXC(J,LM1+2,JATOM,ISPIN)
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             DEXC(LM1,ISPIN) = VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN)*VXC(J,LM1,JATOM,ISPIN) + &
                 C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN)*VXC(J,LM1+1,JATOM,ISPIN) + &
                 C_KUB(L,M+8)*CLMSP(J,LM1+2,JATOM,ISPIN)*VXC(J,LM1+2,JATOM,ISPIN)
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             DVXC(LM1,ISPIN) = VLIN1 + ACORR

             DO J=1, JRI(JATOM)
               VALUE(J) = (C_KUB(L,M)*VXC(J,LM1,JATOM,ISPIN) + &
                 C_KUB(L,M+4)*VXC(J,LM1+1,JATOM,ISPIN) + &
                 C_KUB(L,M+8)*VXC(J,LM1+2,JATOM,ISPIN))*CLMSP(J,1,JATOM,ISPIN)
             ENDDO
             CALL CHARG2(R,DX(JATOM),VALUE,1,JRJ,VLIN1)
             ACORR = VALUE(1)*R(1)*0.5D0
             FXCK(LM1,ISPIN) = (VLIN1 + ACORR)/ELNUMSP(JATOM,ISPIN)

             WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), ESPXC(ISPIN)
             WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), &
               FXCK(LM1,ISPIN), DEXC(LM1,ISPIN), DVXC(LM1,ISPIN)

           ENDIF

 41      CONTINUE

         WRITE(6,655) JATOM, ESPXC(ISPIN)
         WRITE(49+ISPIN,654) JATOM, ESPXC(ISPIN)

       ENDIF    !non-cubic or cubic

      endif
     endif     !  EECE end of potential and energy output
         DO J=1,JRI(JATOM)
            R(J)=R0(JATOM)*EXP(DX(JATOM)*(J-1))
            VTOT(J,1,JATOM)=VTOT(J,1,JATOM)*R(J)/SQFP
         ENDDO
!
      if(switch1.ne.' EECE')then   !EECE
         WRITE(ITAPS,2011) 0,0
         WRITE(ITAPS,2021) ( VTOT(J,1,JATOM), J=1,JRI(JATOM) )
         WRITE(ITAPS,2031)
         WRITE(ITAPS,2030)
!
         DO LM1=2,LMMAX(JATOM)
            WRITE(ITAPN,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
            WRITE(ITAPN,2021) ( VTOT(J,LM1,JATOM), J=1,JRI(JATOM) )
            WRITE(ITAPN,2031)
         ENDDO
!
         WRITE(ITAPN,2030)
       endif                 !EECE
      ENDDO
!
!      CLOSE(itapn)

     if(switch1.ne.' EECE')then   !EECE
      IF(myid.EQ.icoul) THEN
!      IF(myid.EQ.0) THEN
!         open(itapn,file=fnavn(itapn),POSITION='append')
         WRITE(itapn,*)    '    TOTAL POTENTIAL IN INTERSTITIAL'
         WRITE(itapn,*)
         WRITE(itapn,2061) nkk
         if(inversion)then
         DO j=1,nkk
            WRITE(itapn,2071) (kzz(jx,j),jx=1,3), DBLE(cvalue(j)),0.D0
         ENDDO
         else
         DO j=1,nkk
            WRITE(itapn,2071) (kzz(jx,j),jx=1,3), CONJG(cvalue(j))
         ENDDO
         endif
!      CLOSE(itapn)
      ENDIF
      CLOSE(itaps)
      CLOSE(itapn)
     endif                 !EECE
#ifdef Parallel
      IF(myid.LT.icoul) THEN
         CALL MPI_SEND(i,1,MPI_INTEGER,myid+1, 0,comm_coul,ierr)
      ENDIF
#endif
  ENDDO SPIN
!
!      WRITE(6,702)IFF1,IFF2,IFF3
!
!.....calculate used cpu-time
!
!      TIME6=SECOND()
  CALL CPUTIM(time6)
  timtot=time6-time0
  tmulti=time1-time0
  pmulti=tmulti/timtot*100.
  tcouli=timey-time1
  pcouli=tcouli/timtot*100.
  tcoul=time2-time1-tcouli-t947
  pcoul=tcoul/timtot*100.
  txcsph=time3-time2
  pxcsph=txcsph/timtot*100.
     if(switch1.ne.' EECE')then          !EECE
  txcpw=time4-time3
  pxcpw=txcpw/timtot*100.
  tenerg=time5-time4
  penerg=tenerg/timtot*100.
  treana=time6-time5
  preana=treana/timtot*100.
     endif
  WRITE(6,3000) timtot,100.
  WRITE(6,3010) tmulti,pmulti
  WRITE(6,3021) tcouli,pcouli
  WRITE(6,3022) t947,t947/timtot*100.
  WRITE(6,3020) tcoul,pcoul
  WRITE(6,3030) txcsph,pxcsph
     if(switch1.ne.' EECE')then          !EECE
  WRITE(6,3040) txcpw,pxcpw
  WRITE(6,3050) tenerg,penerg
  WRITE(6,*) '            TIME rean0, rean3',time4a-time4,time4b-time4a
  WRITE(6,3060) treana,preana
     endif
        IF(myid.EQ.0)  CALL errclr(errfn)
!
889 CONTINUE
!  write(6,*) 'Stay cool, honey bunny'
#ifdef Parallel
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
#endif
  STOP ' LAPW0 END'
!
!        error handling
!
910 INFO = 1
!
!        lapw0.def couldnt be opened
!
  WRITE (ERRMSG,9000) FNAME
  CALL OUTERR('LAPW0',ERRMSG)
  GOTO 999
920 INFO = 2
!
!        file FNAME couldn''t be opened
!
  WRITE (ERRMSG,9010) IUNIT
  CALL OUTERR('LAPW0',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  CALL OUTERR('LAPW0',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  CALL OUTERR('LAPW0',ERRMSG)
  GOTO 999
930 INFO = 3
!
!        illegal number of equivalent atoms
!
  CALL OUTERR('LAPW0','MULT .EQ. 0')
  GOTO 999
960 INFO = 7
!
!        Error reading file 'lapw2.def'
!
  WRITE (ERRMSG,9040) FNAME
  CALL OUTERR('LAPW0',ERRMSG)
  GOTO 999
970 INFO = 8
  WRITE (ERRMSG,9050) JATOM,LMMAX1(JATOM),NCOM+3
  CALL OUTERR('LAPW0',ERRMSG)
  GOTO 999
971 INFO = 8
  WRITE (ERRMSG,9051) JATOM, l,lmax2
  CALL OUTERR('LAPW0',ERRMSG)
  GOTO 999
901 INFO = 11
  CALL OUTERR('LAPW0','IBES.GT.20')
  GOTO 999
999 STOP 'LAPW0 - Error. Check file lapw0.error.'
!
!
!
77 FORMAT(1X,9(F9.7))
78 FORMAT(1X,A20,A4,3F10.6,I5)
699 FORMAT(':FIT',i3.3,':',1x,'SIGMA OF V-XC FIT FOR ATOM', &
         i3,e15.7)
700 FORMAT(I3,A77)
702 FORMAT(' IFFT PARAMETERS:',3I3)
705 FORMAT(/,/,12X,9('-'),/,':ITE',i3.3,':',I3,'. ITERATION',/,12X,9('-'),/)
710 FORMAT (7x,'HELLMAN-FEYMAN-FORCE IN mRy/a.u. = |F|',2x,'Fx',7x,'Fy',13x,'Fz')
711 FORMAT (':FHF',i3.3,':',1x,i3,'.ATOM',4f15.3,/)
787 FORMAT('     TOTAL SPHERICAL POTENTIAL IN MT SPHERES',5X,I3,'.ITERATION')
800 FORMAT(/,/,/,/,30X,50(1H-),/,33X,'S T R U C T U R A L   ','I N F O R M A T I O N',/,30X,50(1H-),/,/)
805 FORMAT(7X,'SUBSTANCE: ',A60,/)
810 FORMAT(7X,'LATTICE',22X,'= ',A4)
819 FORMAT(':POT  : POTENTIAL OPTION',i5)
820 FORMAT(':LAT  :',1x,'LATTICE CONSTANTS=',3F9.5,3F9.3)
821 FORMAT(':VOL  :',1x,'UNIT CELL VOLUME =',F14.5)
840 FORMAT(7X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)
850 FORMAT(7X,'MODE OF CALCULATION IS',7X,'= ',A4)
860 FORMAT(7X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)
870 FORMAT(1H1,/,30X,30(1H-),/,30X,'M U L T I P O L M O M E N T S',/,30X,30(1H-))
880 FORMAT(1H1,/,30X,40(1H-),/,33X,'C O N V E R G E N C E   PARA','METERS',/,30X,40(1H-),/,/)
991 FORMAT(2E15.7,3I5)
994 FORMAT(3X,'L=',I3,3X,'M=',I2,3X,'PSEUDO MULTIPOLMOMENT     = ',2F13.9)
995 FORMAT(A4,I4)
996 FORMAT(3X,'L=',I3,3X,'M=',I2,3X,'PLANE WAVE MULTIPOLMOMENT = ',2F13.9)
997 FORMAT(A4,2A5)
1011 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1022 FORMAT(5E14.7)
1121 FORMAT(':DEN  :',' DENSITY INTEGRALS =',3x,F20.6,3X,'(Ry)')           
1333 FORMAT(2I5,A10)
1334 FORMAT(6E12.5)
1335 FORMAT(7X,'ATOMNUMBER=',I3,1X,A10,1X,'VCOUL-ZERO = ',E12.5)
1340 FORMAT(':EFG',i3.3,':',24X,'EFG         =',F12.5,'   *10**21  V / m**2')         
1341 FORMAT(31X,'V20  TOT/SRF=',2F12.5)                                 
1342 FORMAT(31X,'V22  TOT/SRF=',2F12.5)                               
1343 FORMAT(31X,'V22M TOT/SRF=',2F12.5)                               
1344 FORMAT(31X,'V21  TOT/SRF=',2F12.5)                               
1345 FORMAT(31X,'V21M TOT/SRF=',2F12.5,/)                               
1350 FORMAT(10X,3F11.5,4X,3F11.5)                                            
1355 FORMAT(/,9X,'MAIN DIRECTIONS OF THE EFG ',3F8.4,/,36X,3F8.4,/,36X,3F8.4,/,&
 & ':ANG',i3.3,':',2X,'ANGLE WITH OLD X-AXIS =    ',f8.1)
1360 FORMAT(/,':ETA',i3.3,':',25X,'ASYMM. ETA =',F12.5,/,/) 
1520 FORMAT(/,/,3X,'CONVERGENCE PARAMETER FOR PSEUDOCHARGE: NCON= ', &
          I2,/,3X,'MAXIMAL VALUE OF RMT(JATOM)*ABSK(NKK) : RK  =', &
          F8.5,/,/)
1540 FORMAT(3X,'MATRIX FACTORISATION  OF LEAST-SQUARE-FITS IN ', &
          'XCPOT1,XCPOT3')
1550 FORMAT(/,7X,'VALUE OF V0-EXCHANGE',11X,'V0 = ',2E13.5, &
          / ,7X,'VALUE OF V -BAR     ',8X,'V-BAR = ',2E13.5)
1560 FORMAT(7X,'VALUE OF MUFFIN-TIN ZERO',10X,'= ',2E13.5)
1040 FORMAT(8I2)
1060 FORMAT(/,/,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')
1070 FORMAT(30X,3F10.5)
1080 FORMAT(13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)
1970 FORMAT(3X,'             exchange correlation potential ', &
          'BY',I3,'. ITERATION (NORM: CLM=CLM*R*R)')
1980 FORMAT(3X)
1990 FORMAT(3X,'ATOMNUMBER =',I3,5X,10A4)
2000 FORMAT(15X,I3,/,/)
2001 FORMAT(3X,'NUMBER OF LM',I3,/,/)
2010 FORMAT(15X,I3,5X,I2,/)
2011 FORMAT(3X,'VLM(R) FOR L',I3,3X,'M=',I2,/)
2020 FORMAT(3X,5E13.7)
2021 FORMAT(3X,4E19.12)
2030 FORMAT(/,/,/)
2031 FORMAT(/)
2040 FORMAT(/,/,3X,'ATOM=',I3,3X,A10,3X,'Z=',F5.2,3X, &
          'LM=',I3,5X,'POSITION=',3F7.3,/)
2044 FORMAT(49X,I3,/,/,/)
2050 FORMAT(3X,'L=',I3,3X,'M=',I2,3X, &
          'SPHERE MM = ',2F13.9,3x,2f6.2)
2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')
2060 FORMAT(/,13X,I6)
2070 FORMAT(3X,3I5,2E15.7)
2071 FORMAT(3X,3I5,2E19.12)
3000 FORMAT(/,/,3X,'=====>>> CPU-TIME SUMMARY',/, &
          12X,'TOTAL CPU-TIME USED : ',F7.1,5X,F4.0,' PERCENT')
3010 FORMAT(12X,'TIME MULTIPOLMOMENTS: ',F7.1,5X,F4.0,' PERCENT')
3020 FORMAT(12X,'TIME COULOMB POT SPH: ',F7.1,5X,F4.0,' PERCENT')
3021 FORMAT(12X,'TIME COULOMB POT INT: ',F7.1,5X,F4.0,' PERCENT')
3022 FORMAT(12X,'TIME COULOMB POT RMT: ',F7.1,5X,F4.0,' PERCENT')
3030 FORMAT(12X,'TIME XCPOT SPHERES  : ',F7.1,5X,F4.0,' PERCENT')
3040 FORMAT(12X,'TIME XCPOT INTERST  : ',F7.1,5X,F4.0,' PERCENT')
3050 FORMAT(12X,'TIME TOTAL ENERGY   : ',F7.1,5X,F4.0,' PERCENT')
3060 FORMAT(12X,'TIME REANALYSE      : ',F7.1,5X,F4.0,' PERCENT')
9000 FORMAT('can''t open definition file ',A40)
9010 FORMAT('can''t open unit: ',I2)
9020 FORMAT('       filename: ',A50)
9030 FORMAT('         status: ',A,'  form: ',A)
9040 FORMAT('Error reading file: ',A47)
9050 FORMAT('lmmax1 for jatom',i3,' of ',i4,' gt ',i3)
9051 FORMAT('l-value for jatom',i3,' =',i3,' gt lmax2=',i2)
9060 FORMAT('NWAV TOO SMALL COMPARED WITH NKK',2i7)
9070 FORMAT('NEEDED VALUES:',3I4)
9080 FORMAT('PRESENT PARAMETERS:',3I4)
!
!        End of 'LAPW0'
!
END
