!cad
!cad         J O I N T  D E N S I T Y  O F  S T A T E S
!cad
!cad            ACTUAL VERSION  with BLOECHL SCHEME
!cad
!cad         written by Robert Abt startind from TETRA
!cad         modifications by cad, November 1998  
!cad         modifications by Jan Kunes, May 1999
!cad         modifications by cad, May-August 1999 
!cad         modifications by cad, August 2002 
!cad
!cad
!cad            FILE  3       case.outmat       MOMENTUM MATRIX ELEMENTS
!cad            FILE  4       case.weight       ENERGY BANDS, WHEIGHTS
!cad            FILE  5       case.injoint      INPUT
!cad            FILE  6       case.outputjoint  OUTPUT
!cad            FILE  7       case.joint        DOS, JDOS IM(EPSILON) 
!cad            FILE  8       case.sigma_intra  intraband contributions
!cad            FILE 14       case.kgen         TETRAHEDRA
!cad            FILE 20       case.struct       STRUCTURAL DATA
!cad
!cad            for band analysis further files are usd
!cad
      PROGRAM JOINT
      use felder
      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'param.inc'
      PARAMETER (NPF=10)
      PARAMETER (RYDeV=13.605698)
      PARAMETER (E    =1.602E-19)
      PARAMETER (H    =6.625E-34)
      PARAMETER (C    =2.99793E8)
!ad
      CHARACTER*2  aif,HINTRA(MG0)
      CHARACTER*4  ECV
      CHARACTER*6  ECV1
      CHARACTER*6  HSIGMA,HELOSS
      CHARACTER*7  HIMEPS,HREEPS
      CHARACTER*9  HEADER(MG0),HHEAD(9)
      CHARACTER*11 FORM,STATUS
      CHARACTER*67 ERRMSG
!ad   CHARACTER*70 SYSTEM
      CHARACTER*80 DEFFN, ERRFN, FNAME
      CHARACTER*80 ffile,f1
      CHARACTER*117 HEAD
!ad
      INTEGER hh
!ad
      REAL*8 imeps(MG0),reeps(MG0),sigma(MG0),eloss(MG0),sumr(MG0)
      REAL*8 gamma(MG0),plasm2(MG0),sig1,eps1,eps2
!      REAL*4 OPMAT
!ad
      LOGICAL SO,SPIN
!ad
!      DIMENSION  ENG(MET),ENG2(MET),gesDENSTY(MET,MG)
      real*8, allocatable ::  ENG(:),ENG2(:),gesDENSTY(:,:) ! met,met,(MET,MG)
      DIMENSION  npcol(MG0)
      integer, allocatable :: i1(:),i2(:) !INUMEden
!ad
      COMMON /NCOUNT/ NNOC
!      COMMON /SN/ DENSTY(INUMEden,MET,MG)
!      COMMON /BST/ EBS(NKPT,NUME), FC(NKPT,NUME)
      COMMON /SWITCH/ ISWITCH
      COMMON /EME/ EEF,EMIN, EMAX, EFACTR, DE, NFIRST, NCOL, NLAST
!
      COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
!
!      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
!                      PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO), &
!                      JRI(NATO),R0(NATO)
!      COMMON /CHAR/   TITLE,LATTIC,MODUS,ANAME(NATO)
!
!      COMMON /OPME/ OPMAT(NKPT,INUME,MG), &
!                    EMINo,EMAXo,OML,OM1,MIMA(NKPT,2),NK,KRA
      COMMON /OPME/ EMINo,EMAXo,OML,OM1,NK,KRA
                    
      DATA PI /3.141592654/,  ONE/1.0D+0/
!ad
!ad
!ad ________________________ INITIALIZE VARIABLES ____________________
!ad
!ad
      SPIN=.FALSE.
      SO=.FALSE.
      ISWITCH=0
!ad
!ad
!ad ____________________________ OPEN FILES ____________________________
!ad
      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in JOINT')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
   
      if(iunit.eq.3) then
        f1=fname
        do i=80,1,-1
        if(fname(i:i).ne.' ') then
          if(fname(i-7:i).eq.'outmatup') SPIN=.true.
          if(fname(i-7:i).eq.'outmatdn') SPIN=.true.
          goto 10
        endif
        enddo
      endif
     
      GOTO 10
   20 CONTINUE
      CLOSE (1)

!---------------------------------------------------------------------
!
!     NST          number of k points
!     NYMIN        LOWER BAND INDEX
!     NYMAX        UPPER BAND INDEX
!fb
!     NYOCC        LAST OCCUPIED BAND INDEX
!     EMIN,DE,EMAX ENERGY GRID FOR CDOS
!
!---------------------------------------------------------------------
      CALL CPUTIM(ONTIME)
!ad
!ad ____________________________ READ INPUT __________________________
!ad
!ad   READ(5,1) SYSTEM 
!ad   if(spin) write(6,*) 'spinpolarized'
      nst=1
!fb      READ(5,*) NYMIN,NYMAX
      READ(5,*,err=101) NYMIN,NYMAX,NYOCC
      GOTO 102
 101  NYOCC=NYMAX-1
 102  CONTINUE     
      READ(5,*) EMIN,DE,EMAX
!      if(nymax.gt.NUME) then
!        write(*,*) ' NUME',NUME,' .lt. nymax=',nymax
!        stop ' nume .lt. nymax '
!      end if                           
      JEND=1 + (EMAX-EMIN)/DE
!      IF(JEND.GT.MET) STOP 'JEND GT. MET'
      met=jend
!      allocate (ENG(jend),ENG2(jend),gesDENSTY(jend,mg)) 
      gesdensty=0.d0
      READ(5,11) ECV
      EC=ONE
      IF (ECV(1:2).EQ.'eV'.or.ECV(2:3).eq.'eV'.or.ECV(3:4).eq.'eV' &
        .or.ECV(1:2).EQ.'EV'.or.ECV(2:3).eq.'EV'.or.ECV(3:4).eq.'EV') &
        then 
        EC=RYDeV
        ecv1='  [eV]' 
      ELSEIF (ECV.EQ.'cm-1') then 
        EC=RYDeV*E/H/C/100.   
        ecv1='[cm-1]'
      ELSE
        ecv1='  [Ry]'
!ad
        EC=1.0d0
!ad 
      ENDIF
      READ(5,*) ISWITCH
      READ(5,*) NCOL
     
      if (iswitch.lt.4) NCOL=1
      mg=ncol
      allocate (ENG(jend),ENG2(jend),gesDENSTY(jend,mg)) 
      if(iswitch.eq.6.or.iswitch.eq.7) READ(5,*) (gamma(i),i=1,NCOL)
!ad
!ad          0...JOINTDOS FOR EACH BAND COMBINATION  '
!ad          1...JOINTDOS AS SUM OVER ALL BAND COMBINATIONS'
!ad          2...DOS FOR EACH BAND '
!ad          3...DOS AS SUM OVER ALL BANDS'
!ad          4...Im(EPSILON) '
!ad          5...Im(EPSILON) for each band combination'
!ad          6...INTRABAND contributions'
!ad          7...INTRABAND contributions including band analysis'
!ad
!ad
!ad _____________________ READ STRUCTURAL INFORMATION __________________
!ad
      CALL RSTRU
!ad
!ad ___________________________ DEFINE HEADERS _________________________
!ad
!ad
      IF (ISWITCH.GT.3) then
      hhead(1)='Re <x><x>'
      hhead(2)='Re <y><y>'
      hhead(3)='Re <z><z>'
      hhead(4)='Re <x><y>'
      hhead(5)='Re <x><z>'
      hhead(6)='Re <y><z>'
      hhead(7)='Im <x><y>'
      hhead(8)='Im <x><z>'
      hhead(9)='Im <y><z>'
      himeps='Im(eps)'
      hreeps='Re(eps)'
      hsigma='sigma_'
      heloss='eloss_'
!cad
      READ(3,330) ncol1,head
      if(ncol.gt.ncol1) ncol=ncol1
         lcol=0
      if(iswitch.eq.6.or.iswitch.eq.7) then
      READ(9,330) ncol1,head
         nplcol=ncol
         lcol=0
      endif

      do 222 icol=1,ncol
        ih1=3+13*(icol-1)
!cad
       if(head(ih1+4:ih1+7).eq.'x><x') then
         header(icol)='Im_eps_xx'
         lcol=lcol+1
         hintra(lcol)='xx'
         npcol(lcol)=icol
         goto 222
       endif
!ad
       if(head(ih1+4:ih1+7).eq.'y><y') then 
         header(icol)='Im_eps_yy'
         lcol=lcol+1
         hintra(lcol)='yy'
         npcol(lcol)=icol
         goto 222
       endif
!ad
      if(head(ih1+4:ih1+7).eq.'z><z') then 
         header(icol)='Im_eps_zz'
         lcol=lcol+1
         hintra(lcol)='zz'
         npcol(lcol)=icol
         goto 222
       endif
!ad
       if((head(ih1:ih1+1).eq.'Re').and.(head(ih1+4:ih1+7).eq.'x><y'))  &
       then 
         header(icol)='Im_eps_xy'
         if(iswitch.eq.6.or.iswitch.eq.7) then
         write(6,444) icol
         endif
         goto 222
       endif
!ad
       if((head(ih1:ih1+1).eq.'Re').and.(head(ih1+4:ih1+7).eq.'x><z')) &
       then 
         header(icol)='Im_eps_xz'
         if(iswitch.eq.6.or.iswitch.eq.7) then
         write(6,444) icol
         endif
         goto 222
       endif
!ad
       if((head(ih1:ih1+1).eq.'Re').and.(head(ih1+4:ih1+7).eq.'y><z')) &
       then 
         header(icol)='Im_eps_yz'
         if(iswitch.eq.6.or.iswitch.eq.7) then
         write(6,444) icol
         endif
         goto 222
       endif
!ad
       if((head(ih1:ih1+1).eq.'Im').and.(head(ih1+4:ih1+7).eq.'x><y')) &
       then 
         header(icol)='Re_eps_xy'
         if(iswitch.eq.6.or.iswitch.eq.7) then
         write(6,444) icol
         endif
         goto 222
       endif
!ad
       if((head(ih1:ih1+1).eq.'Im').and.(head(ih1+4:ih1+7).eq.'x><z')) &
       then 
         header(icol)='Re_eps_xz'
         if(iswitch.eq.6.or.iswitch.eq.7) then
         write(6,444) icol
         endif
         goto 222
       endif
!ad
       if((head(ih1:ih1+1).eq.'Im').and.(head(ih1+4:ih1+7).eq.'y><z')) &
        then 
         header(icol)='Re_eps_yz'
         if(iswitch.eq.6.or.iswitch.eq.7) then
         write(6,444) icol
         endif
        endif
!ad
 222    continue

!       opmat=0.0
!ad
!ad ____________________ READ MOMENTUM MATRIX ELEMENTS _________________
!ad
!ad
      if(iswitch.eq.6.or.iswitch.eq.7) then
      call read_diag(ncol)
      else
      CALL READOPMAT(NCOL,NYOCC)
      endif
      END IF
!ad
!ad
!ad ____________ READ EIGENVALUES AND WEIGHTS FOR ALL K-POINTS _________
!ad
!ad
!ad            calculate sum of weights for lowest energy level
!ad
      sumw=0.0
!ad
      READ(4,*)
      READ(4,809) EEF, NST 
!      if(nst.gt.NKPT) then
!        write(6,79) NST,NKPT
!        stop 'nst. gt. nkpt'
!      end if
!     write(*,*)'EF: ',eef,'K points: ',nst
!  
!....initialize bandenergies and weights.............
!     
!      IF (NYMAX.GT.NUME) STOP ' NUME .lt. NEMAX '
!ad
      allocate (EBS(nst,nymax), FC(nst,nymax))     ! NKPT,NUME
      DO K=1,NST
       DO I=1,NYMAX
        EBS(K,I)=0.0
        FC(K,I) =0.0
       ENDDO
      ENDDO
!
!...NYMAXi controls NYMAX ..........
!
      NYMAXi=0
      DO 190 K=1,NST
      READ(4,789) KK,NEK 
!fb
      READ(4,799) EBS(K,1),fci
      if (fci.gt.1.0d-10) then
        WKONE=fci
        sumw=sumw + WKONE
      else
        WKONE=1.0d0
      endif	 
      FC(K,1)=fci/WKONE
      DO 190 I=2,NEK
!      IF (NUME.GE.I)  THEN
      IF (Nymax.GE.I)  THEN
         READ(4,799) EBS(K,I),fci
         FC(K,I)=fci/WKONE
      ELSE
         READ(4,799) EBSDUMMY,FCDUMMY
      END IF
!ad
      if(SO) then
        dummy=real(nymin)/2.d0-nymin/2
        if(abs(dummy).lt.1.d-5) then
          NYMIN=NYMIN-1
          write(*,*)
          write(*,*) '  NYMIN decreased by one!' 
        endif
      endif
!ad

      NYMAXi=max(NYMAXi,NEK)
 190  CONTINUE
      NYMAX=min(NYMAX,NYMAXi)
!fb
      NYOCC=min(NYOCC,NYMAX-1)      
!     write(*,*) ' after weightfile : (nemin,nemax) ',NYMIN,NYMAX
      IF (NYMIN.GE.NYMAX) THEN
        write(*,*) ' ERROR: no proper joice of band indices! '
        write(*,*) ' NEMIN,NEMAX: ',NYMIN,NYMAX
        stop 'see band indices ! '
      ENDIF
!fb
!ad
!ad   WRITE(7,3) SYSTEM
!fb      WRITE(6,5) NYMIN,NYMAX,EMIN,DE,EMAX
      
      WRITE(6,5) NYMIN,NYMAX,NYOCC,EMIN,DE,EMAX
!ad
!ad
!      IF (ISWITCH.ne.5.AND.ISWITCH.ne.0) THEN
!        HH=NYMAX-NYMIN+1
!      ELSE
!        HH=NYMAX-NYMIN+1
!        HH=(HH*(HH-1))/2
!      ENDIF
!fb
      IF (ISWITCH.eq.5.or.ISWITCH.eq.0) THEN
!fb        HH=NYMAX-NYMIN+1
!fb        HH=(HH*(HH-1))/2
!          HH=(NYMAX+1-(NYMIN+NYOCC)/2)*(NYOCC-NYMIN+1)
        HH=(NYMAX+1)*(NYOCC-NYMIN+1)-(NYMIN+NYOCC)*(NYOCC-NYMIN+1)/2
      ELSE IF (ISWITCH.eq.1.or.ISWITCH.eq.4.or.ISWITCH.eq.6) THEN
        HH=1
      ELSE
        HH=NYMAX-NYMIN+1
      ENDIF
!fb
!ad
!ad
      allocate (DENSTY(hh,jend,MG),i1(hh),i2(hh))
          densty=0.d0
!      if(iswitch.eq.6.or.iswitch.eq.7) then
!        IF (HH.GT.INUMEden) then
!           WRITE(6,77) HH,INUMEden
!           STOP 'INUMEDEN TOO SMALL'
!        ENDIF
!      endif
      if(iswitch.eq.0.or.iswitch.eq.2.or.iswitch.eq.5) then
!        IF (HH.GT.INUMEden) THEN
!           WRITE(6,77) HH,INUMEden
!           if(iswitch.eq.0) iswitch=1
!           if(iswitch.eq.2) iswitch=3
!           if(iswitch.eq.5) iswitch=4
!           goto 777
!        ENDIF
           nfiles=hh/npf
           nrest=hh-npf*nfiles
           if(nrest.gt.0) nfiles=nfiles+1
           write(6,78) hh,nfiles
           call filename(f1,length,ffile)
        endif
 777  continue

!ad
!ad ________________ ENERGY MESH FOR PLASMA FREQUENCY ____________________
!ad
!ad
      IF(ISWITCH.EQ.6.OR.ISWITCH.EQ.7) THEN
      npl=10
      EMINpl=EMIN
      EMAXpl=EMAX
      EMIN=EEF-npl*DE
      EMAX=EEF+npl*DE
      JEND=1 + (EMAX-EMIN)/DE
!     write(*,*) jend,met,emin,emax,de
      IF(JEND.GT.MET) STOP 'JEND GT. MET'
      ENDIF
!...
      NFIRST=1
      NLAST=JEND
      write(*,*) ' SUM OF WEIGHTS: ',sumw
!ad
!ad
!ad ____________________________ PREFACTORS ____________________________
!ad

      EFACTR=1.0D+0/DE
!ad   ESF =sumw / 2.d0 / EC
      ESF =sumw        / EC
      EPLF=sumw * 2.d0
      EJDF=sumw / 2.d0 / EC * 64*PI**2/ VOL
      EPSF=sumw / 2.d0 * EC**2 * 64*PI**2/ VOL
!ad   if(SO) then
!ad   EJDF=EJDF/2.d0
!ad   if(.not.SPIN) EPSF=EPSF*2.d0
!ad   if(.not.SPIN) EPLF=EPLF*2.d0
!ad   endif
!ad
!ad   PLASMA FREQUENCY:
!ad
!.....(hbar w(pl) )^2 [ryd^2] = 16 Pi / V(a.u.) * n (electrons per cell)
      plasm2fac= 16.d0 * Pi / VOL * EC**2
      sumfac=8.d0 * Pi * Pi / VOL * EC
!ad
!.....sig= i w / 4Pi * ( 1 - eps ) [cgs] but here in 1/(ohm cm).........
!
      sigfac = 134.59 * RYDeV
!ad
!ad _______________________ GENERATE ENERGY MESH _______________________
!ad
                 DO J=1,JEND
                 ENG(J)=EMIN+(J-1)*DE
                 ENG(J)=ENG(J)*EC
                 ENG2(j)=ENG(j)*ENG(j)
                 enddo
!ad
!ad ____________________________________________________________________
!ad
!fb      CALL ARBDOS(NYMIN,NYMAX)
      CALL ARBDOS(NYMIN,NYMAX,NYOCC)
!
         CALL CPUTIM(DETIME)
         DETIME=DETIME-ONTIME
         WRITE(6,13) DETIME  
!ad
!ad ________________________ WRITE OUT HEADERS _________________________
!ad
      if(iswitch.eq.4.or.iswitch.eq.5) then
         write(7,110) ncol,VOL
         write(6,111) ecv1,(header(icol),icol=1,ncol)
         write(7,111) ecv1,(header(icol),icol=1,ncol)
         write(6,*)
         write(7,114)
      endif
      if(iswitch.eq.0.or.iswitch.eq.1) then
         write(6,211) ecv1
         write(7,211) ecv1
      endif
      if(iswitch.eq.2.or.iswitch.eq.3) then
         write(6,311) ecv1
         write(7,311) ecv1
      endif
!ad
      if(iswitch.eq.0.or.iswitch.eq.1) goto 100
      if(iswitch.eq.2.or.iswitch.eq.3) goto 200
      if(iswitch.eq.4.or.iswitch.eq.5) goto 400
      if(iswitch.eq.6.or.iswitch.eq.7) goto 600
!ad
!ad ____________________________ case JDOS _____________________________
!ad

 100   CONTINUE

       NF=1
!ad
      IF (ISWITCH.EQ.1) THEN
        IB=1
        DO 151 J=2,JEND
        WRITE(7,116) ENG(j), &
             EJDF*DENSTY(IB,J,1),EJDF*DENSTY(IB,J,1)/ENG2(j)
 151    WRITE(6,116) ENG(j), &
             EJDF*DENSTY(IB,J,1),EJDF*DENSTY(IB,J,1)/ENG2(j)
      ENDIF
!ad
      IF (ISWITCH.EQ.0) THEN
       IB=0
!fb         DO 141 II=NYMIN,NYMAX-1
         DO 141 II=NYMIN,NYOCC
         DO 141 JJ=II+1,NYMAX
        IB=IB+1
          WRITE(6,130) II,JJ
          i1(ib)=ii
          i2(ib)=jj
!ad
        DO 141 J=2,JEND
!ad
!ad............................BAND ANALYSIS............................
!ad
        WRITE(6,116)  &
        ENG(j),EJDF*DENSTY(IB,J,NF),EJDF*DENSTY(IB,J,NF)/eng2(j)
!ad
!ad.......................written to output file........................
!ad
        gesDENSTY(J,NF)=gesDENSTY(J,NF)+EJDF*DENSTY(IB,J,NF)
 141    CONTINUE
 142    continue
!ad
        write(6,132)
!ad
        DO 251 J=2,JEND
        WRITE(7,116) ENG(j), &
             (gesDENSTY(J,Nd),gesDENSTY(J,Nd)/ENG2(j),Nd=1,NCOL)
 251    WRITE(6,116) ENG(j), &
             (gesDENSTY(J,Nd),gesDENSTY(J,Nd)/ENG2(j),Nd=1,NCOL)
!ad
!ad
!ad............................BAND ANALYSIS............................
!ad
         do if=1,nfiles
         ifile=40+if
         ffile(length+1:length+6)='.jdos_'
         call dumm(if,aif)
         if(if.lt.10) then
         ffile=ffile(1:length+6)//aif(2:2)
         else
         ffile=ffile(1:length+6)//aif
         endif
         open(ifile,file=ffile,status='unknown',form='formatted') 
         write(ifile,266)
         ib1=(if-1)*npf+1
         ib2=ib1+npf-1
         if(ib2.gt.hh) ib2=hh
         write(ifile,360) (i1(ib),i2(ib),ib=ib1,ib2)

          do J=2,JEND
          write(ifile,166) eng(j),(EJDF*DENSTY(IB,J,NF),ib=ib1,ib2)
          enddo
         enddo
!ad
!ad........................written to plot file.........................
!ad
       endif
       goto 900
!ad
!ad  __________________________ case JDOS end __________________________
!ad
!ad
!ad _____________________________ case DOS _____________________________
!ad

200      CONTINUE

         IB=0      
         DO 240 II=NYMIN,NYMAX
         IB=IB+1  
         IF (IB.GE.2.AND.ISWITCH.EQ.3) GOTO 243
         IF (ISWITCH.EQ.2) WRITE(6,131) II 
         DO 240 J=1,JEND
!ad
!ad............................BAND ANALYSIS............................
!ad
         IF (ISWITCH.EQ.2) WRITE(6,116) ENG(j)-eef*ec,ESF*DENSTY(IB,J,1)
!ad
!ad.......................written to output file........................
!ad
         gesDENSTY(J,1)=gesDENSTY(J,1)+DENSTY(IB,J,1)
 240     CONTINUE 
         nband=ib
         if(ib.ne.hh) write(*,*) 'nband, hh: ',nband,hh
!ad
!ad............................BAND ANALYSIS............................
!ad
         if(ISWITCH.eq.2) then
!ad
         nfiles=nband/npf
         nrest=nband-npf*nfiles
         if(nrest.gt.0) nfiles=nfiles+1
!ad
        do if=1,nfiles
         ifile=40+if
         ffile(length+1:length+5)='.dos_'
         call dumm(if,aif)
         if(if.lt.10) then
         ffile=ffile(1:length+5)//aif(2:2)
         else
         ffile=ffile(1:length+5)//aif
         endif
         open(ifile,file=ffile,status='unknown',form='formatted')
         write(ifile,267)
         ib1=(if-1)*npf+1
         ib2=ib1+npf-1
         if(ib2.gt.hh) ib2=hh
      if(ib2-ib1.eq.0) write(ifile,461) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.1) write(ifile,462) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.2) write(ifile,463) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.3) write(ifile,464) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.4) write(ifile,465) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.5) write(ifile,466) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.6) write(ifile,467) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.7) write(ifile,468) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.8) write(ifile,469) (ib,ib=ib1,ib2)
      if(ib2-ib1.eq.9) write(ifile,460) (ib,ib=ib1,ib2)
!ad
         DO J=1,JEND
         write(ifile,166) ENG(j)-eef*ec, (ESF*DENSTY(IB,J,1),ib=ib1,ib2)
         enddo
        enddo
         endif
!ad
!ad........................written to plot file.........................
!ad
 243     CONTINUE
!ad
        IF (ISWITCH.EQ.2) write(6,232)
!ad
        DO 351 J=2,JEND
        WRITE(7,116) ENG(j)-eef*ec,ESF*gesDENSTY(J,1)
 351    WRITE(6,116) ENG(j)-eef*ec,ESF*gesDENSTY(J,1)
!ad
      goto 900
!ad
!ad  ___________________________ case DOS end __________________________
!ad
!ad
!ad  ____________________ case INTRABAND CONTRIBUTIONS _________________
!ad

 600     CONTINUE

         IB=0     
         DO 246 II=NYMIN,NYMAX
         IB=IB+1
         DO 246 J=1,JEND
         do i=1,lcol
         ip=npcol(i)
         gesDENSTY(J,ip)=gesDENSTY(J,ip)+DENSTY(IB,J,ip)
         enddo
 246     CONTINUE
         IBAND=IB
!ad
           do i=1,LCOL
            ip=npcol(i)
            plasm2(ip)= plasm2fac*(EPLF*gesDENSTY(npl+1,ip))
           enddo
           write(6,416)
           if(lcol.eq.1) write(6,5161) hintra(1),ecv1
           if(lcol.eq.2) write(6,5162) (hintra(i),i=1,2),ecv1
           if(lcol.eq.3) write(6,5163) (hintra(i),i=1,3),ecv1
           WRITE(6,618) (gamma(npcol(i)),i=1,lcol)
           if(lcol.eq.1) write(6,6101) ecv1,hintra(1)
           if(lcol.eq.2) write(6,6102) ecv1,(hintra(i),i=1,2)
           if(lcol.eq.3) write(6,6103) ecv1,(hintra(i),i=1,3)
 
!          endif
!ad
      DO J=1,JEND
      WRITE(6,116) ENG(j)-eef*ec,(EPLF*gesDENSTY(J,npcol(ni)),ni=1,lcol)
      ENDDO
!ad
!ad............................BAND ANALYSIS............................
!ad
      IF (ISWITCH.EQ.7) then
      write(6,315)
      do ib=1,iband
      if(DENSTY(IB,npl+1,1).gt.1.d-10) &
         WRITE(6,115) IB,(EPLF*DENSTY(IB,npl+1,npcol(i)),i=1,lcol)
      enddo
      ENDIF
!ad
!ad.......................written to output file........................
!ad
           if (SPIN) then
             write(6,317)
           else
             write(6,316)
           endif
!ad
           if(lcol.eq.1) write(6,6161) hintra(1),ecv1
           if(lcol.eq.2) write(6,6162) (hintra(i),i=1,2),ecv1
           if(lcol.eq.3) write(6,6163) (hintra(i),i=1,3),ecv1
           WRITE(6,618) (SQRT(plasm2(npcol(i))),i=1,LCOL)
!ad
!............restore energyscale for drude spectra................
!ad
      EMIN=EMINpl
      EMAX=EMAXpl
        write(6,612) ecv1, &
        (himeps,hintra(i),i=1,lcol),(hreeps,hintra(i),i=1,lcol)
        write(7,612) ecv1, &
        (himeps,hintra(i),i=1,lcol),(hreeps,hintra(i),i=1,lcol)
        write(6,*)
        write(7,114)
        write(8,614) (hsigma,hintra(i),i=1,lcol),(heloss,hintra(i),i=1,lcol)
        write(8,615) ecv1
      JEND=1 + (EMAX-EMIN)/DE
      DO  J=2,JEND 
        ENG(j)=EMIN+(J-1)*DE
        eng(j)=eng(j)*ec
!ad
        do i=1,lcol
        ip=npcol(i)
        imeps(i)=eps2(plasm2(ip),gamma(ip),ENG(j))
        reeps(i)=eps1(plasm2(ip),gamma(ip),ENG(j))
        sigma(i)=sigfac*EC*sig1(plasm2(ip),gamma(ip),ENG(j)) 
        eloss(i)=imeps(i)/(imeps(i)**2+reeps(i)**2)
        enddo
        WRITE(6,619) ENG(j),(imeps(i),i=1,lcol),(reeps(i),i=1,lcol)
        WRITE(7,619) ENG(j),(imeps(i),i=1,lcol),(reeps(i),i=1,lcol)
        WRITE(8,619) ENG(j),(sigma(i),i=1,lcol),(eloss(i),i=1,lcol)
      ENDDO

         goto 900
!ad
!ad  __________________ case INTRABAND CONTRIBUTIONS end _______________
!ad
!ad
!ad  ______________________ case DIELECTRIC TENSOR _____________________
!ad

  400   CONTINUE
         IF (ISWITCH.EQ.4) THEN
          IB=1 
          DO J=2,JEND   
          DO  241 icol=1,ncol
            VV=EPSF*DENSTY(IB,J,icol)/ENG2(j) 
            DENSTY(IB,J,icol)=VV
            gesDENSTY(J,icol)=gesDENSTY(J,icol)+VV
  241     CONTINUE
          END DO
         ENDIF

         IF(ISWITCH.EQ.5) THEN
         IB=0      
!fb         DO 341 II=NYMIN,NYMAX-1  
         DO 341 II=NYMIN,NYOCC 
         DO  341 JJ=II+1,NYMAX
         IB=IB+1
          WRITE(6,1130) II,JJ
          i1(ib)=ii
          i2(ib)=jj
         VM=0 

         DO J=2,JEND   
         DO 1241 icol=1,ncol 
           VV=EPSF*DENSTY(IB,J,icol)/ENG2(j) 
           DENSTY(IB,J,icol)=VV
           gesDENSTY(J,icol)=gesDENSTY(J,icol)+VV
1241    CONTINUE
!ad
!ad............................BAND ANALYSIS............................
!ad
        WRITE(6,116) ENG(j),(DENSTY(IB,J,icol),icol=1,ncol)
!ad
!ad.......................written to output file........................
!ad
          END DO
 341     CONTINUE                
!ad
!ad............................BAND ANALYSIS............................
!ad
         if(nfiles.gt.20) write(6,333)
!ad
         do icol=1,ncol
         do if=1,nfiles
         ifile=40+if
         ffile(length+1:length+1)='.'
	 ffile(length+2:length+10)=header(icol)
	 ffile(length+11:length+11)='_'
         call dumm(if,aif)
         if(if.lt.10) then
         ffile=ffile(1:length+11)//aif(2:2)
         else
         ffile=ffile(1:length+11)//aif
         endif

         open(ifile,file=ffile,status='unknown',form='formatted')
!ad
         write(ifile,1061) header(icol)
!ad
         ib1=(if-1)*npf+1
         ib2=ib1+npf-1
         if(ib2.gt.hh) ib2=hh
         write(ifile,360)  (i1(ib),i2(ib),ib=ib1,ib2)
	   do J=2,JEND
           write(ifile,166) eng(j),(DENSTY(IB,J,icol),ib=ib1,ib2)
           enddo
          enddo
	  close(ifile)
       	 enddo
!ad
         ENDIF
!ad
!ad........................written to plot file.........................
!ad
!ad
        IF (ISWITCH.EQ.5) write(6,332)
!ad
         DO 451 J=1,JEND
         WRITE(7,116) ENG(j),(gesDENSTY(J,icol),icol=1,ncol)
 451     WRITE(6,116) ENG(j),(gesDENSTY(J,icol),icol=1,ncol)
!ad
!ad............................sum rule check...........................
!ad
       write(6,401)
       do icol=1,ncol
       sumr(icol)=(gesDENSTY(1,icol)*ENG(1) + &
                  gesDENSTY(JEND,icol)*ENG(jend))/2.d0/sumfac
       enddo
!ad
       do j=2,jend-1
       do icol=1,ncol
       shelp=gesDENSTY(j,icol)*ENG(j)*DE/sumfac
       sumr(icol)=sumr(icol)+shelp
       enddo
!ad    write(55,*) eng(j),(sumr(icol),icol=1,ncol)
       enddo
       write(6,402)  (sumr(icol),icol=1,ncol)  
!ad
!ad...........................end sum rule check........................
!ad
!ad
!ad  ____________________ case DIELECTRIC TENSOR end ___________________
!ad
  900  CONTINUE
!      STOP ' JOINT: LEGAL END'
!
    1 FORMAT(A70)
    2 FORMAT(2I5/3F10.5)
    3 FORMAT(1X,A70/)
    4 FORMAT(I5)
!fb    5 FORMAT(/,2X,'LOWER AND UPPER BAND-INDEX',2X,':',2I5, &
!fb             /,2X,'EMIN, DE, EMAX',14X,':',3F10.5,/) 
    5 FORMAT(/,2X,'LOWER AND UPPER BAND-INDEX',2X,':',2I5,' LAST OCCUPIED BAND-INDEX  :',I5,&
             /,2X,'EMIN, DE, EMAX',14X,':',3F10.5,/) 
   11 FORMAT(A4)
   12 FORMAT(A90)                       
   13 FORMAT(/,1X,' CPU - TIME needed: ',f8.1,//)
   15 FORMAT(1X,A70,/,'NENRG=',I5,//)
   16 FORMAT(f10.5,7(F10.2,f10.4))
  110 FORMAT('#',I1,30x,'Vol = ',F18.10)
  111 format('#',2x,'Energy',A6,5x,8(A9,10x),A9)
  211 format('#',2x,'Energy',A6,7x,'JDOS',14x,'JDOS/E^2',/)
  311 format('#',2x,'Energy',A6,5x,'   DOS    ',/)
  114 FORMAT('#')
  115 FORMAT(3x,'BAND:',I5,3(1X,e18.8))
  116 FORMAT(f13.5,9(1X,e18.8))
  117 FORMAT(2X,'ENERGY',2X,7(5x,I2,1X,A6,6X))
  130 FORMAT(/,3X,'JOINT DOS OF BANDS : ',2I5)
 1130 FORMAT(/,3X,'DIELECTRIC TENSOR COMPONENTS OF BANDS : ',2I5)
  131 FORMAT(/,3X,'DOS OF BAND : ',I5)
  132 FORMAT(//,3X,'TOTAL JOINT DOS: ',/)
  232 FORMAT(//,3X,'TOTAL DOS: ',/)
  332 FORMAT(//,3X,'TOTAL DIELECTRIC TENSOR COMPONENTS: ',/)
  333 format(/,3x,'WARNING: more than 20 files per case ',/, &
       '  data will be overwritten!',/)
  444 format(/,' WARNING: no plasma frequency calculated for column',i2)
  401 format(/,'SUM RULE CHECK: number of electrons'/)
  402 format(3x,9(f8.3))
 9999 CONTINUE
      CALL ERRCLR(ERRFN)
      STOP 'JOINT DOS END'
!
!        error handling
!
  910 INFO = 1
!
!     joint.def couldn't be opened
!
      WRITE (ERRMSG,9000) DEFFN
      CALL OUTERR('JOINT',ERRMSG)
      GOTO 999
!
  920 INFO = 2
!
!     file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('JOINT',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('JOINT',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('JOINT',ERRMSG)
      GOTO 999
!
  960 INFO = 7
!
!     error reading file *.def
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('JOINT',ERRMSG)
  999 STOP 'JOINT - ERROR'
!ad
  77  format(2x,'NUMBER OF BANDS OR BAND COMBINATIONS (', &
                 I4,') BIGGER THAN PARAMETER INMUEden (',I8,')'/, &
       '  WARNING: no band analysis possible! check parameters!',/)
  78  format(2x,'BAND ANALYSIS POSSIBLE: ',I4, ' BANDS/BAND COMBINATIONS: ',i3,' files will be written')
  79  format(/,3x,'parameter NKPT (',I5,') smaller than number of k-points (',I5,')')
!
 330  format(10x,I1,A117)
 166  format(f10.5,10(2X,e10.4))
 266  format(' BAND ANALYSIS FOR JOINT DENSITY OF STATES: ')
 267  format(' BAND ANALYSIS FOR DENSITY OF STATES: ')
 360  format(/,3x,'Energy',10(2x,2i4,2x),/)
 460  format(/,3x,'Energy',10(3x,'BAND:',i3),/)
 461  format(/,3x,'Energy',1(3x,'BAND:',i4),/)
 462  format(/,3x,'Energy',2(3x,'BAND:',i4),/)
 463  format(/,3x,'Energy',3(3x,'BAND:',i4),/)
 464  format(/,3x,'Energy',4(3x,'BAND:',i4),/)
 465  format(/,3x,'Energy',5(3x,'BAND:',i4),/)
 466  format(/,3x,'Energy',6(3x,'BAND:',i4),/)
 467  format(/,3x,'Energy',7(3x,'BAND:',i4),/)
 468  format(/,3x,'Energy',8(3x,'BAND:',i4),/)
 469  format(/,3x,'Energy',9(3x,'BAND:',i4),/)
1061  format(' BAND ANALYSIS FOR DIELECTRIC TENSOR COMPONENT ',A9,':')
6101  format(/,3x,'Energy',A6,5x,'charge_',A2,/)
6102  format(/,3x,'Energy',A6,5x,'charge_',A2,10x,'charge_',A2,/)
6103  format(/,3x,'Energy',A6,5x,'charge_',A2,10x,'charge_',A2, &
       10x,'charge_',A2,/)
 315  format(//,'   BAND ANALYSIS (contributions at Fermi energy): ',/)
 316  format(//,'   Plasma frequencies: ')
 317  format(//,'   Plasma frequencies: ',//, &
       '   !!! WARNING: ',  &
       'w_pl = sqrt( w_pl^2(up-spin) + w_pl^2(dn-spin) )!!!')
 416  format(/,'   Lifetime broadening: ')
 5161 format(/,6x,'gamma_',A2,A6,/)
 5162 format(/,6x,'gamma_',A2,4x,'gamma_',A2,2x,A6,/)
 5163 format(/,6x,'gamma_',A2,4x,'gamma_',A2,4x,'gamma_',A2,2x,A6,/)
 6161 format(/,5x,'  w_p_',A2,2x,A6,/)
 6162 format(/,5x,'  w_p_',A2,5x,' w_p_',A2,2x,A6,/)
 6163 format(/,5x,'  w_p_',A2,5x,' w_p_',A2,5x,' w_p_',A2,2x,A6,/)
 618  format(1x,3f12.4,/)
 612  format('#',2x,'Energy',A6,3x,3(A7,A2,6x))
 614  format(6x,'Energy',6x,6(A6,A2,7x))
 615  format(5x,A6,5x,'[(Ohm cm)-1]')
 619  format(f13.6,12(1x,e14.6))
 789  FORMAT(7x,I4,3x,I4)
 799  FORMAT(2(F16.12))
 809  FORMAT(14x,F12.9,I8)
 819  FORMAT(A30)

 9000 FORMAT(' can''t open definition file ',A40)
 9010 FORMAT(' can''t open unit: ',I2)
 9020 FORMAT(' filename: ',A50)
 9030 FORMAT(' status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)

      END    
