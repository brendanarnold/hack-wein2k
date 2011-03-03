      PROGRAM LAPW2                                                     
!                                                                       
!                  D E N S I T Y   G E N E R A T I O N                  
!                                                                       
      USE param; USE defs
      USE parallel
      USE reallocate
      USE ams
      USE bandm
      USE kpp1
      USE struk; USE xa2; USE com
      USE char
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef Parallel
  INCLUDE 'mpif.h'
#endif
      REAL*8, allocatable       ::  ZZ(:)
      REAL*8           GMAX
      CHARACTER*4      IREL,RELA,cform,RCFILE
      CHARACTER*5      COORD
      CHARACTER*10    KNAME
      CHARACTER*11     STATUS,FORM                                      
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN
      CHARACTER*77     RMARGN                                           
      CHARACTER*80     FNAME,FNAME1,fnamehelp
!para begin
      CHARACTER*80     VECFN, VECFND,WEIFN,WEIFND,enefn, enefnd
      INTEGER          IPROC
!para end

      LOGICAL         SO
      common /mkef/delef,ts2
!para begin
      COMMON /PROC/    VECFN,VECFND,enefn,enefnd
      COMMON /IPROC/   IPROC 
!para end
      DATA             RELA/'RELA'/,SO/.false./                           
      save TTOTAL,TFERMI,TCLM,TFOUR,tstart,errfn                          
      
      nspin1=2
      nlo=0
      RCFILE='NOFI'
!-----------------------------------------------------------------------  
!                                                                       
      CALL init_ams
      CALL GTFNAM(DEFFN,ERRFN,iproc)
      CALL ERRFLG(ERRFN,'Error in LAPW2')
!      write(*,*)DEFFN,ERRFN,IPROC
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         IF(iunit.NE.6.AND.iunit.NE.8) THEN
            OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         ELSE
            IF(myid.EQ.0) OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         ENDIF
         if(iunit.eq.18) then
           do i=80,5,-1
             if(fname(i:i).ne.' ') then
               if(fname(i-2:i).eq.'vsp') nspin1=1
               goto 10
             endif
           enddo
         endif
         if(iunit.eq.12) then
           do i=80,6,-1
             if(fname(i:i).ne.' ') then
               if(fname(i-5:i).eq.'normso') SO=.true.
!             goto 10
             endif
           enddo
         endif
         if(iunit==31) fnamehelp=fname
         if(iunit.eq.13) fname1=fname
         close(13)
!para begin
! remember head of vector file for subsequent opening
! of splitted case.vector_*
         if(iunit.eq.9)  VECFND=FNAME
         if(iunit.eq.10) VECFN=FNAME
         if(iunit.eq.26) WEIFN=FNAME
         if(iunit.eq.27) WEIFND=FNAME
         if(iunit.eq.29)  eneFND=FNAME
         if(iunit.eq.30) eneFN=FNAME         
!para end
      GOTO 10
   20 CONTINUE
      CLOSE (1)

!para begin
      IF(myid.EQ.0) THEN
         IF(IPROC.GT.0)THEN
            write(6,*)'Running LAPW2 on ',IPROC,' processors'
            write(6,*)' '
         ELSE
            write(6,*)'Running LAPW2 in single processor mode'
            write(6,*)' '
         END IF
      ENDIF
!para end

!                                                                       
!.....READ STRUCT                                                       
      READ(20,1000) TITLE                                               
      READ(20,1010) LATTIC,NAT,cform,IREL    
      allocate (rmt(nat),v(nat),iatnr(nat),mult(nat),isplit(nat))
      allocate (rotloc(3,3,nat))
      allocate (r0(nat),dx(nat),jri(nat))
      allocate (aname(nat))
      allocate (ZZ(nat))
      allocate ( pos(3,48*nat) )
!     READ IN LATTICE CONSTANTS                                         
      READ(20,1020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      IF(ABS(ALPHA(1)).LT.TEST) ALPHA(1)=NINETY
      IF(ABS(ALPHA(2)).LT.TEST) ALPHA(2)=NINETY
      IF(ABS(ALPHA(3)).LT.TEST) ALPHA(3)=NINETY
 IF(IREL.EQ.'RELA') REL=.TRUE.                                     
      IF(IREL.EQ.'NREL') REL=.FALSE.                                    
!
!        read crystal-structure (atompositions, symmetry-parameters,
!                                muffin-tin radius, ...)
!        'INDEX' counts all atoms in the unit cell,
!        'JATOM' counts only the non-equivalent atoms
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1
!
         READ(20,1030) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM),  &
                       ISPLIT(JATOM)                                    
         if(isplit(jatom).eq.99.and.lxdos.ne.3) goto 956
         if(isplit(jatom).eq.88.and.lxdos.ne.3) goto 956
         IF (MULT(JATOM) .EQ. 0) THEN
!
!           illegal number of equivalent atoms
!
            IF(myid.EQ.0) WRITE (6,6000) JATOM, INDEX, MULT(JATOM)
            GOTO 930
         ENDIF
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
!
               READ(20,1031) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
 55         CONTINUE                                                    
         READ(20,1050) ANAME(JATOM),JRI(JATOM),R0(JATOM),RMT(JATOM), &
                       ZZ(JATOM)
         DX(JATOM)=LOG(RMT(JATOM)/R0(JATOM)) / (JRI(JATOM)-1)           
         RMT(JATOM)=R0(JATOM)*EXP( DX(JATOM)*(JRI(JATOM)-1) )           
         READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)                
 1051    FORMAT(20X,3F10.8)                                             
 50   CONTINUE                                                          
      ALLOCATE (rotij(3,3,index),tauij(3,index))
      ndif=index
      CALL doreallocate(pos, 3, ndif)
!                                                                       
!.....READ INPUT AND POTE                                               
      READ(5,1003)  MODUS,MODUS1,coord                                        
      IF(myid.EQ.0) write(6,*)  ' Modus: ', MODUS
      IF(myid.EQ.0) then
         IF(MODUS1.eq.'EECE ') write(6,*)' Density for selected electrons only'
      endif
      esepermin=0.5d0
      eseper0=0.05d0
      READ(5,*,err=4716)  EMIN,ELECN,esepermin,eseper0
      if(myid.eq.0) write(6,*) 'no read error'
 4716 if(esepermin.lt.1.d-3)   esepermin=0.5d0
      READ(5,1004) efmod,ef
      IF(MODUS1.eq.'EECE ') then
       read(5,*)natu         !EECE
        do i=1,natu          !EECE
         read(5,*)jatu,nlatu !EECE
        enddo                !EECE
       endif
      DO I=1,NAT
         READ(5,*) IDUMMY
      ENDDO
      READ(5,*,END=4715) GMAX
      READ(5,'(A4)',END=4715) RCFILE
 4715 CONTINUE
      REWIND 5
 1234 FORMAT(//,1A)
      READ(5,1234) DUMMY
       if (myid.eq.0) write(21,720)gmax
  720  FORMAT(':GMA  :',' POTENTIAL AND CHARGE CUT-OFF',f7.2,' Ry**.5')

      IF(myid.EQ.0) WRITE(6,*) 'RECPR LIST: ',RCFILE
     
      ts2=0.d0
!                                                                       
      IF(myid.EQ.0) THEN
         WRITE(6,800)                                                      
         WRITE(6,805)  TITLE                                               
         WRITE(6,810)  LATTIC                                              
         WRITE(6,820)  AA,BB,CC                                            
         WRITE(6,840)  NAT                                                 
         WRITE(6,850)  IREL                                                
!     WRITE(6,860)  ISCF                                                
         WRITE(6,870)  COORD
      ENDIF
      REL=.FALSE.                                                       
      IF(IREL.EQ.RELA) REL=.TRUE.                                       
      CALL LATGEN (NAT)                                                 
!     TSTART=SECOND()                                                   
      CALL CPUTIM(TTIME)
      TSTART=TTIME                                        

      CALL RECFIL(FNAME1,GMAX,RCFILE, &
                   KXMAX,KYMAX,KZMAX,NWAVE)

      if((IPROC.EQ.0).OR.(MODUS.EQ.'FERMI') &
         .OR.(MODUS.EQ.'QTL  ').OR.(MODUS.EQ.'EFG  '))then
!para end
!.....CALCULATE FERMI-ENERGY AND WEIGHT OF K-POINTS                     
!         IF(MODUS1.eq.'SO   ') elecn=elecn*2.0d0
         nbmax=-5
         CALL FERMI(nbmax)
         eseper=-999.
         eftest=ef 
         if(efmod.eq.'TETRA') eftest=ef-0.5d0
         if(efmod.eq.'GAUSS') eftest=ef-0.5d0
!     Added TEMPS mode
         if(efmod(1:4).eq.'TEMP') eftest=ef+delef

         IF(myid.eq.0) THEN
            write(6,*)
            write(6,*) 'Bandranges (emin - emax):'
            write(21,*) '        Bandranges (emin - emax):'
            do i=1,nume
               if(ebmax(i).lt.-990.d0) exit
               write(6,*) 'band',i,ebmin(i),ebmax(i)
               if(i.le.nbmax+5.and.i.ge.nbmax-10) write(21,'(a9,i4,2f12.6)') ':BAND   :',i,ebmin(i),ebmax(i)
            enddo
            write(6,*)
            write(6,*) 'eseper below EF',esepermin
            write(6,*) 'eseper minimum gap',eseper0
         ENDIF
         do i=nume-1,2,-1
         if(ebmax(i).lt.-990.d0) cycle
!         write(6,*) 'band',i,ebmin(i),ebmax(i)
!         if(ebmin(i).lt.ebmax(i)) then
           if(ebmax(i-1).lt.eftest - esepermin ) then
!... requires band gaps of more than eseper0 and Esepermin below EF for e-seperation
             if(ebmin(i)-eseper0 .gt. ebmax(i-1)) then
              eseper=ebmin(i)-eseper0
              goto 119
             endif
           endif
!         endif
         enddo
119      if(myid.eq.0) write(6,*) 'Energy to separate semicore and valence states:' &
                    ,eseper
         if(myid.eq.0) write(21,121) eseper
121      format('        Energy to separate low and high energy' &
               ,'states: ',f10.5)
         IF(MODUS1.eq.'AF   ') then
            do 120 i1=1,nkpt-1
               do 120 i2=1,nume
                  WEIGH(i1,i2)=WEIGH(i1,i2)/2.0d0
 120        continue
         endif
         IF(SO) then
            do 110 i1=1,nkpt-1
               do 110 i2=1,nume
                  WEIGH(i1,i2)=WEIGH(i1,i2)*2.0d0/NSPIN1
 110        continue
         endif
!         IF(MODUS1.eq.'SO   ') then
!            elecn=elecn/2.0d0
!            do 110 i1=1,nkpt-1
!               do 110 i2=1,nume
!                  WEIGH(i1,i2)=WEIGH(i1,i2)/2.0d0
! 110        continue
!         endif
!         write(6,*) 'zuruck von fermi'
!     write(6,*) elecn,ef,delef,efmod                                     
!     IF(MODUS1.eq.'AF   ') elecn=elecn/2.0d0
         
         if(efmod.eq.'TETRA') then
            if(myid.eq.0) WRITE(6,1061)  ELECN,EF-0.5d0 
            if(myid.eq.0) WRITE(21,1061) ELECN,EF-0.5d0
            DELEF=-0.5d0
         elseif(efmod.eq.'GAUSS') then
            if(myid.eq.0) WRITE(6,1062)  ELECN,EF-0.5d0 
            if(myid.eq.0) WRITE(21,1062) ELECN,EF-0.5d0 
            DELEF=-0.5d0
!       Added TEMPS mode
         elseif(efmod(1:4).eq.'TEMP') then
            if(myid.eq.0) WRITE(6,1063)  ELECN,EF+delef
            if(myid.eq.0) WRITE(21,1063) ELECN,EF+delef 
         else
            if(myid.eq.0) WRITE(6,1060)  ELECN,EF
            if(myid.eq.0) WRITE(21,1060) ELECN,EF    
            DELEF=0.d0                                       
         endif
         

!para begin   
!     running on single processor?
         if(IPROC.EQ.0) then
!            write(6,*)'single processor'
               write(26)EF,eseper
                      wmist=0.d0
                do j=1,nkpt-1
                  write(26)NE(J)
                  write(26)(WEIGH(J,N), N=1,NE(J))
                      do n=1,ne(j)
                      wmist=wmist+WEIGH(J,N)
                      enddo
                enddo
                IF(myid.EQ.0) write(6,*) 'wmist',wmist
            GOTO 4711  
         endif
         if(iproc.gt.0) then
            
!     write weighs for spin up to files
            write(6,*)WEIFN,iproc
            DO I=1,IPROC
               close(26)
               call mknam(FNAME,WEIFN,I)
               OPEN(26,FILE=FNAME,STATUS='unknown',FORM='unformatted')
               write(26)EF,eseper
               write(6,*)'writing ',fname
               if(I.eq.1) then
                  ifrom=1
               else
                  ifrom=KPP(i-1)
               endif
!               write(6,*)'Range: ',ifrom,kpp(i)
               DO J=IFROM,KPP(I)-1
!                  write(6,*) j,': from:',ifrom,'to:',kpp(i),
!     *                                       'N=',NE(J)
                  write(26)NE(J)
                  write(26)(WEIGH(J,N), N=1,NE(J))
               ENDDO
            ENDDO

!     test if we're  running spin polarized; if so, write the
!      weighs for spin down 
            if(KPP(IPROC+1).GT.0) then
               write(6,*)'Spin polarized'
               DO I=IPROC+1,2*IPROC
                  close(27)
                  call mknam(FNAME,WEIFND,I)
                  OPEN(27,FILE=FNAME,STATUS='unknown', &
                       FORM='unformatted')
                  write(27)EF,eseper
                  write(6,*)'writing ',fname
                  ifrom=KPP(i-1)+1
                  DO J=IFROM,KPP(I)
                     write(27)NE(J)
                     write(27)(WEIGH(J,N), N=1,NE(J))
                  ENDDO
               ENDDO
            else
!               write(6,*)'Not Spin polarized'
            endif
	    if(MODUS.EQ.'QTL  ') GOTO 4711
	    if(MODUS.EQ.'EFG  ') GOTO 4711

! stop removed here for optic parallel weights
!            CALL ERRCLR(ERRFN)
!            STOP'LAPW2 - FERMI; weighs written'     

         endif
      else
! running on multiple processors, read in weighs
! here!           
      DO I=1,NAT 
         READ(30,'(f9.5)') EMIST
         READ(30,'(f9.5)') EMIST
      ENDDO
     k=0  
       DO
         READ(30,'(3e19.12,a10,2i6)',IOSTAT=ios) S,T,Z,KNAME,N,NEn
         IF (ios /= 0) EXIT
         k=k+1
         nmat=MAX(n,nmat)
         nume=MAX(nen,nume)
         DO ii=1,nen
            READ(30,*) NUM,E1
         ENDDO
      ENDDO
      nkpt=k+1
      REWIND(30)
      CALL init_bandm(nume)
      CALL init_com(nkpt,nume)
      CALL init_kpp(nkpt)
      CALL init_xa2(nume,nkpt)

         if(myid.eq.0) write(6,*)'reading weighs from files'
         k=0
         READ(26) EF,eseper
         if(myid.eq.0) write(6,*)'EF=',EF
 4712    k=k+1
         READ(26,END=4713) NN
         if(myid.eq.0) WRITE(6,*)'K=',k,' N=',NN
         READ(26) (weigh(k,n),N=1,NN)
         GOTO 4712
 4713    CONTINUE
!         write(6,*)'done reading weighs'
      endif
! 
 4711 CONTINUE                                                          

!para end

!
!*******roa************** needed for weight-file ***********************************
!
      IF(MODUS.EQ.'FERMI') THEN
         IF(myid.eq.0) write(6,*) 'NK',nk

      IF(EFMOD.EQ.'TETRA') THEN
      nemax=0
      do ik=1,2*nk
        if(ne(ik).gt.nemax) nemax=ne(ik)
      enddo
!
!
      WRITE(11,81) TITLE
!      if(nspin.eq.2) WRITE(92,81) TITLE
      if(efmod.eq.'TETRA') then 
        write(11,80) ef-0.5d0,nk
!        if(nspin.eq.2) write(92,80) ef-0.5d0,nk
      else
        WRITE(11,80) EF,NK
!        if(nspin.eq.2) WRITE(92,80) EF,NK
      endif
         DO K=1,NK

           WRITE(11,78) K,NE(K)
!           if(nspin.eq.2) WRITE(92,78) K,NE(K+NKPT)
           DO IE=1,NE(K)
               WRITE(11,79)  E((k-1)*nemax*nspin+ie),  &
               WEIGHT((k-1)*nemax*nspin+ie)*2/nspin
           END DO
!           if (nspin.eq.2) then
!             DO IE=1,NE(K+NKPT)
!             WRITE(92,79)  E((k-1)*nemax*nspin+ie+nemax), 
!     &       WEIGHT((k-1)*nemax*nspin+ie+nemax)*2/nspin
!             END DO
!           endif
         END DO
!.....end do of nk .......
      ENDIF
!i...endif of tetra...
      
      IF((EFMOD.EQ.'ROOT ').OR.(EFMOD.EQ.'ALL  ')) THEN
      WRITE(11,81) TITLE
      WRITE(11,80) EF,NK
         DO K=1,NK
           WRITE(11,78) K,NE(K)
           eeeo=-1000.
           shift=-10.
           DO IE=1,NE(K)
               wei=WEIGH(k,ie)
               eee= E(k+(ie-1)*2*nkpt)
               IF (eee.lt.eeeo) shift=0. 
               eeeo=eee
               eee=eee+shift
               IF (eee.gt.EF) wei=0.
               WRITE(11,79) eee,wei 
           END DO
         END DO
      ENDIF
!....endif of root and all..........
            if (myid == 0) CALL ERRCLR(ERRFN)
            STOP 'LAPW2 - FERMI; weighs written'     
      ENDIF 
 78   FORMAT('K-PNT',I6,':',2x,I4,' BANDS')
 79   FORMAT(2(F16.12))
 80   FORMAT('FERMI ENERGY: ',F12.9,I8,' K-POINTS')
 81   FORMAT(A30)

      CALL CPUTIM(TTIME)
      TFERMI=TTIME                                        
!.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES      
      IF(MODUS.ne.'FOURI') call l2main(cform,nlo,coord,zz,so,NSPIN1,fnamehelp)
      CALL CPUTIM(TTIME)
      TCLM=TTIME

!.....CALCULATE CHARGE DENSITY RHOK(K) IN INTERSTITIAL                  
      IF(MODUS.EQ.'QTL  ') GOTO 1213
      IF(MODUS1.EQ.'EECE ')then        !EECE
       NWAVE=0                         !EECE
  WRITE(8,*) '   VALENCE CHARGE DENSITY IN INTERSTITIAL '
  WRITE(8,*)
       WRITE(8,2061)  NWAVE            !EECE
2061   FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')    !EECE
      GOTO 1213                        !EECE
      endif
      IF(MODUS.NE.'CLM  ')  &
             CALL FOURIR(nlo,nwave)
 1213 CALL CPUTIM(TTIME)
      TFOUR=TTIME                                          

!.....CALCULATE CPUTIME REQUIRED                                        
      TTOTAL=TFOUR-TSTART                                               
      TFOUR=TFOUR-TCLM                                                  
      PFOUR=TFOUR/TTOTAL*100.                                           
      TCLM=TCLM-TFERMI                                                  
      PCLM=TCLM/TTOTAL*100.                                             
      TFERMI=TFERMI-TSTART                                              
      PFERMI=TFERMI/TTOTAL*100.                                         
      IF(myid.EQ.0) THEN
         WRITE(6,2000)                                                     
         WRITE(6,2010) TTOTAL,100.0
         WRITE(6,2020) TFERMI,PFERMI
         WRITE(6,2030) TCLM,PCLM
         WRITE(6,2040) TFOUR,PFOUR
      ENDIF
      if (myid == 0) CALL ERRCLR(ERRFN)
#ifdef Parallel
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_FINALIZE(ierr)
#endif
      STOP ' LAPW2 END'                                                 
!
!        error handling
!
  910 INFO = 1
!
!        'lapw2.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('LAPW2',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('LAPW2',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('LAPW2',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('LAPW2',ERRMSG)
      GOTO 999
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      CALL OUTERR('LAPW2','MULT .EQ. 0')
      GOTO 999
  956 INFO = 56
      CALL OUTERR('LAPW2','LXDOS must be 3 for ISPLIT=99 or 88')
      CALL OUTERR('LAPW2','LXDOS must be 3 for ISPLIT=99 or 88')
      GOTO 999
  960 INFO = 7
!
!        Error reading file 'lapw2.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('LAPW2',ERRMSG)
      GOTO 999
  999 STOP 'LAPW2 - Error. Check file lapw2.error'
!                                                                       
!                                                                       
  43  FORMAT(3X,A77)                                                    
  44  FORMAT(I3,A77)                                                    
 700  FORMAT(I3,A77)                                                    
 800  FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ',            &
             'I N F O R M A T I O N',/,30X,50(1H-),//)                  
 805  FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
 810  FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
 820  FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
 830  FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)                      
 840  FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
 850  FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    
 860  FORMAT(3X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)                 
 870  FORMAT(3X,'TYPE OF COORDINATES IN DSPLIT= ',A5)                   
 1000 FORMAT(A80)                                                       
 1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
 1020 FORMAT(6F10.7,10X,F10.7)                                          
 1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1040 FORMAT(///,3X,'ERROR IN LAPW2 : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1002 FORMAT(3F10.5,I5)                                                 
 1003 FORMAT(A5,a5,a5)                                                        
 1004 FORMAT(A5,f10.5)                                                        
 1060 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,            &
             //   ,':FER  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)            
 1061 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(TETRAH.M.)','= ',F9.5)           
 1062 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(GAUSS-.M.)','= ',F9.5)           
 1063 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(FERMI-SM.)','= ',F9.5)           
 2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)                       
 2010 FORMAT(12X,'TOTAL       : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2020 FORMAT(12X,'PART FERMI  : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2030 FORMAT(12X,'PART CLM    : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2040 FORMAT(12X,'PART FOURIR : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 6000 FORMAT(///,3X,'ERROR IN LAPW2 : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
      END                                                               
