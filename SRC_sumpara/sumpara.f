      PROGRAM SUMPARA                                                     
!                                                                       
!     SUMMATION OF FILES FROM PARALLEL PROCESSING
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!
      CHARACTER*19     nkktext
      CHARACTER*4    IREL
      CHARACTER*5    MIX                                                
      CHARACTER*11   STATUS,FORM                                        
      CHARACTER*20   TITLE                                              
      CHARACTER*60   COMM1,COMM2,COMM3,COMM4                            
      CHARACTER*79   MARGN                                             
      CHARACTER*80   FNAME,CLMFN,SCFFN,SCFALL,scfdmatFN
      CHARACTER*80      :: dmatfn,dmatfn10,dmatfn11,dmatfn12

      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN

!                                                                       
      COMMON /SYMETR/  TAU(3,NSYM),IZ(3,3,NSYM),             &
                       INUM(NSYM),IORD        
      integer,allocatable :: KVCOLD(:,:)                  ! 3,NWAV)
      COMMON /GENER/   BR1(3,3),BR2(3,3),avec(3,3)                    
!para begin

      real*8,allocatable :: CLMOLD(:,:,:),CLMNEW(:,:,:)   ! NRAD,NCOM,NATO)
      complex*16,allocatable :: ROKNEW(:),ROKOLD(:)       ! NWAV)              
      real*8,allocatable :: RNOT(:),DX(:),ZZ(:)           ! NATO)
      integer,allocatable :: JRI(:),LMMAX(:)              ! NATO) 
      integer,allocatable :: LM(:,:,:)                    ! 2,NCOM,NATO)  
      INTEGER             ::idu

      dimension R(NRAD)
!para end  
!
!        Common blocks
!
!
!        LATTIC  - lattice type
!        NAME(i) - name of atom i in the unit-cell (e.g. 'Titanium')
!
      CHARACTER*4        LATTIC
      CHARACTER*10,allocatable ::       NAME(:)            ! NATO)
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
      INTEGER,allocatable ::          IATNR(:), MULT(:)    ! NATO)
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3)
      real*8,allocatable ::               POS(:,:)             ! 3,NDIF)
      real*8,allocatable ::           RMT(:), V(:)         ! NATO)


      COMMON /IPROC/     IPROC
      INTEGER,allocatable ::   ISPLIT(:)                   ! NATO) 
!      COMMON  /STRUK/    POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
!      SAVE    /STRUK/
      COMMON  /STRUK/    ALAT, ALPHA,  PIA, VI 
      SAVE    /STRUK/
      dimension qpw(2)
      real*8,allocatable :: qel(:,:)     ! nato,2)
      real*8,allocatable :: FORSUM(:,:)  ! NATO,4)       
      CHARACTER*10,allocatable :: ANAME(:)     ! NATO)       
!
      real*8,allocatable :: xlsum(:,:,:),ylsum(:,:,:),zlsum(:,:,:)
      complex*16 ds1,dp1(9),dd1(25),df1(49)
      complex*16,allocatable :: dp(:,:,:,:)        ! 9,nato),
      complex*16,allocatable :: ds(:,:,:)          ! NATO),
      complex*16,allocatable :: dd(:,:,:,:)        ! 25,nato)
      complex*16,allocatable :: df(:,:,:,:)        ! 49,NATO)
      INTEGER,allocatable :: lvalue(:,:,:),jatom1(:,:,:)     ! nato)

!---------------------------------------------------------------------  
!                                                                       
!para begin
      CALL GTFNAM(DEFFN,ERRFN,IPROC)

!para end
      CALL ERRFLG(ERRFN,'Error in SUMPARA')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
      idu=0
      dmatFN10 = "This_file_should_be_empty"
      dmatFN11 = "This_file_should_be_empty"
      dmatFN12 = "This_file_should_be_empty"
 783  CONTINUE
         READ (1,*,END=784,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         IF(IUNIT.EQ.8)  scfdmatFN=FNAME
         IF(IUNIT.EQ.17) CLMFN=FNAME
         IF(IUNIT.EQ.21) SCFFN=FNAME
         IF(IUNIT.EQ.22) SCFALL=FNAME

         IF(IUNIT.EQ.10) dmatFN10=FNAME
         IF(IUNIT.EQ.11) dmatFN11=FNAME
         IF(IUNIT.EQ.12) dmatFN12=FNAME
         IF(IUNIT.EQ.11) idu=2
      GOTO 783
 784  CONTINUE
      CLOSE (1)
!test
      write(6,*)'SUMPARA running from ',IPROC,'input files'
      write(6,*)' '
      PI=ACOS(-1.0d0)                                                     
!                                                                       
!.....READ TAPE20=POTE                                                  
!                                                                       
      READ(20,1000) TITLE                                               
      READ (20,5010) LATTIC, NAT, IREL
!
      allocate ( CLMOLD(NRAD,NCOM,NAT),CLMNEW(NRAD,NCOM,NAT))
      allocate ( RNOT(nat),DX(nat),ZZ(NAT))
      allocate ( JRI(nat),LMMAX(nat))
      allocate ( LM(2,NCOM,NAT))  
      allocate ( NAME(NAT))
      allocate ( IATNR(nat), MULT(nat))
      allocate ( RMT(nat), V(nat))
      allocate ( ISPLIT(nat))
      allocate ( qel(nat,2))
      allocate ( FORSUM(NAT,4))       
      allocate ( ANAME(NAt))       

      allocate ( dp(9,nat,4,0:2))
      allocate ( ds(NAT,4,0:2),xlsum(nat,4,0:2),ylsum(nat,4,0:2),zlsum(nat,4,0:2))
      allocate ( dd(25,nat,4,0:2))
      allocate ( df(49,NAT,4,0:2))
      allocate ( lvalue(nat,4,0:2),jatom1(nat,4,0:2))
!

      jatom1=0
      xlsum=0.d0
      ylsum=0.d0
      zlsum=0.d0
      ds=(0.d0,0.d0)                                          
      dp=(0.d0,0.d0)                                          
      dd=(0.d0,0.d0)                                          
      df=(0.d0,0.d0)                                          

      CLMOLD(1:nrad,1:ncom,1:nat)=0.0d0 
      CLMNEW(1:nrad,1:ncom,1:nat)=0.0d0 

!
!        read in lattice constants
!
      READ (20,5020) ALAT(1), ALAT(2), ALAT(3), &
                     ALPHA(1), ALPHA(2), ALPHA(3)
      IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
      IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
      IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
      call latgen(nat)
        assign 2021 to iform1
        assign 2071 to iform2
        assign 2076 to iform3

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
      allocate (POS(3,48*nat))
      DO 5 JATOM = 1,NAT                                                
         INDEX=INDEX+1                                                  
!
         READ(20,1040) IATNR(JATOM),POS(1,INDEX),POS(2,INDEX), &
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
         READ(20,1050) ANAME(JATOM),JRI(JATOM),RNOT(JATOM),RMT(jatom), &
                 ZZ(JATOM)                                        
!para begin
         DX(JATOM)=LOG(RMT(jatom)/RNOT(JATOM)) / (JRI(JATOM)-1)               
!para end
         READ(20,*)                                                     
         READ(20,*)                                                     
         READ(20,*)                                                     
         V(JATOM)=4.D0*PI/3.D0*RMT(JATOM)**3
 5    CONTINUE                                                          




!para
!     read scf file and sums up total energy and forces
!
      write(6,*)'Summing ',IPROC,' SCF Files:'
      CALL SCFSUM(SCFFN,SCFALL,nat)

!
!                                                                       


!para loop over all clm_$i files
      
      DO II=1,IPROC
         write(6,*)'CLMVAL summation of file ',II
         close(17)
         call mknam(FNAME,CLMFN,II)
         OPEN(17,FILE=FNAME,STATUS='old',FORM='formatted',err=180)
! inlude Error!!!!!
!         write(*,*)FNAME

         READ(17,2032,err=180,end=180)                                        
                                     
!                                                                       
!.....Read  DENSITY IN SPHERES                                            
!
!para begin 
         DO 15 JATOM=1,NAT   
            do 16 jrj=1,jri(jatom)
 16            r(jrj)=rnot(jatom)*exp(dx(jatom)*(jrj-1))
!para end
!     
            READ(17,1980)                                               
            READ(17,2000) LMMAX(JATOM) 
            if(lmmax(jatom).gt.ncom) stop 'ncom too small'
            DO 20 LM1=1,LMMAX(JATOM)                                    
               READ(17,2010) LM(1,LM1,JATOM),LM(2,LM1,JATOM)               
               READ(17,iform1) (CLMOLD(J,LM1,JATOM),J=1,JRI(JATOM))
! sum densities from single processor files
               DO J=1,JRI(JATOM)
                  CLMNEW(J,LM1,JATOM)=CLMNEW(J,LM1,JATOM)+ &
                       CLMOLD(J,LM1,JATOM)
               ENDDO
!

 20            READ(17,2031)                                               
               READ(17,2033)                                               
                                                                      
 15      CONTINUE 

!                                                                  
!.....READ DENSITY IN INTERSTITIAL                                       
!                                                                       
         READ(17,1980)                                                  
!         READ(17,2060) NKKVL  
  read(17,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkkvl
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkkvl
 6768 continue
         if(nkkvl.eq.0) cycle
         if(ii.eq.1) then
! allocate only for the first time
           allocate ( KVCOLD(3,Nkkvl))
           allocate ( ROKNEW(nkkvl),ROKOLD(nkkvl))
           DO 2 J    =1,Nkkvl                                                 
 2         ROKNEW(J)=(0.0d0,0.0d0)
         endif                                      
         DO 140 J=1,NKKVL                                               
            READ(17,iform2)  (KVCOLD(JX,J),JX=1,3),ROKOLD(J)                
            ROKNEW(J)=ROKNEW(J)+ROKOLD(J)
 140     CONTINUE
! loop over all processes                                                    
      ENDDO 

! Now write to sums to the master clmsum-file
      close(17)
!      write(*,*)CLMFN
      OPEN(17,FILE=CLMFN,STATUS='unknown',FORM='formatted')
!                                                                       
      WRITE(17,1970)                                               
      WRITE(17,78) TITLE,LATTIC,ALAT(1),ALAT(2),ALAT(3),JRI(1)
      WRITE(17,77)                               
      DO 60 JATOM=1,NAT                                                 
         WRITE(17,1990) JATOM                                           
         WRITE(6,1990) JATOM                                           
         WRITE(17,2001) LMMAX(JATOM)                                    
         DO 65 LM1=1,LMMAX(JATOM)                                       
            WRITE(17,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
            WRITE(6,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
            WRITE(17,iform1) ( CLMNEW(J,LM1,JATOM), J=1,JRI(JATOM) )    
 65         WRITE(17,2031)                                              
 60   WRITE(17,2033)                                                    
!                                                                       
      WRITE(17,2051)                                                    
      WRITE(17,1980)                                                    
      WRITE(17,2061) NKKVL                                             
      DO 155 J=1,NKKVL                                                 
 155  WRITE(17,iform3) (KVCOLD(JX,J),JX=1,3),ROKNEW(J)      

 180  continue
!
!  sum dmat for LDA+U
!
      DO II=1,IPROC
         write(6,*)'DMAT summation of file ',II
         DO ifile=0,idu
            close(10+ifile)
            IF(ifile==0) dmatfn=dmatfn10
            IF(ifile==1) dmatfn=dmatfn11
            IF(ifile==2) dmatfn=dmatfn12
            call mknam(FNAME,dmatFN,II)
            OPEN(10+ifile,FILE=FNAME,STATUS='unknown',FORM='formatted',err=200)
            write(6,*)FNAME

            inatmold=0
            jatom=0 
9999        continue
            READ(10+ifile,*,end=199) indatm
            if (indatm.eq.0) then
               goto 199   
            elseif (indatm.eq.inatmold) then
               iorb=iorb+1
            else
               iorb=1
               indatmold=indatm
               jatom=jatom+1
            endif
            jatom1(jatom,iorb,ifile)=indatm          
            READ(10+ifile,*)  lvalue(jatom,iorb,ifile),xl,yl,zl
            xlsum(jatom,iorb,ifile)=xlsum(jatom,iorb,ifile)+xl
            ylsum(jatom,iorb,ifile)=ylsum(jatom,iorb,ifile)+yl
            zlsum(jatom,iorb,ifile)=zlsum(jatom,iorb,ifile)+zl
            if(lvalue(jatom,iorb,ifile).eq.0) then                           
               read(10+ifile,'(2e16.8)') ds1
               write(6,*) jatom,ds1
               ds(jatom,iorb,ifile)=ds(jatom,iorb,ifile)+ds1
            else if(lvalue(jatom,iorb,ifile).eq.1) then  
               read(10+ifile,201) dp1
               if(jatom.eq.1) write(6,*) dp1
               do i=1,9
                  dp(i,jatom,iorb,ifile)=dp(i,jatom,iorb,ifile)+dp1(i)
               enddo
            else if(lvalue(jatom,iorb,ifile).eq.2) then 
               read(10+ifile,204) dd1
               do i=1,25
                  dd(i,jatom,iorb,ifile)=dd(i,jatom,iorb,ifile)+dd1(i)
               enddo
            else if(lvalue(jatom,iorb,ifile).eq.3) then 
               read(10+ifile,205) df1
               do i=1,49
                  df(i,jatom,iorb,ifile)=df(i,jatom,iorb,ifile)+df1(i)
               enddo
            endif
201         format(2e16.8,2x,2e16.8,/,2e16.8)
204         format(2(2e16.8,2x,2e16.8,/),2e16.8)
205         format(3(2e16.8,2x,2e16.8,/),2e16.8)
            goto 9999
199       continue
            enddo
! loop over all processes                                                    
         ENDDO
                                                         
! Now write to sums to the master dmat-file
         DO ifile=0,idu
            close(10+ifile)
            IF(ifile==0) dmatfn=dmatfn10
            IF(ifile==1) dmatfn=dmatfn11
            IF(ifile==2) dmatfn=dmatfn12
            OPEN(10+ifile,FILE=dmatFN,STATUS='unknown',FORM='formatted')
            DO  JATOM=1,NAT
               if(jatom1(jatom,1,ifile).eq.0) exit   
               do iorb=1,4
                  if(jatom1(jatom,iorb,ifile).eq.0) exit   
                  WRITE(10+ifile,203) jatom1(jatom,iorb,ifile),mult(jatom1(jatom,iorb,ifile))
                  write(10+ifile,202) lvalue(jatom,iorb,ifile),xlsum(jatom,iorb,ifile), &
                       ylsum(jatom,iorb,ifile),zlsum(jatom,iorb,ifile)
                  if(lvalue(jatom,iorb,ifile).eq.0) then
                     write(10+ifile,201) ds(jatom,iorb,ifile)
                  else if(lvalue(jatom,iorb,ifile).eq.1) then
                     write(10+ifile,201) (dp(i,jatom,iorb,ifile),i=1,9)
                  else if(lvalue(jatom,iorb,ifile).eq.2) then
                     write(10+ifile,204) (dd(i,jatom,iorb,ifile),i=1,25)
                  else if(lvalue(jatom,iorb,ifile).eq.3) then
                     write(10+ifile,205) (df(i,jatom,iorb,ifile),i=1,49)
                  endif
202               format(i5,3f10.6,' L, Lx,Ly,Lz in global orthogonal system')
203               format(2i5,' atom density matrix, multiplicity')
        enddo
     enddo
enddo
 200     continue
!para
!     read scfdmat file and sums up orb moment + density matrix
!
      write(6,*)'Summing ',IPROC,' SCFdm Files:'
      CALL dmatSCF(SCFdmatFN,nat)

                                     
! ..................................................................
!                                                                       
!      WRITE(21,2074) ESUM                                               
!      WRITE(6,2074) ESUM                                                
!      if(force(jatom)) WRITE(21,720)
!      if(force(jatom)) WRITE(6,720)
!      DO 5720 JATOM=1,NAT
!      if(force(jatom)) then
!         WRITE(21,710) JATOM,JATOM,(FORSUM(JATOM,j),j=1,4)
!         WRITE(6,710) JATOM,JATOM,(FORSUM(JATOM,j),j=1,4)
!      endif
! 5720 CONTINUE
      CALL ERRCLR(ERRFN)
      STOP ' SUMPARA END'                                                  
!
!        error handling
!
  910 INFO = 1
!
!        'SUMPARA.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('SUMPARA',ERRMSG)
      GOTO 999
!
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('SUMPARA',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('SUMPARA',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('SUMPARA',ERRMSG)
      GOTO 999
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      WRITE (ERRMSG,9050) JATOM, MULT(JATOM), INDEX
      CALL OUTERR('SUMPARA',ERRMSG)
      CALL OUTERR('SUMPARA','MULT .EQ. 0')
      GOTO 999
  940 INFO = 4
      CALL OUTERR('SUMPARA','Maximum number of K-points exceeded.')
      GOTO 999
!                                                                       
  960 INFO = 7
!
!        Error reading file 'sumpara.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('SUMPARA',ERRMSG)
      GOTO 999
  999 STOP 'SUMPARA - Error'
!                                                                       
 77   FORMAT(1X,9(F9.7))                                                
 78   FORMAT(1X,A20,A4,3F10.6,I5)                                       
 700  FORMAT(I3,A76)                                                    
 701  FORMAT(I3,A4)                                                     
 702  FORMAT(A5,13X,A60)                                                
 703  FORMAT(A8,F8.5,A60)                                               
 704  FORMAT('CONTINUE',F8.5,A60)                                       
 705  FORMAT(F10.3,6X,A60)                                              
 710  FORMAT(':FOR',i2.2,':',1x,i2,'.ATOM',4F15.3)
 720  FORMAT (7x,'TOTAL-FORCE IN mRy/a.u. = |F|',5x,'Fx',13x, &
              'Fy',13x,'Fz')
 743  FORMAT(3X,A76)                                                    
 1000 FORMAT(A20)                                                       
 1040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1041 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1055 FORMAT(I4)                                                        
 1056 FORMAT(3(3I2,F5.2,/),I8)                                          
 1061 FORMAT(7X,A5,' MIXING SCHEME WITH',F6.3)                          
 1070 FORMAT(I3,'. ITERATION',//)                                       
 1970 FORMAT(3X,'             TOTAL CHARGE DENSITY GENERATED ',        &
             'BY',I3,'. ITERATION ')                 
 1980 FORMAT(3X)                                                        
 1990 FORMAT(3X,'ATOM NUMBER=',I3,5X,10A4)                             
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
 2074 FORMAT(//,':ENE  :',1X,'********** TOTAL ENERGY IN Ry =',F20.6/) 
 5010 FORMAT(A4,23X,I3,/,13X,A4)
 5020 FORMAT(6F10.7)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
 9050 FORMAT('MULT(',I3,'=',I3,', INDEX=',I3)
      END                                                               
