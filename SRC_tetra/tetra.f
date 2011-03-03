     PROGRAM TETRA                                                     
!                                                                       
!               D E N S I T Y    O F    S T A T E S                     
!               ACTUAL VERSION  with BLOECHL SCHEME                
!                                                                       
      use reallocate
      INCLUDE 'param.inc'
      CHARACTER *7  A7                                                  
      CHARACTER *9  A9,B9                                               
      CHARACTER *20 TITLE                                               
      CHARACTER*45,allocatable :: TEXT(:)                             
      CHARACTER *70 SYSTEM                                              
      CHARACTER*11  FORM,STATUS                                         
      CHARACTER*80  FNAME ,errfn,fname1,fname2
      character*2   fnameupdn,fnum2
      character*1   fnum1                                        
!                                                                       
      real*4,allocatable :: RNUMB(:,:), DENSTY(:,:)
      real*4  fermden(MG)
      COMMON /NCOUNT/ NNOC                                              
      real*4,pointer :: EBS(:,:), FC(:,:,:)
      CHARACTER *6  DOSTYP(MG)                                          
      DIMENSION IDOS(MG,2)
      real*4,allocatable :: QTL(:,:),xqtl(:,:)  
      real*4,allocatable :: ehelp(:),denbr(:,:)
!                                                                       
      DATA PI /3.141592654/, EEF /0.0/                                      
!--------------------------------------------------------------------   
!      call getarg(2,fname)
      CALL GTFNAM(fname,ERRFN)
!      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=8002)
         if(iunit.eq.7) then
           do i=80,5,-1
             if(fname(i:i).eq.'1') then
               nnchar=i-1
               fname1(1:i-1)=fname(1:i-1)
               fnameupdn=''
               GOTO 8003
             else if(fname(i:i).eq.'p') then
               nnchar=i-3
               fname1(1:i-3)=fname(1:i-3)
               fnameupdn='up'
               GOTO 8003
             else if(fname(i:i).eq.'n') then
               nnchar=i-3
               fname1(1:i-3)=fname(1:i-3)
               fnameupdn='dn'
               GOTO 8003
             endif
           enddo
         endif
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING TETRA.DEF !!!!'
      STOP 'TETRA.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS,'  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)
!---------------------------------------------------------------------- 
      READ(5,1) TITLE
      IAV=0
      NPRINT=1
!     MI=1000                                                           
      eef=999. 
!      EEF=0.  
!     NST          number of k points                                   
!     IAV          0 NO PLOT OTHERWISE DOS AVERAGED OVER IAV VALUES     
!     NPRINT       0 NO DOS PRINTED                                     
!                  1 DOS AND INTEGRATED DOS PRINTED                     
!     NYMIN        LOWER BAND INDEX                                     
!     NYMAX        UPPER BAND INDEX                                     
!     ZEL          NUMBER OF ELECTRONS                                  
!     EMIN,DE,EMAX ENERGY GRID FOR DOS                                  
!                                                                       
!     DOS1  DOS: GES, S, P, D, EG, T2G, F, FOR ATOM NR.1                
!     DOS2  DOS: S, ... F, OUT FOR ATOM NR.2                            
!                                                                       
!*****   FIRST PRINTED ENERGY ON DOS PLOT IS                            
!*****   E0 + (IAV-1) * DE / 2 + 4 * DE * IAV                           
!*****   E.G.IAV=5    FIRST ENERGY  E0 + 22*DE                          
!*****                NEXT  AT 25*DE HIGHER                             
!                                                                       
!      read(5,*) nst
      read(3,*) nst
      rewind 3
!      if(nst.gt.nktot) then
!        write(*,*) ' nktot',nktot,' .lt. nst=',nst
!        stop
!      end if
!      READ(5,*) NYMIN,NYMAX
!      read(5,*) ZEL
      nymin=1
!      nymax=nbtot
      prev=0.d0
!     prev: Gaussian broadening parameter in Ry!!!
      read(5,*,ERR=4848) EMIN,DE,EMAX,prev
 4848 continue
!      if(nymax.gt.nbtot) then
!        write(*,*) ' nbtot',nbtot,' .lt. nymax=',nymax
!        stop
!      end if
      READ(5,*) NNDOS 
      if(nndos.gt.mg) stop 'mg too small'                                    
      WRITE(6,3) TITLE                                           
      WRITE(6,4) IAV,NPRINT                                        
!      WRITE(6,5) NYMIN,NYMAX                              
      WRITE(6,7) NNDOS                                                  
      DO 90 JCASE=1,NNDOS                                               
         READ(5,*) IATOM,ICOLUM,DOSTYP(JCASE)                          
         IDOS(JCASE,1)=IATOM                                            
         IDOS(JCASE,2)=ICOLUM                                           
         IF(IATOM.EQ.0) THEN                                            
            IDOS(JCASE,2)=0                                             
         ENDIF                                                          
         WRITE(6,25) JCASE,IDOS(JCASE,1),IDOS(JCASE,2),DOSTYP(JCASE)    
   90 CONTINUE                                                          
!                                                                       
!c      DO 200 NDOS=1,NNDOS                                               
!                                                                       
!c         IAT=IDOS(NDOS,1)                                               
!c         JAT=IDOS(NDOS,2)                                               
         REWIND 4                                                       
            READ(4,1000) SYSTEM                                         
            READ(4,1010) ALAT1,ALAT2,ALAT3,EFERM                        
            READ(4,1020) MINWAV,MAXWAV,KSPIN,NSORT,kso
        allocate( TEXT(Nsort+1),QTL(Nsort+1,15),xqtl(nsort,lxdos2*(lxdos2+1)) )
!           INCREASE NSORT BY ONE FOR THE INTERSTITIAL DOS (OUT)        
            NSORT=NSORT+1                                               
!cc         IF(NDOS.GT.1) GOTO 106                                         
            WRITE(6,2000) SYSTEM                                        
            WRITE(6,2010) ALAT1,ALAT2,ALAT3,EFERM                       
            WRITE(6,2020) MINWAV,MAXWAV,KSPIN,NSORT
!         if(kspin.eq.2.and.eferm.gt.1.5d0) then
!            write(6,*) ' Spin-polarized case dedected!!'
!            write(6,*) ' Fermi-energy must be set properly in QTL-file'
!            STOP ' FERMI must be set in QTL-file!!'
!         endif                     
         eef=eferm                                                      
         SPIN=1.d0                                                        
         IF(KSPIN.EQ.2) SPIN=0.5d0                                        
         so=1.d0
         if(kso.eq.2) then
               so=2.d0
               spin=spin/2.d0
         endif
! in spin.pol SO case set total DOS to DOS of interstital (=normso)
         if(kso.eq.2.and.kspin.eq.2) then 
          DO  JCASE=1,NNDOS                                               
           if(IDOS(JCASE,1).eq.0) then
             IDOS(JCASE,1)=NSORT                                          
             IDOS(JCASE,2)=1                                           
           ENDIF                                                          
          enddo
         endif
!
 106     KAT=NST*NSORT                                                  
!                                                                       
         DO 107 ISORT=1,NSORT                                           
!              NATO=NSORT-1  BECAUSE OF THE INTERSTITIAL DOS            
               IF(ISORT.EQ.NSORT) GOTO 107                              
               READ(4,1030) JATOM,NEQU,isplit,text(isort)                      
!cc            IF(NDOS.GT.1) GOTO 107                                      
               WRITE(6,2030) JATOM,NEQU,isplit,text(isort)                    
 107     CONTINUE  
         text(nsort)='                       '                               
!         IF(NYMIN.LE.1) GOTO 115                                        
!         NLOW=NYMIN-1                                                   
!         DO 112 NY=1,NLOW                                               
!            READ(4,1040) NRBAND                     
!            DO 112 I=1,KAT                                              
!            READ(4,21) EVAL                                             
! 112     CONTINUE                                                       
!                                                                       
         emax1=-999.d0
         nbtot=1000
         allocate (EBS(Nst,NBTOT), FC(MG,Nst,NBTOT))
         ny=0
 115     NY=NY+1
         if(ny.gt.nbtot) then
            nbtot=nbtot*2
            call doreallocate(ebs, nst, nbtot)                                
            call doreallocate(fc, mg, nst, nbtot)     
         endif                        
            READ(4,1040,end=141) NRBAND                     
            BMIN=10000.                                                 
            BMAX=-10000.                                                
            DO 120 I=1,NST                                              
               DO 113 ISORT=1,NSORT                                     
                  READ(4,1050) EVAL,ISRT,( QTL(ISORT,J), J=1,15 )
!                  WRITE(6,1050) EVAL,ISRT,( QTL(ISORT,J), J=1,13 )
                  if(text(isort)(13:20).eq.'xdos(i,j') then
                    if(lxdos.ne.3) stop 'LXDOS must be 3 for ISPLIT 99 or 88'
                  read(4,208)(xqtl(isort,i1),i1=1,(lxdos2)*(lxdos2+1)) 
                  else if(text(isort)(13:20).eq.'xdos(i,i') then
                  read(4,208)(xqtl(isort,i1),i1=1,lxdos2) 
                  endif
 113           CONTINUE            
               EBS(I,NY)=EVAL                                           
               BMIN=AMIN1(EVAL,BMIN)                                    
               BMAX=AMAX1(EVAL,BMAX)                                    
               DO 201 NDOS=1,NNDOS                                        
               IAT=IDOS(NDOS,1)                                               
               JAT=IDOS(NDOS,2)                                               
               if(iat.gt.NSORT) stop 'iat in input too large'
               IF(IAT.EQ.0) FC(NDOS,I,NY)=1.d0*SPIN                            
!               IF(IAT.GT.0.and.jat.le.14)then
               IF(IAT.GT.0.and.jat.le.15)then
                  FC(NDOS,I,NY)=QTL(IAT,JAT)*SPIN*so    
               else if(IAT.GT.0) then
                  FC(NDOS,I,NY)=xQTL(IAT,JAT-100)*SPIN
               endif
 201           continue
 120        CONTINUE                                                    
!                                                                       
            WRITE(6,14) NRBAND,BMIN,BMAX
            emax1=amax1(bmin,emax1)   
 140     goto 115
!            emax1=amin1(emax1,emax)
!         goto 200
 141     nymax=ny-1
!cc            IF(NDOS.EQ.1) then
               if(emax.gt.emax1) then
            write(6,*) ' EMAX reduced due to lower HIGHEST BAND-minimum'
            emax=emax1
               endif
!cc            endif
 200  CONTINUE                                                          
!
      JEND=1 + (EMAX-EMIN)/DE    
      allocate (RNUMB(jend,MG), DENSTY(jend,MG))
      allocate (ehelp(jend),denbr(jend,mg))
      DENSTY=0.0
      RNUMB=0.0        
      write(6,55) EMIN,DE,EMAX 
      CALL ARBDOS (0,NNDOS,1,JEND,EMIN,EMAX,NYMIN,NYMAX,DE,ebs,fc,nst,rnumb,densty)
!                                                                       
!.....WRITE OUT RESULTS TO DOS FILES                                    
      IFILE=6  
      ilmin=1  
      if(idos(1,1).eq.0) then
      do 295 j=1,jend
      if(emin+(j-1)*de.ge.eef) then
      zel=rnumb(j,1)
      goto 300
      end if
 295  continue  
      end if                                                            
 300  CONTINUE                                                          
         IFILE=IFILE+1
           if(ifile.gt.7) then
            close(ifile-1)
            close(ifile-1+50)
            fname=''
            fname2=''
            if(ifile-6.lt.10) then
                 write(fnum1,'(i1)') ifile-6
                 fname=fname1(1:nnchar)//fnum1//fnameupdn 
                 fname2=fname1(1:nnchar)//fnum1//'ev'//fnameupdn 
            else 
                 write(fnum2,'(i2)') ifile-6
                 fname=fname1(1:nnchar)//fnum2//fnameupdn
                 fname2=fname1(1:nnchar)//fnum2//'ev'//fnameupdn 
            endif 
            print*,'openfilename:',fname
            OPEN(ifile,FILE=fname,STATUS='unknown')
            OPEN(ifile+50,FILE=fname2,STATUS='unknown')
           endif                            
         IL=MIN0(7,NNDOS)  
         ilmax=il+ilmin-1                                               
         WRITE(IFILE,15) TITLE,EEF,IL,JEND,prev                              
         WRITE(IFILE,217) (DOSTYP(IX),IX=ilmin,ILmax)         
         WRITE(IFILE+50,15) TITLE,0.0,IL,JEND,prev                           
         WRITE(IFILE+50,217) (DOSTYP(IX),IX=ilmin,ILmax)      
         WRITE(6,15) TITLE,EEF,IL,JEND,prev                                  
         if(idos(1,1).eq.0) then
             write(6,44) zel
             if(spin.eq.1.d0) then
                 write(6,*) 'DOS in states/Ry'   
             else
                 write(6,*) 'DOS in states/Ry/spin'   
             endif
         endif     
         WRITE(6,117) (IDOS(IX,1),DOSTYP(IX),IX=ilmin,ILmax)
!
         DO  J=1,JEND                                                
         Ehelp(j)=EMIN+(J-1)*DE
         enddo                                             
         do IND=ILMIN,ILMAX
         call Broad(prev,jend,ehelp,DENSTY(1,IND),denbr(1,ind))
         enddo
!            
         DO J=1,JEND                                                
            ENRG=EMIN+(J-1)*DE                                             
            WRITE(6,16) ENRG,(DENSTY(J,IND)*2.d0,RNUMB(J,IND),IND=ILMIN,ILMAX) 
            WRITE(IFILE,116) ENRG,(DENbr(J,IND)*2.d0,IND=ILMIN,ILMAX)          
            WRITE(IFILE+50,116) (ENRG-eef)*13.6058,(DENbr(J,IND)/13.6058*2.d0,IND=ILMIN,ILMAX)                 
            IF(ENRG.GE.eferm.AND.enrg-DE.LE.eferm) THEN
               DO ind=ilmin,ilmax
                  fermden(ind)=2.d0*((eferm-enrg+de)*DENSTY(J,IND)+(enrg-eferm)*DENSTY(J-1,IND))/de
               ENDDO
            ENDIF
            ENDDO            
 240     CONTINUE                                                       
         NNDOS=NNDOS-7   
         WRITE(6,*) '******** EF and DOS at fermi level *******'
         WRITE(6,'(f9.5,7(F9.2,9x))') eferm,(fermDEN(ind),ind=ILMIN,ilmax)
! AM (2.80)
         DO ind=ilmin,ilmax
!            fermDEN(ind)=fermDEN(ind)/2.d0*(8.617D-5)**2*(pi**2)/3.d0*2.6255D6/27.212d0*1000.d0
            fermDEN(ind)=fermDEN(ind)*(pi**2)/3.d0*(6.3336303d-6)**2*6.0221d23*13.6058d0*1.60219d-19*1000.d0
         ENDDO
         WRITE(6,*) 'Gamma in mJ/(mol cell K**2). (Divide by number of formula units in cell to get it per mole only)'
         WRITE(6,'(A9,7(F9.2,9x))') 'Cv/T     ',(fermDEN(ind),ind=ILMIN,ilmax)

         ilmin=ilmin+7                                                  
      IF(NNDOS.GT.0) GOTO 300                                           

!DOS at Fermi
      STOP ' LEGAL END TETRA'                                           
!                                                                       
    1 FORMAT(A20)                                              
!    1 FORMAT(A20/A4,/I2/I2)                                            
    2 FORMAT(2I3/F10.5/3F10.5)                                          
    3 FORMAT(1X,A20/)    
    4 FORMAT(1X,'IAV',25X,':',I3,/,1X,'NPRINT',22X,':',I3)
 44   Format(1X,'NUMBER OF ELECTRONS UP TO EF',9X,':',F10.4,/)
    5 FORMAT(1X,'LOWER AND UPPER BAND-INDEX',2X,':',2I5)                
 55   format(1X,'EMIN, DE, EMAX',14X,':',3F10.5,/)                    
    6 FORMAT(I5)                                                        
    7 FORMAT(1X,I2,' CASES FOR DOS',12X,':  ATOM   L')                  
   13 FORMAT(8X,I2)                                                     
   14 FORMAT(' BAND LIMITS OF BAND',I3,' ARE',2F10.5)                   
   15 FORMAT('# ',A20,/,'#EF=',F10.5,5X,'NDOS=',I2,5X,'NENRG=',I5,'    Gaussian bradening:',f8.5)
   16 FORMAT(f9.5,7(F9.2,f9.4))                                      
  116 FORMAT(f10.5,7f14.8)                                              
  117 FORMAT('# ENERGY',1X,7(4x,I2,1X,A6,6X))                          
  217 FORMAT('# ENERGY',2X,7(7X,A6,1X))                             
   18 FORMAT(' ATOM NO.',I3,4X,A6,'-CHARGE',F10.5)                      
   19 FORMAT(A72/17X,3F9.5,18X,F7.4/5X,I3,9X,I2,11X,I3)                 
   20 FORMAT(A7,F4.1,1X,A9,F8.5,A9,I3)                                  
  221 FORMAT(F10.5,I3,1X,F9.5,3X,6F9.5)                                 
   21 FORMAT(F10.5,I4,1X,F10.5,5X,12F10.5)                              
   22 FORMAT(/,1X,A72/17X,3F9.5,18X,F7.4/5X,I3,9X,I2,11X,I3)            
   23 FORMAT(1X,A7,F4.1,1X,A9,F8.5,A9,I3)                               
   24 FORMAT(3X,I2,2X,I3,3X,A6)                                         
   25 FORMAT(1X,'CASE',I2,' :   ATOM NUMBER',I3,'   COLUMN READ',I3,'   DOSTYPE=',A6)
  208 FORMAT((10F10.5))                             
 1000 FORMAT(1X,A70,/)                                                  
 1010 FORMAT(1X,15X,3F8.5,16X,F10.5)                                    
 1020 FORMAT(1X,I4,10X,I4,8X,I1,7X,I3,8x,i2)                                  
 1030 FORMAT(6X,I3,7X,I2,9x,i2,2x,a45)                                       
 1040 FORMAT(7X,I3)                                                     
 1050 FORMAT(F10.5,I3,F8.5,3X,14F8.5)                                
 2000 FORMAT(/,1X,A70)                                                  
 2010 FORMAT(1X,'LATTICE CONST.=',3F8.5,3X,'FERMI ENERGY=',F10.5)       
 2020 FORMAT(1X,I4,' < NMAT < ',I4,3X,'SPIN=',I1,3X,'NATO=',I2)         
 2030 FORMAT(1X,'JATOM',I3,2X,'MULT=',I2,2x,'ISPLIT=',i2,2x,a45)            
      END                                                               
