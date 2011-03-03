      SUBROUTINE SCFANA(NAT,JSPIN,FORSUM,ESUM,force,iscf,cut,itertwice,mult,fninc)
!
!     read scf file and sums up total energy and forces
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!
      CHARACTER*79   MARGN                                             
      character*67   ERRMSG
      CHARACTER*80   fninc,fninc1
!
!        Common blocks
!
!                                                                       
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        VI          - volume of the direct lattice unit-cell
!
      INTEGER            MULT(NAT),iscftest(999)
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3)
      COMMON  /STRUK/    ALAT, ALPHA, PIA, VI
      SAVE    /STRUK/
!
!
!
      LOGICAL   FORCE(NAT+1),cut,itertwice
      DIMENSION FORSUM(NAT,4),FHELP(4)       
      real*8,allocatable ::  ESX(:,:)       
!
! parameter for the true occupation read from case.inc
!
      INTEGER            IDUMMY,jdummy,index(1:6,-4:3)
      INTEGER,allocatable ::            NORB(:)
      real*8,allocatable ::             NEL(:,:)
      REAL*8,allocatable ::             OCC(:,:,:)
!
!---------------------------------------------------------------------  
!                                                                       
      allocate (  ESX(26,NAT))       
      allocate ( NORB(NAT))
      allocate ( NEL(30,NAT))
      allocate ( OCC(30,NAT,2))
!
! INITIALIZE FORCE AND TOTAL ENERGY PARAMETER
!
      index(1,-1)=1
      index(2,-1)=2
      index(2, 1)=3
      index(2,-2)=4
      index(3,-1)=5
      index(3, 1)=6
      index(3,-2)=7
      index(3, 2)=8
      index(3,-3)=9
      index(4,-1)=10
      index(4, 1)=11
      index(4,-2)=12
      index(4, 2)=13
      index(4,-3)=14
      index(5,-1)=15
      index(5, 1)=16
      index(5,-2)=17
      index(4, 3)=18
      index(4,-4)=19
      index(5, 2)=20
      index(5,-3)=21
      index(6,-1)=22
      index(6, 1)=23
      index(6,-2)=24
      index(5, 3)=25
      index(5,-4)=26
!
      ESUM=0.d0                                                           
      EDENS=0.d0                                                          
      FORCE=.false.
!      FORCE(JATOM)=.false.
      FORSUM=0.d0 
      ESX=0.d0
      NEL=0.d0
      OCC=0.d0
      NORB=0
      iscftest=0
! .................................................................
!
! READ true occupation from case.inc/up/dn
!
      if(jspin.eq.2) then
       do jj=80,1,-1
         if(fninc(jj:jj).ne.' ') exit
       enddo
       fninc1=' '
       fninc1(1:jj+2)=fninc(1:jj)//'up'
       write(6,*) 'filename of case.inc: ',fninc1(1:jj+2)    
       close (7)
       OPEN (7,FILE=fninc1,STATUS='unknown',FORM='formatted',ERR=920)
       DO  JATOM=1,NAT
        READ(7,*,end=32) NORB(JATOM)
        IF(NORB(JATOM).EQ.0) cycle
          DO  J=1,NORB(JATOM)
            READ(7,*,end=32) IDUMMY,jDUMMY,NEL(J,JATOM)
            OCC(index(idummy,Jdummy),JATOM,1) = NEL(J,JATOM) 
!           WRITE(6,*) JATOM,J,NEL(J,JATOM),OCC(J,JATOM,1)
          enddo
       enddo
       fninc1(1:jj+2)=fninc(1:jj)//'dn'
       write(6,*) 'filename of case.inc: ',fninc1(1:jj+2)    
       close (7)
       OPEN (7,FILE=fninc1,STATUS='unknown',FORM='formatted',ERR=920)
       DO  JATOM=1,NAT
        READ(7,*,end=32) NORB(JATOM)
        IF(NORB(JATOM).EQ.0) cycle
          DO  J=1,NORB(JATOM)
            READ(7,*,end=32) IDUMMY,jDUMMY,NEL(J,JATOM)
            OCC(index(idummy,Jdummy),JATOM,2) = NEL(J,JATOM) 
!           WRITE(6,*) JATOM,J,NEL(J,JATOM),OCC(J,JATOM,2)
          enddo
       enddo
      goto 33
      endif
! non-spinpolarized inc file
 32   continue 
       close (7)
       OPEN (7,FILE=fninc,STATUS='old',FORM='formatted',ERR=920)
      DO 31 JATOM=1,NAT
        READ(7,*) NORB(JATOM)
        IF(NORB(JATOM).EQ.0) GOTO 31
          DO J=1,NORB(JATOM)
            READ(7,*) IDUMMY,jDUMMY,NEL(J,JATOM)
            OCC(index(idummy,Jdummy),JATOM,1) = NEL(J,JATOM) 
            OCC(index(idummy,Jdummy),JATOM,2) = NEL(J,JATOM) 
!           WRITE(6,*) JATOM,J,NEL(J,JATOM),OCC(J,JATOM,1)
          enddo
 31   CONTINUE
 33   CONTINUE
!

!                                                                       
!.....READ TAPE22=SCFDATA UNTIL BOTTOM                                  
!     
      ispin=0                                                                  
      ISCF=0                                                            
  10  READ(22,700,err=20,END=20) MARGN
      IF(MARGN(2:4).EQ.'ITE')GO TO 30                                          
!      WRITE(21,700)MARGN                                               
      GO TO 10                                                          
  30  READ(MARGN,710) ISCF 
      iscftest(iscf)=iscftest(iscf)+1
!      WRITE(21,700)MARGN
      GOTO 10                                                           
  20  CONTINUE                                                          
!
!..... READ RIGHT SCF SECTION
!
      REWIND 22                                                         
 60   READ(22,700,END=40,err=40) MARGN
      IF (MARGN(2:4).NE.'ITE') GO TO 60                                        
      READ(MARGN,710) ISCFT 
      IF (ISCFT.NE.ISCF) GO TO 60                                        
      if(iscftest(iscf).ne.1) then
         write(6,*) ':WARN : Iteration ',iscf,'occurs more than once in this scf file!!!'
         iscftest(iscf)=iscftest(iscf)-1
         goto 60
      endif
 50   READ(22,700,END=40) MARGN
!
!..... READ TOTAL ENERGY PARTS
!
      IF (MARGN(2:4).EQ.'WAR') cut=.true.
      IF (MARGN(2:4).EQ.'DEN') READ(MARGN,720)EDENS
      IF (MARGN(2:4).EQ.'SUM') THEN                                         
         READ(MARGN,720) ES1 
         ESUM=ESUM+ES1                                                     
      ENDIF                                                             
!      IF (MARGN(2:6).EQ.'LDAUE') THEN                                         
!         READ(MARGN,725) ECORR, EDC
!         ESUM = ESUM + ECORR -EDC
!              write(21,'(10x,5HECORR,3f10.6)') ECORR, -EDC,ECORR-EDC
!      ENDIF                                                             
      IF (MARGN(2:5).EQ.'EORB') THEN                                         
         READ(MARGN,725) EORB
         ESUM = ESUM + EORB
              write(6,'(10x,5HEORB ,3f10.6)') EORB 
      ENDIF                                                             
      IF (MARGN(2:4).EQ.'1S ') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX
         if(jatom.eq.1) then
             ispin=ispin+1     
!             print*, 'ispin set to',ispin 
             if(ispin.gt.2) ispin=1
         endif                                          
         ESX(1,JATOM)=ESX(1,JATOM)+ESXX*OCC(1,JATOM,ispin)*MULT(JATOM)/JSPIN 
      ELSE IF(MARGN(2:4).EQ.'2S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(2,JATOM) =ESX(2,JATOM)+ESXX*OCC(2,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'2PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(3,JATOM) =ESX(3,JATOM)+ESXX*OCC(3,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'2P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(4,JATOM) =ESX(4,JATOM)+ESXX*OCC(4,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'3S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(5,JATOM) =ESX(5,JATOM)+ESXX*OCC(5,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'3PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(6,JATOM) =ESX(6,JATOM)+ESXX*OCC(6,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'3P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(7,JATOM) =ESX(7,JATOM)+ESXX*OCC(7,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'3DD') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(8,JATOM) =ESX(8,JATOM)+ESXX*OCC(8,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'3D ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(9,JATOM) =ESX(9,JATOM)+ESXX*OCC(9,JATOM,ispin)*MULT(JATOM)/JSPIN                                 
      ELSE IF(MARGN(2:4).EQ.'4S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(10,JATOM) =ESX(10,JATOM)+ESXX*OCC(10,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'4PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(11,JATOM) =ESX(11,JATOM)+ESXX*OCC(11,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'4P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(12,JATOM) =ESX(12,JATOM)+ESXX*OCC(12,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'4DD') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(13,JATOM) =ESX(13,JATOM)+ESXX*OCC(13,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'4D ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(14,JATOM) =ESX(14,JATOM)+ESXX*OCC(14,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'5S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(15,JATOM) =ESX(15,JATOM)+ESXX*OCC(15,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'5PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(16,JATOM) =ESX(16,JATOM)+ESXX*OCC(16,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'5P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(17,JATOM) =ESX(17,JATOM)+ESXX*OCC(17,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'4FF') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(18,JATOM) =ESX(18,JATOM)+ESXX*OCC(18,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'4F ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(19,JATOM) =ESX(19,JATOM)+ESXX*OCC(19,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'5DD') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(20,JATOM) =ESX(20,JATOM)+ESXX*OCC(20,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'5D ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(21,JATOM) =ESX(21,JATOM)+ESXX*OCC(21,JATOM,ispin)*MULT(JATOM)/JSPIN                              
      ELSE IF(MARGN(2:4).EQ.'6S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(22,JATOM) =ESX(22,JATOM)+ESXX*OCC(22,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'6PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(23,JATOM) =ESX(23,JATOM)+ESXX*OCC(23,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'6P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(24,JATOM) =ESX(24,JATOM)+ESXX*OCC(24,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'5FF') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(25,JATOM) =ESX(25,JATOM)+ESXX*OCC(25,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      ELSE IF(MARGN(2:4).EQ.'5F ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(26,JATOM) =ESX(26,JATOM)+ESXX*OCC(26,JATOM,ispin)*MULT(JATOM)/JSPIN                               
      END IF                                                            
!
!..... READ FORCE PARTS
!
      IF (MARGN(2:4).EQ.'FHF') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,740) (FHELP(j),j=1,4)
         FORCE(JATOM)=.true.
         FORCE(NAT+1)=.true.
      do 200 j=2,4
         FORSUM(JATOM,j)=FORSUM(JATOM,j)+FHELP(j)
 200  CONTINUE
      END IF                                                            
      IF (MARGN(2:4).EQ.'FVA') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,740) (FHELP(j),j=1,4)
      do 210 j=2,4
         FORSUM(JATOM,j)=FORSUM(JATOM,j)+FHELP(j)                               
 210  CONTINUE
      END IF                                                            
! GM 13/8 - 00
      IF (MARGN(2:4).EQ.'FSU') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,740) (FHELP(j),j=1,4)
      do j=2,4
         FORSUM(JATOM,j)=FORSUM(JATOM,j)+FHELP(j)
      enddo
      END IF                                                            
! .................
      IF (MARGN(2:4).EQ.'FCO') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,740) (FHELP(j),j=1,4)
      do 220 j=2,4
         FORSUM(JATOM,j)=FORSUM(JATOM,j)+FHELP(j)/jspin
 220  CONTINUE
      END IF                                                            
      IF (MARGN(2:4).NE.'ITE') GO TO 50                                        
      cut=.true.
      itertwice=.true.
      write(6,888) iscf
 888  format(':WARN : Iteration ',i3,'occurs more than once in this scf file!!!', &
      /,7x,' ETOT is not calculated correctly !!!!!', &
      /,7x,' You should save your calculations using save_lapw and continue')
      GO TO 50                                                          
 40   CONTINUE                                                          
! DEN ENERGIE DE L ORBITALE EN UNITE ATOMIQUE ET NEGATIVE
! NQN NOMBRE QUANTIQUE PRINCIPAL   NK NOMBRE QUANTIQUE KAPPA            
! NEL OCCUPATION DE L ORBITALE                      
!
! ..................................................................
!
! CALCULATE TOTAL ENERGY
!
      ESUM=ESUM+EDENS                                                   
      DO 11 JATOM=1,NAT                                                 
      FORSUM(JATOM,1)=sqrt &
            (forsum(jatom,2)**2 &
            +forsum(jatom,3)**2 &
            +forsum(jatom,4)**2)
      DO 9 J=1,26 
      if((occ(j,jatom,1).gt.0.d0).or.(occ(j,jatom,2).gt.0.d0)) then
        WRITE(6,*) 'JATOM = ',JATOM,'  J = ',J
        WRITE(6,*) '  ESX = ',ESX(J,JATOM),'  OCC = ',OCC(J,JATOM,1),OCC(J,JATOM,2)
        WRITE(6,*) '  ESUM = ',ESUM
      endif
      ESUM=ESUM+ESX(J,JATOM)                      
!      ESUM=ESUM+ESX(J,JATOM)*OCC(J,JATOM)*MULT(JATOM)/JSPIN                       
 9    continue
 11   CONTINUE
      WRITE(6,*) '  ESUM = ',ESUM
!                                                                       
      RETURN
 920  INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) 7
      CALL OUTERR('MIXER',ERRMSG)
      WRITE (ERRMSG,9020) FNinc//' in scfana'
      CALL OUTERR('MIXER',ERRMSG)
      STOP 'MIXER - Error'
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 
!
!
!
 700  FORMAT(A79)                                                    
 710  FORMAT(4X,i3)                                              
 720  FORMAT(30X,f20.6)                                              
 725  FORMAT(7X,2f13.6)                                              
 730  FORMAT(21X,f20.6)                                             
 740  FORMAT(17X,4f15.3)                                             
!
      END 

