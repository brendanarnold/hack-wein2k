      SUBROUTINE SCFSUM(SCFFN,SCFALL,nat)
!
!     read scf file and sums up total energy and forces
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      CHARACTER*100   MARGN                                             
!
!        Common blocks
!
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
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3)
      COMMON  /STRUK/    ALAT, ALPHA,  PIA, VI
      SAVE    /STRUK/
!
!
!
      LOGICAL,allocatable ::   FORCE(:)    !nato
      real*8,allocatable :: QTL(:,:)       !NATO,12)
      real*8,allocatable :: epL(:,:),FORSUM(:,:),FSUSUM(:,:),QpL(:,:)  !nato,4
      real*8,allocatable :: eph(:,:),Qph(:,:),VZZ(:,:)   !NATO,4)
      real*8,allocatable :: CHA(:)      ! NATO)
      DIMENSION QTLIN(12),epLIN(8),FHELP(4),ephIN(8),EFG(4)
 
      CHARACTER*80 FNAME,SCFFN,SCFALL
      COMMON /IPROC/   IPROC
!---------------------------------------------------------------------  
      allocate ( FORCE(nat))
      allocate ( QTL(NAT,12))
      allocate ( epL(nat,4),FORSUM(nat,4),FSUSUM(nat,4),QpL(nat,4))
      allocate ( eph(nat,4),Qph(nat,4),VZZ(nat,4))
      allocate ( CHA(NAT))

!      WRITE(6,*)'In SCFSUM'
!                                                                       
!
! INITIALIZE 
!
      EVALS=0.d0
      QTLB =0.d0
      CHAC =0.d0
      CHACS=0.d0
      ESUM =0.d0                                                           
      EDENS=0.d0
      do j=1,8
         eplin(j)=0.d0
         ephin(j)=0.d0
      enddo                                                        
      DO 101 JATOM=1,NAT
      do 101 j=1,4
         epl(jatom,j)=0.d0
         qpl(jatom,j)=0.d0
         eph(jatom,j)=0.d0
         qph(jatom,j)=0.d0
         FORSUM(JATOM,j)=0.d0
         FSUSUM(JATOM,j)=0.d0
         FHELP(J)=0.d0
         FORCE(JATOM)=.false.
         EFG(J)=0
 101  CONTINUE
      DO i=1,NAT
         DO N=1,12
            QTL(i,n)=0.d0
            QTLIN(n)=0.d0
         ENDDO
         DO N=1,4
            VZZ(I,N)=0.d0
         ENDDO
         CHA(I)=0.d0
      ENDDO
  100 CONTINUE
!....read bandenergy correction
      ts2=0.d0
 9    READ(21,700,END=8) MARGN
        IF(MARGN(10:19).EQ.'BANDENERGY') THEN
            READ(MARGN,'(31x,f10.6)') ts2
            write(*,*) '      BANDENERGY CORRECTION:',ts2
            goto 8
        endif
        IF(MARGN(11:18).EQ.'-(T*S)/2') THEN
            READ(MARGN,'(30x,f10.6)') ts2
            write(*,*) '      BANDENERGY CORRECTION:',ts2
            goto 8
        endif
      goto 9
 8    continue

!                                                                       
!.....READ TAPE21=SCFDATA UNTIL BOTTOM                                  
!
!      close(22)
!      OPEN(22,FILE=SCFALL,STATUS='unknown',FORM='formatted')

      DO ILOOP=1,IPROC
         WRITE(6,*)'SCF: DOING PROCESSOR FILE: ',ILOOP
         CALL mknam(FNAME,SCFFN,ILOOP)
         CLOSE(21)
!         write(6,*)'opening ',FNAME
         OPEN(21,FILE=FNAME,STATUS='old',FORM='formatted')
 10      READ(21,700,END=20) MARGN
         write(6,*)MARGN
        IF(MARGN(2:4).EQ.'CHA') THEN
            IF(MARGN(2:6).EQ.'CHA  ')THEN
               READ(MARGN,881)CHAC
               CHACS=CHACS+CHAC
               IF(ILOOP.EQ.IPROC)THEN
                  write(22,880)CHACS
               ENDIF
            ELSE
               READ(MARGN,721)IATOM, IATOM,X
               CHA(IATOM)=CHA(IATOM)+X
               IF(ILOOP.EQ.IPROC)THEN
                  WRITE(22,720)IATOM,IATOM,CHA(IATOM)
               ENDIF
            ENDIF
         ELSE IF(MARGN(2:4).EQ.'QTL') THEN
            READ(MARGN,251)IATOM,(QTLIN(I),I=1,12)
            DO i=1,12
               QTL(IATOM,I)=QTL(IATOM,I)+QTLIN(I)
            ENDDO
            IF(ILOOP.EQ.IPROC)THEN
               WRITE(22,250)IATOM,(QTL(IATOM,I),I=1,12)
            ENDIF
         ELSE IF(MARGN(2:4).EQ.'EPL') THEN
            READ(MARGN,260)IATOM,(eplin(I),I=1,8)
            DO i=1,4
               QpL(IATOM,I)=QpL(IATOM,I)+epLIN(I*2-1)
               epL(IATOM,I)=epL(IATOM,I)+eplin(i*2)*epLIN(I*2-1)
            ENDDO
            IF(ILOOP.EQ.IPROC)THEN
               do i=1,4
!RED   Abfrage aus SRC_lapw2/l2main hierher verschoben
!              if(qpl(iatom,i).lt.0.0001d0) then
               if(qpl(iatom,i).lt.0.00005d0) then
!RED   Abfrage aus SRC_lapw2/l2main hierher verschoben
                 qpl(iatom,i)=0.00001d0
                 epl(iatom,i)=10.d0*qpl(iatom,i)
               endif
               enddo
               WRITE(22,261)IATOM,(QpL(IATOM,I),epl(iatom,i)/ &
                                                QpL(IATOM,I),I=1,4)
            ENDIF
         ELSE IF(MARGN(2:4).EQ.'EPH') THEN
            READ(MARGN,260)IATOM,(ephin(I),I=1,8)
            DO i=1,4
               Qph(IATOM,I)=Qph(IATOM,I)+ephIN(I*2-1)
               eph(IATOM,I)=eph(IATOM,I)+ephin(i*2)*ephIN(I*2-1)
            ENDDO
            IF(ILOOP.EQ.IPROC)THEN
               do i=1,4
               if(qph(iatom,i).lt.0.0001d0) then
                 qph(iatom,i)=0.00001d0
                 eph(iatom,i)=10.d0*qph(iatom,i)
               endif
               enddo
               WRITE(22,263)IATOM,(Qph(IATOM,I),eph(iatom,i)/ &
                                                Qph(IATOM,I),I=1,4)
            ENDIF
         ELSE IF(MARGN(2:4).EQ.'SUM')THEN
            READ(MARGN,731)EVAL
            EVALS=EVALS+EVAL
            IF(ILOOP.EQ.IPROC)THEN
               WRITE(22,730)EVALS+ts2
            ENDIF            
         ELSE IF(MARGN(2:4).EQ.'VZZ')THEN
            IF(MARGN(52:56).EQ.'UP TO')THEN
               iform=1
            ELSE
!RED nur zu Kontrollzwecken geaendert
!              IF(MARGN(73:74).EQ.'.')THEN
               IF(MARGN(73:73).EQ.'.')THEN
!RED nur zu Kontrollzwecken geaendert
                  iform=4 
               ELSE
                  iform=3
               ENDIF
            ENDIF
            IF(IFORM.EQ.1)THEN
               READ(MARGN,7721)i,i,EFG(1),r
            ELSE IF(IFORM.EQ.3)THEN
               READ(MARGN,3021)i,EFG(1),EFG(2), &
                                 EFG(3),r              
            ELSE IF(IFORM.EQ.4)THEN
               READ(MARGN,3121)i,EFG(1),EFG(2), &
                                 EFG(3),EFG(4),r
            ENDIF
            DO L=1,iform
               VZZ(i,L)=VZZ(i,L)+EFG(L)
            ENDDO
            IF(ILOOP.EQ.IPROC)THEN
               IF(IFORM.EQ.1)THEN
                  WRITE(22,7720)i,i,VZZ(i,1),r
               ELSE IF(IFORM.EQ.3)THEN
                  WRITE(22,3020)i,VZZ(i,1),VZZ(i,2), &
                       VZZ(i,3),r              
               ELSE IF(IFORM.EQ.4)THEN
                  WRITE(22,3120)i,VZZ(i,1),VZZ(i,2), &
                       VZZ(i,3),VZZ(i,4),r
               ENDIF               
            ENDIF
                      
         ELSE IF(MARGN(2:4).EQ.'FVA')THEN
            READ(MARGN,81)IATOM,IMIST,(FHELP(I),I=1,4)
            DO i=2,4
               FORSUM(IATOM,I)=FORSUM(IATOM,I)+FHELP(I)
            ENDDO
            FORSUM(IATOM,1)=SQRT(FORSUM(IATOM,2)**2 &
                 +FORSUM(IATOM,3)**2+FORSUM(IATOM,4)**2)
            IF(ILOOP.EQ.IPROC)THEN
               WRITE(22,80)IATOM,IATOM,(FORSUM(IATOM,I),I=1,4)
            ENDIF
         ELSE IF(MARGN(2:4).EQ.'FSU')THEN
            READ(MARGN,81)IATOM,IMIST,(FHELP(I),I=1,4)
            DO i=2,4
               FSUSUM(IATOM,I)=FSUSUM(IATOM,I)+FHELP(I)
            ENDDO
            FSUSUM(IATOM,1)=SQRT(FSUSUM(IATOM,2)**2 &
                 +FSUSUM(IATOM,3)**2+FSUSUM(IATOM,4)**2)
            IF(ILOOP.EQ.IPROC)THEN
               WRITE(22,78)IATOM,IATOM,(FSUSUM(IATOM,I),I=1,4)
            ENDIF
         ELSE IF(MARGN(4:8).EQ.'QTL-B') THEN
            READ(MARGN,1001)HELP
            IF(HELP.GT.QTLB) QTLB=HELP
            WRITE(*,*)'QTL-B=',QTLB
         ELSE
            IF(ILOOP.EQ.IPROC)THEN
!               call strln(MARGN,LEN)
               len=len_trim(margn)
               IF(LEN.LT.1) LEN=1 
               WRITE(22,'(a)')MARGN(1:LEN)
            ENDIF
         ENDIF
         GOTO 10
 20      CONTINUE
      ENDDO
      IF(QTLB.NE.0)THEN
         WRITE(22,'(//)')
         WRITE(22,1000)QTLB
      ENDIF

                  
!                                                                       
      RETURN
!
!
!
                                            
 78   FORMAT(':FSU',i3.3,':',1x,i3,'.ATOM',4f15.3)
  80  FORMAT(':FVA',i3.3,':',1x,i3,'.ATOM',4f15.3)
  81  FORMAT(4x,i3.3,2x,i3,5x,4f15.3)
 250  FORMAT(':QTL',i3.3,':',12F7.4)               
 251  FORMAT(4x,i3.3,1x,12F7.4)
 261  FORMAT(':EPL',i3.3,':',4(2F8.4,2x))
!RED   Erhöhung der Genauigkeit wie in SRC_lapw2/l2main.F
!260  FORMAT(4x,i3.3,1x,4(2F8.4,2x))
 260  FORMAT(4x,i3.3,1x,4(2F11.7,2x))
!RED   Erhöhung der Genauigkeit wie in SRC_lapw2/l2main.F
 263  FORMAT(':EPH',i3.3,':',4(2F8.4,2x))
 700  FORMAT(A100)        
 720  FORMAT(':CHA',i3.3,':',1X,'TOTAL CHARGE INSIDE SPHERE ',I3, &
           ' = ',F12.6)         
 721  FORMAT(4x,i3.3,2X,27x,I3,3x,F12.6)
 730  FORMAT(':SUM  :',1X,'SUM OF EIGENVALUES =  ',F20.6)
 731  FORMAT(7X,1X,22x,F20.6)
 880  FORMAT(':CHA  :',' TOTAL CHARGE INSIDE UNIT CELL =',F15.6) 
 881  FORMAT(7x,32x,F15.6)
 1000 FORMAT('   QTL-B VALUE .EQ. ',F10.5,'   !!!!!!')
 1001 FORMAT(20X,F10.5)
 7720 FORMAT(':VZZ',i3.3,':',1X,'EFG INSIDE SPHERE ',I3,' = ',F12.6,5X, &
      ' UP TO R =',F10.5)
 7721 FORMAT(4x,i3.3,2X,18x,I3,3x,F12.6,5X,10x,F10.5)
 3120 FORMAT(':VZZ',i3.3,':',8x,4F12.5,F12.3)
 3121 FORMAT(4x,i3.3,1x,8x,4F12.5,F12.3)
 3020 FORMAT(':VZZ',i3.3,':',8x,3F12.5,F12.3)
 3021 FORMAT(4x,i3.3,1x,8x,3F12.5,F12.3)
                       
!     
      END 
