!$hp9000_800 intrinsics on
      PROGRAM clmcopy                                                     
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      CHARACTER*19     nkktext
      CHARACTER*4    IREL,LATTIC,cform                                 
      CHARACTER*5    vech(20,3)                                                
      CHARACTER*10   ANAME                                              
      CHARACTER*11   STATUS,FORM                                        
      CHARACTER*20   TITLE                                              
      CHARACTER*77   FNAME
      CHARACTER*79   MARGN
      CHARACTER*92   MARGN92
      LOGICAL        jatch2,forces,posprint
      LOGICAL        kdelta
      external       kdelta
!
      integer,allocatable :: jatch(:,:),lmch(:,:,:,:),lmchm(:)             
      COMPLEX*16          TAUP(nsym)
      complex*16 tk,imag                                
      COMPLEX*16,allocatable ::  ROKVL(:,:)                
      INTEGer,allocatable :: KVECVL(:,:)
      real*8,allocatable :: CLMNEW(:,:,:,:),fhelp(:,:),fsuhelp(:,:),clmfac(:,:)
      real*8,allocatable :: symmat(:,:,:)
      INTEGer,allocatable :: LM(:,:,:),LMMAX(:),JRI(:),MULT(:)
!
      integer msym(3,3),krot(3),krotst(3,nsym)
      real*8  taum(3)
      DIMENSION fakc(20,2)
      COMMON /SYMETR/  TAU(3,NSYM),lsym(3,3,NSYM),             &
                       INUM(NSYM),IORD                                  
!      COMMON /SYMM/  LSYM(3,3,NSYM),IORD,INUM(NSYM)                     
      DATA      tiny/1.0d-15/
!---------------------------------------------------------------------  
!                                                                       
      call getarg(2,fname)
      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING MIXER.DEF !!!!'
      STOP 'MIXER.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)

      tpi=2.d0*acos(-1.d0)
      imag=(0.d0,1.d0)

!.....READ TAPE20=POTE                                                  
!                                                                       
      READ(20,1000) TITLE                                               
      READ(20,1010) LATTIC,NAT,cform,IREL              
      nato=nat                       
      allocate ( CLMNEW(NRAD,NCOM,NATO,2),fhelp(nato,4),fsuhelp(nato,4))
      allocate ( LM(2,NCOM,NATO),LMMAX(NATO),JRI(NATO),MULT(NATO))
      allocate ( symmat(3,3,nato))
      allocate ( jatch(nato,2),lmch(nato,ncom,2,2),lmchm(nato),clmfac(nato,ncom))                                              
      READ(20,1030) AA,BB,CC                                            
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
         READ(20,1040) IATNR,POS1,POS2,             &
                       POS3,MULT(JATOM),ISPLIT            
            MULTND=MULT(JATOM)-1                                        
            DO 10 M=1,MULTND                                            
            INDEX=INDEX+1                                               
            READ(20,1041) IATNR,POS1,POS2,          &
                          POS3                                   
 10         CONTINUE                                                    
         READ(20,1050) ANAME,JRI(JATOM),RNOT,RMTT,         &
                       ZZ                                        
         READ(20,*)                                                     
         READ(20,*)                                                     
         READ(20,*)                                                     
 5    CONTINUE                                                          
!                                                                       
!     READ SYMMETRY-OPERATIONS FROM TAPE20=POTE                         
      READ(20,1055) IORD                                                
      DO 12 J=1,IORD                                                    
 12   READ(20,1056) ( (LSYM(J1,J2,J),J1=1,3),TAU(J2,J), J2=1,3 ),INUM(J)
!     ALL OF POTE IS READ                                               
!                                                                       
!.....LOOK, WHICH CHARGE-DENSITY TAPES EXIST                            
!                                                                       
      READ(17,2032,IOSTAT=ITAPE)                                        
      write(18,2032,IOSTAT=ITAPE)                                        
!                                                                       
      DO 15 JATOM=1,NAT                                                 
            READ(17,1980)                                               
            READ(17,2000) LMMAX(JATOM)                                  
               if(lmmax(jatom).gt.ncom) stop 'ncom too small'
            DO 20 LM1=1,LMMAX(JATOM)                                    
            READ(17,2010) LM(1,LM1,JATOM),LM(2,LM1,JATOM)               
            READ(17,iform1) (CLMNEW(J,LM1,JATOM,1),J=1,JRI(JATOM))
 20         READ(17,2031)                                               
            READ(17,2033)                                               
 15      continue
! reading input parameters : changing atoms
      read (5,*) jmax
      do 30 j=1,jmax
        read(5,*) (jatch(j,jx),jx=1,2) 
        read(5,*) (symmat(1,i,j),i=1,3)
        read(5,*) (symmat(2,i,j),i=1,3)
        read(5,*) (symmat(3,i,j),i=1,3)
        read(5,*) lmchm (j)
           do iclm=1,lmchm (j)
!                    old-L,M   new-L,M  factor
           read(5,*) lmch(j,iclm,1,1),lmch(j,iclm,2,1), &
                     lmch(j,iclm,1,2),lmch(j,iclm,2,2),clmfac(j,iclm) 
           enddo
 30   continue 
      do 32 j=1,jmax
        write(6,3000) j,(jatch(j,jx),jx=1,2) 
        write(6,3001) lmchm (j)
        write(6,3002) (lmch (j,jx,1,1),lmch (j,jx,2,1),jx=1,lmchm(j)) 
        write(6,3003) (lmch (j,jx,1,2),lmch (j,jx,2,2),jx=1,lmchm(j)) 
        write(6,3004) (clmfac (j,jx),jx=1,lmchm(j)) 
 32   continue 
 3000 format(i2,' case: atom',i3,'   interchanged with',i3)
 3001 format(3x,i2,' LM values to change')
 3002 format(3x,'orig LM',24i3)
 3003 format(3x,'new  LM',24i3)
 3004 format(3x,'factor ',12f6.1)
!
! calculate clmfiles
!
      DO 16 JATOM=1,NAT    
        jato=jatom
        do 35 j=1,jmax
          if(jatch(j,1).eq.jatom ) then
               jx=1
               jato=jatch(j,2)
!          endif              
          else if(jatch(j,2).eq.jatom ) then
               jx=2
               jato=jatch(j,1)
          endif       
 35     continue                                             
!        print *, 'atom',jatom,' generated from atom',jato
        write(18,1990) jatom                                              
        write(18,2001) LMMAX(JATO)                                  
        DO 21 LM1=1,LMMAX(JATO)                                    
          write(18,2011) LM(1,LM1,JATO),LM(2,LM1,JATO)               
!
          lm2=lm1
          FFF=1.0D0
          do 40 j=1,jmax
          if ((JATO.EQ.jatch(j,1)).or.(jato.eq.jatch(j,2))) then
            do 45 lj=1,lmchm (j)
            IF((LM(1,LM1,JATO).EQ.lmch(j,lj,1,1)).and. &
               (LM(2,LM1,JATO).EQ.lmch(j,lj,2,1))) then
                  do lm2=1,lmmax(jato)
                  IF((LM(1,LM2,JATO).EQ.lmch(j,lj,1,2)).and. &
                     (LM(2,LM2,JATO).EQ.lmch(j,lj,2,2))) then
                     fff=clmfac(j,lj)
!                    print*,'case1',lj,lm2,fff
                     goto 41
                  endif
                  enddo
            else if((LM(1,LM1,JATO).EQ.lmch(j,lj,1,2)).and. &
               (LM(2,LM1,JATO).EQ.lmch(j,lj,2,2))) then
                  do lm2=1,lmmax(jato)
                  IF((LM(1,LM2,JATO).EQ.lmch(j,lj,1,1)).and. &
                     (LM(2,LM2,JATO).EQ.lmch(j,lj,2,1))) then
                     fff=clmfac(j,lj)
!                    print*,'case2',lj,lm2,fff
                     goto 41
                  endif
                  enddo
            endif
 45         continue
                  lm2=lm1
                  fff=1.d0
!                    print*,'case3',lm2,fff
          endif
 40       continue           
 41       continue           
!
        if(abs(fff-1.D0).lt.1d-7)then
!                    print *, 'LM',lm1,' copied from',lm2
                    WRITE(18,IFORM1) (CLMNEW(J,LM2,JATO,1),J=1,JRI(JATO))
        else if(abs(fff+1.D0).lt.1d-7)then
!                    print *, 'LM',lm1,' negative from',lm2
                    WRITE(18,IFORM1) (-CLMNEW(J,LM2,JATO,1),J=1,JRI(JATO))
        else
!              print *, 'LM',lm1,' generated from',lm2,' with',fff
            WRITE(18,IFORM1) (FFF*CLMNEW(J,LM2,JATO,1),J=1,JRI(JATO))
        endif
            write(18,2031)                                               
 21         continue
            write(18,2033)                                               
 16   continue                      
!
         READ(17,1980)                                                  
!        READ(17,2060) NKKVL                                            
  read(17,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkkvl
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkkvl
 6768 continue
         write(18,1980)                                                  
         write(18,2061) NKKVL   
         nwav=nkkvl
      allocate (  ROKVL(NWAV,2))                
      allocate ( KVECVL(3,NWAV))

! reading input parameters : changing vectors 
      do j2=1,3
      READ(5,*) (mSYM(J1,J2),J1=1,3),TAUm(J2)
      enddo
! reading fourier coeff.
         DO  J=1,NKKVL                                               
           READ(17,iform2)  (KVECVL(JX,J),JX=1,3),ROKVL(J,1)              
         enddo       
! writing:
         DO  J=1,NKKVL  
             TKK=0.                                                          
             DO  J1=1,3
               TKK=TKK+TAUm(J1)*KVECVL(J1,J)
!              TK=TK+TAUm(J1)*KVECVL(J1,J)*TPI        
               Krot(J1)=0                                      
             DO L=1,3                                                  
               Krot(J1)=msym(J1,L)*KVECVL(L,J)+Krot(J1)                   
             ENDDO
!!             Krot(j1)=(msym(J1,1)+msym(J1,2)+msym(J1,3))*KVECVL(J1,J)
             ENDDO
!             tk=EXP(IMAG*TK)
             TKK=TKK*TPI
             ct1=cos(tkk)
             st1=sin(tkk)   
             tk=cmplx(ct1,st1)
!
     call STERN (1,NST,krotst,TAUP,krot,1)  
     if(i.lt.10) write(*,'(i4,/(3i4,2f10.5))') nst, ((krotst(i5,i6),i5=1,3),taup(i6),i6=1,nst)                             
!                                                                       
!     STERN GENERATES THE STAR OF REC LATTICE VECTOR KZZ(1,JJA)         
!     THE STAR VECTORS ARE STORED IN STG, THE STAR-SIZE IN  NST         
!     IZ CONTAINS THE SYMMETRY-MATRICES                                 
!             CALL STERN1(krot,NST,krotst)
! .... test to which kvecvl the new krot belongs
! Try first near current location
        ILow=max(1,j-40)
        Ihigh=min(NKKVL,j+40)
             DO  I=ILow,Ihigh            
               IF(KDELTA(Kvecvl(1,I),NST,krotst)) GOTO 50
             enddo                   
             DO  I=1,NKkvl            
               IF(KDELTA(Kvecvl(1,I),NST,krotst)) GOTO 50
             enddo
             stop 'error: krot not found'
 50          CONTINUE                                                          
!      test if there is a phasef between kvecvl and the new krot
!
             do i1=1,nst
               IF(KDELTA(Kvecvl(1,i),1,krotst(1,i1))) GOTO 55
             enddo
             stop 'error: krotst not found'
 55          continue
!             taup(i1) found
             if(j.lt.50) write(6,51) j,(KVECVL(JX,J),JX=1,3),i,tk,taup(i1)
 51          format(i3,' K-vector:',3i4,' equal',i3,'th K-vector with factor',4f10.5)
              WRITE(18,IFORM2)  (KVECVL(JX,J),JX=1,3),ROKVL(i,1)*tk*taup(i1)    
         enddo       
!        print *,'Done plane waves'
!
! correct scfdn
!
 210  READ(21,700,END=200) MARGN
      IF(MARGN(2:4).EQ.'POS')GO TO 220
      WRITE(22,700)MARGN
      GO TO 210
 220  DO 230 JATO=1,NAT 
         posprint=.false.   
         jatch2=.true.
         do 240 j=1,jmax
         if ((JATO.EQ.jatch(j,1)).or.(JATO.EQ.jatch(j,2))) then
            jatch2=.false.
            REWIND (21)
 260        READ(21,700,END=200) MARGN
            if(margn(2:4).ne.'POS') goto 260
               READ(MARGN,710) JATOM
!               print*, 'jatom,jato,j',jatom,jato,j
            if (jatom.eq.jato) then
               write (22,720) MARGN(1:4),jato,MARGN(8:15),jato,MARGN(20:79)
               if(jatom.eq.jatch(j,1)) jatomdo=jatch(j,2)
               if(jatom.eq.jatch(j,2)) jatomdo=jatch(j,1)
               rewind 21
               goto 262
            endif
            goto 260
 262        READ(21,700,END=200) MARGN
            if(margn(2:4).ne.'POS') goto 262
               READ(MARGN,710) JATOM
!               print*, 'jatom,jatomdo,j',jatom,jatomdo,j
            if (jatom.ne.jatomdo) goto 262

!            if (jatom.eq.jato) then
!               IF(MARGN(2:4).EQ.'POS') write (22,720) MARGN(1:4),jato,MARGN(8:15),jato,MARGN(20:79)
!               posprint=.true. 
!               if (jatch(j,1).ne.jatch(j,2)) goto 260
!            else 
!               if ((jatom.ne.jatch(j,1)).and. &
!                (JATOm.ne.jatch(j,2))) goto 260
!                       if(.not.posprint) goto 260
!            endif
! read posXX, then copy LM lists
!            write (22,720) MARGN(1:4),jato,MARGN(8:15),jato,MARGN(20:79)
 261        READ(21,700)MARGN
            if(Margn(1:1).ne.':') then
              WRITE(22,700)MARGN
              goto 261
            endif
!read CHA,PCS,QTL(91 characters)
            write (22,730) MARGN(1:4),jato,MARGN(8:36),jato,MARGN(40:79)
            READ(21,700)MARGN
            write (22,740) MARGN(1:4),jato,MARGN(8:32),jato,MARGN(37:79)
            READ(21,791)MARGN92
            write (22,751) MARGN(1:4),jato,MARGN92(8:92)
!read EPL, EPH
            READ(21,700)MARGN
            WRITE(22,700)MARGN
            READ(21,700)MARGN
            write (22,750) MARGN(1:4),jato,MARGN(8:79)
            READ(21,700)MARGN
            WRITE(22,700)MARGN
            READ(21,700)MARGN
            write (22,750) MARGN(1:4),jato,MARGN(8:79)
            READ(21,700)MARGN
            IF(MARGN(2:4).EQ.'VZZ')   &
            write (22,760) MARGN(1:4),jato,MARGN(8:26),jato,MARGN(31:79)
!            write(22,2033)
            READ(21,700)MARGN
            write(22,700)MARGN
            IF(MARGN(23:25).EQ.'QXX')   then
              READ(21,700)MARGN
              write(22,700)MARGN
              READ(21,700)MARGN
              write(22,700)MARGN
            endif
            
         endif
 240     continue   
 !        if (jatch2) then
 !250        READ(21,700,END=200) MARGN
 !           WRITE(22,700)MARGN
 !           IF(MARGN(2:4).EQ.'QTL') then
 !           READ(21,700)MARGN
 !           IF(MARGN(2:4).EQ.'VZZ')WRITE(22,700)MARGN 
 !           write(22,2033)
 !           GO TO 230
 !           ENDIF 
 !           GO TO 250
 !        ENDIF 
 230  CONTINUE
 360  READ(21,700,END=200) MARGN
      WRITE(22,700)MARGN
      IF(MARGN(1:7).NE.':CHA  :' ) goto 360
 270  READ(21,700,END=200) MARGN
      WRITE(22,700)MARGN
      IF(MARGN(1:7).NE.':SUM  :' ) goto 270
      forces=.false.
 280  READ(21,700,END=290) MARGN
      IF(MARGN(2:4).EQ.'FVA' ) then
         forces=.true.
         READ(MARGN,710) JATOM
         READ(MARGN,770) (FHELP(jatom,j),j=1,4)
      else if(MARGN(2:4).EQ.'FSU' ) then
         forces=.true.
         READ(MARGN,710) JATOM
         READ(MARGN,770) (FsuHELP(jatom,j),j=1,4)
      else if(MARGN(8:10).NE.'VAL' ) then
         WRITE(22,700)MARGN
      endif
      goto 280
 290  if(.not.forces) goto 200
      DO 300 JATOM=1,NAT    
      write (*,*) (fhelp(jatom,j),j=1,4)
      jato=jatom
        do 310 j=1,jmax
          if(jatch(j,1).eq.jatom ) jato=jatch(j,2)
          if(jatch(j,2).eq.jatom ) jato=jatch(j,1)
 310     continue                                             
!
          FFFX=1.0D0
          FFFY=1.0D0
          FFFZ=1.0D0
          do 320 j=1,jmax
          if ((JATO.EQ.jatch(j,1)).or.(jato.eq.jatch(j,2))) then
            do 325 lj=1,lmchm (j)
            IF((1.EQ.lmch(j,lj,1,1)).and. &
              (0.EQ.lmch(j,lj,2,1))) fffz=-1.
            IF((1.EQ.lmch(j,lj,1,1)).and. &
              (1.EQ.lmch(j,lj,2,1))) fffx=-1.
            IF((-1.EQ.lmch(j,lj,1,1)).and. &
              (1.EQ.lmch(j,lj,2,1))) fffy=-1.
 325      continue
          endif
 320   continue           
      write (*,*) 'fffx,fffy,fffz',fffx,fffy,fffz
!
            WRITE(22,780) 
            write(22,790)jatom,jatom, &
        fhelp(jato,1),fhelp(jato,2)*fffx,fhelp(jato,3)*fffy, &
        fhelp(jato,4)*fffz
 300     continue                      
! now do the same for fsu
      DO 1300 JATOM=1,NAT    
      write (*,*) (fsuhelp(jatom,j),j=1,4)
      jato=jatom
        do 1310 j=1,jmax
          if(jatch(j,1).eq.jatom ) jato=jatch(j,2)
          if(jatch(j,2).eq.jatom ) jato=jatch(j,1)
 1310     continue                                             
!
          FFFX=1.0D0
          FFFY=1.0D0
          FFFZ=1.0D0
          do 1320 j=1,jmax
          if ((JATO.EQ.jatch(j,1)).or.(jato.eq.jatch(j,2))) then
            do 1325 lj=1,lmchm (j)
            IF((1.EQ.lmch(j,lj,1,1)).and. &
              (0.EQ.lmch(j,lj,2,1))) fffz=-1.
            IF((1.EQ.lmch(j,lj,1,1)).and. &
              (1.EQ.lmch(j,lj,2,1))) fffx=-1.
            IF((-1.EQ.lmch(j,lj,1,1)).and. &
              (1.EQ.lmch(j,lj,2,1))) fffy=-1.
 1325      continue
          endif
 1320   continue           
      write (*,*) 'fffx,fffy,fffz',fffx,fffy,fffz
!
            write(22,1790)jatom,jatom, &
        fsuhelp(jato,1),fsuhelp(jato,2)*fffx,fsuhelp(jato,3)*fffy, &
        fsuhelp(jato,4)*fffz
 1300     continue                      
!
  200 CONTINUE

      call dmatcopy(jmax,jatch,symmat,nato)  
      STOP 'clmcopy END'                                                  
!                                                                       
!                                                                       
 77   FORMAT(1X,9(F9.7))                                                
 78   FORMAT(1X,A20,A4,3F10.6,I5)                                       
 700  FORMAT(A79)                                                    
 791  FORMAT(A92)                                                    
 710  FORMAT(4x,i3)                                                    
 720  FORMAT(a4,i3.3,a8,i4,a60)
 730  FORMAT(a4,i3.3,a29,i3,a39) 
 740  FORMAT(a4,i3.3,a25,i4,a43) 
 750  FORMAT(a4,i3.3,a72) 
 751  FORMAT(a4,i3.3,a85) 
 760  FORMAT(a4,i3.3,a19,i4,a46) 
 770  FORMAT(17X,4f15.3)
 780  FORMAT (7x,'VALENCE-FORCE IN mRy/a.u. = |F|',4x,'Fx',13x, &
              'Fy',13x,'Fz')
 790  FORMAT (':FVA',i3.3,': ',i3,'.ATOM',4f15.3)
 1790  FORMAT (':FSU',i3.3,': ',i3,'.ATOM',4f15.3)
 1000 FORMAT(A20)                                                       
 1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4)                                        
 1030 FORMAT(3F10.7)                                                    
 1040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1041 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1055 FORMAT(I4)                                                        
 1056 FORMAT(3(3I2,F5.2,/),I8)                                          
 1980 FORMAT(3X)                                                        
 1990 FORMAT(3X,'ATOM NUMBER=',I3,5X,10A4)                             
 2000 FORMAT(15X,I3//)                                                  
 2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
 2010 FORMAT(15X,I3,5X,I2/)                                             
 2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                         
 2031 FORMAT(/)                                                         
 2032 FORMAT(//)                                                        
 2033 FORMAT(///)                                                       
 2060 FORMAT(/,13X,I6)                                                  
 2061 FORMAT(/,1X,'        ',I10,1X,'NUMBER OF PW')                                     
 3010 format(3a4)
 3011 format(3a4,2f4.1)
      END                                                               
      LOGICAL FUNCTION KDELTA(K,NST,STG)
      include 'param.inc'
!                                                                       
!.... TEST, IF K IS IN STAR OF stg  (GENERATED IN STERN)                  
!                                                                       
!
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER        G,NST,STG(3,nsym)

      DIMENSION      K(3)                                               
!---------------------------------------------------------------------  
      DO 1 I=1,NST                                                      
      DO 2 J=1,3 
      IF(STG(J,I).NE.K(J)) GOTO 1                                       
   2  CONTINUE                                                          
      KDELTA=.TRUE.                                                     
      RETURN                                                            
  1   CONTINUE                                                          
      KDELTA=.FALSE.                                                    
      RETURN                                                            
      END                                                               
      SUBROUTINE STERN (JJA,NST,STG,TAUP,kzz,nkk)                               
!                                                                       
!     STERN GENERATES THE STAR OF REC LATTICE VECTOR KZZ(1,JJA)         
!     THE STAR VECTORS ARE STORED IN STG, THE STAR-SIZE IN  NST         
!     IZ CONTAINS THE SYMMETRY-MATRICES                                 
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16          TAUP,IMAG                                     
      INTEGER          G,STG,INDEX(NSYM)                                
      dimension KZZ(3,Nkk)
      COMMON /SYMETR/  TAU(3,NSYM),IZ(3,3,NSYM),             &
                       INUM(NSYM),IORD                                  
!                                                                       
      DIMENSION        STG(3,NSYM),TAUP(NSYM),G(3)                      
      DATA IMAG        /(0.D0,1.D0)/                                    
!-------------------------------------------------------------------    
!                                                                       
      TPI=2.D0*ACOS(-1.D0)                                              
      G(1)=KZZ(1,JJA)                                                   
      G(2)=KZZ(2,JJA)                                                   
      G(3)=KZZ(3,JJA)                                                   
      NST=0                                                             
!                                                                       
!.....START LOOP OVER ALL SYMMETRY OPERATIONS                           
      DO 1 I=1,IORD                                                     
         TK=0.                                                          
         DO 2 J=1,3                                                     
         TK=TK+TAU(J,I)*G(J)*TPI                                     
           K=0                                                         
           DO 3 L=1,3                                                  
   3       K=IZ(J,L,I)*G(L)+K                                          
           STG(J,I)=K                                                     
 2      continue

         IF(NST.EQ.0) GOTO 7                                            
!                                                                       
!........PROOF, IF THE VECTOR STG(J,I) IS A NEW STARMEMBER OR NOT       
         DO 4 M=1,NST                                                   
            DO 5 J=1,3                                                  
               IF(STG(J,M).NE.STG(J,I)) GOTO 4                          
   5        CONTINUE                                                    
!           STG(J,I) IS NOT A NEW STARMEMBER, IT IS EQUIV TO STG(J,M).  
!           BUT THE TAUPHASE OF STG(J,I) CAN BE NEW.  THEREFORE THE     
!           ALREADY DEFINED PHASE TAUP(M) IS AVERAGED WITH THE PHASE    
!           OF STG(J,M)                                                 
!            TAUP(M)=TAUP(M)+EXP(TK*IMAG)                                
            TAUP(M)=TAUP(M)+CMPLX(cos(TK),sin(TK))
            INDEX(M)=INDEX(M)+1                                         
            GOTO 1                                                      
   4     CONTINUE                                                       
!                                                                       
!........NEW VECTOR FOUND ]]]                                           
   7     NST=NST+1                                                      
         DO 6 J=1,3                                                     
   6     STG(J,NST)=STG(J,I)                                            
!         TAUP(NST)=EXP(TK*IMAG)                                         
         TAUP(NST)=CMPLX(cos(TK),sin(TK))
         INDEX(NST)=1                                                   
   1  CONTINUE                                                          
!                                                                       
      DO 10 I=1,NST           
      TAUP(I)=TAUP(I)/INDEX(I)                                          
  10  continue
      RETURN                                                            
      END                                                               
      subroutine dmatcopy(jmax,jatch,symmat,nato)  
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      dimension symmat(3,3,nato),jatch(nato,2)

      dimension jatom_dmat(nato),jatom2_dmat(nato),l_dmat(nato),sym(3,3)
      complex*16 save_dmat(-3:3,-3:3,nato)
      complex*16,allocatable :: dmatup(:,:),dmatdn(:,:),mat2(:,:) 
!

      jatom=1
 100  read(7,*,end=999,err=999) jatom_dmat(jatom)
      read(7,*) l_dmat(jatom)
      l=l_dmat(jatom)
      allocate (dmatup(-l:l,-l:l), dmatdn(-l:l,-l:l),mat2(-l:l,-l:l))
      do i=-l,l
      read(7,'(2e16.8,2x,2e16.8)') (dmatup(i,j),j=-l,l)
      enddo
!
      sym=0.d0
      do i=1,jmax
      if(jatch(i,1).eq.jatom_dmat(jatom)) then
        sym(1:3,1:3)=symmat(1:3,1:3,i)
        jatom2_dmat(jatom)=jatch(i,2)
      else if(jatch(i,2).eq.jatom_dmat(jatom)) then
        sym(1:3,1:3)=transpose(symmat(1:3,1:3,i))
        jatom2_dmat(jatom)=jatch(i,1)
      endif
      enddo
!
      call determinant(sym,det)
      if (det.lt.0d0) then
        sym(1:3,1:3)=-sym(1:3,1:3)
        write(6,'(a)') 'sym is  reflection, so multiply by inversion:'
        write(6,'(3f15.8)') transpose(sym)
        write(6,*)      
      endif
      call euler(1,sym,a2,b2,c2)
      call find_rot_mat(l,a2,b2,c2,mat2,l)
      if (det.lt.0) call apply_inversion_Ylm(l,l,mat2)
                 dmatdn=conjg(dmatup)
                 dmatdn=matmul(mat2,dmatdn)
                 mat2=transpose(conjg(mat2))
                 dmatdn=matmul(dmatdn,mat2)
                 dmatdn=conjg(dmatdn)
                 save_dmat(-l:l,-l:l,jatom_dmat(jatom))=dmatdn
                 deallocate (dmatup, dmatdn, mat2)

     jatom=jatom+1
     goto 100
!
 999    if(jatom.eq.1) then
          write(*,*) 'case.dmatup not present'
          return
        else
          write(*,*) 'case.dmatdn created'
          do jatom1=1,jatom-1
          write(8,'(i5," atom density matrix" )') jatom_dmat(jatom1)
          write(8,'(i5,"  0.000000  0.000000  0.000000 L, Lx,Ly,Lz in global orthogonal system")') l_dmat(jatom1)
           do i=-l_dmat(jatom1),l_dmat(jatom1)
           write(8,'(2e16.8,2x,2e16.8)') (save_dmat(i,j,jatom2_dmat(jatom1)),j=-l_dmat(jatom1),l_dmat(jatom1))
           enddo
          enddo
        endif
        end


 subroutine determinant(sym,det)
  
  real*8   sym(3,3),det

  det=sym(1,1)*(sym(2,2)*sym(3,3)-sym(3,2)*sym(2,3))-&
      sym(1,2)*(sym(2,1)*sym(3,3)-sym(3,1)*sym(2,3))+&
      sym(1,3)*(sym(2,1)*sym(3,2)-sym(2,2)*sym(3,1))

 end subroutine determinant
 subroutine vecprod(a,b,c)
   implicit real*8 (a-h)
   dimension a(3),b(3),c(3)
   
   c(1)=a(2)*b(3)-a(3)*b(2)
   c(2)=a(3)*b(1)-a(1)*b(3)
   c(3)=a(1)*b(2)-a(2)*b(1)
 end subroutine vecprod


 
 real*8 function dot(a,b)
   implicit real*8 (a-h)
   dimension a(3),b(3)
   dot=0
   do i=1,3
      dot=dot+a(i)*b(i)
   end do
 end function dot
subroutine euler(mod,rot_new,a,b,c)
  implicit real*8 (a-h,o-z)
  
  integer mod
  real tmp

  dimension rot_old(3,3),rot_new(3,3),z(3),zz(3)
  dimension y(3),yy(3),yyy(3),pom(3),x(3),xx(3)
  
  zero=1.0d-20
  pi=2.0d0*acos(0.0d0)
  
  do i=1,3
     do j=1,3
        rot_old(i,j)=0
        if (i.eq.j) rot_old(i,i)=1
     end do
  end do
3 format(20x,3(f10.8,2x))
  
  do j=1,3
     y(j)=rot_old(j,2)
     yyy(j)=rot_new(j,2)
     z(j)=rot_old(j,3)
     zz(j)=rot_new(j,3)
  end do
  
  call vecprod(z,zz,yy)
  y_norm=dsqrt(dot(yy,yy))
  
  if (y_norm.lt.1d-20) then
     
     if (abs(dot(y,yyy)).gt.1d0) then
        aa=dot(y,yyy)/abs(dot(y,yyy))
        a=acos(aa)
     else
        a=acos(dot(y,yyy))
     end if
     
     if (dot(z,zz).gt.zero) then
        c=zero
        b=zero
        if (yyy(1).gt.zero) a=2.0d0*pi-a
     else
        c=a
        a=zero
        b=pi
        if (yyy(1).lt.zero) c=2.0d0*pi-c
     end if
     
  else
     
     do j=1,3
        yy(j)=yy(j)/y_norm
     end do
     
     aa=dot(y,yy)
     bb=dot(z,zz)
     cc=dot(yy,yyy)
     if (abs(aa).gt.1d0) aa=aa/abs(aa)
     if (abs(bb).gt.1d0) bb=bb/abs(bb)
     if (abs(cc).gt.1d0) cc=cc/abs(cc)
     b=acos(bb)
     a=acos(aa)
     c=acos(cc)
     if (yy(1).gt.zero) a=2.0d0*pi-a
     call vecprod(yy,yyy,pom)
     if (dot(pom,zz).lt.-zero) c=2.0d0*pi-c
  end if

 if (mod.eq.1) then
!   space fixed axes (exchange a,c)
    tmp=a
    a=c
    c=tmp
 endif 


end subroutine euler
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

 subroutine find_rot_mat(l,a,b,c,mat,lmax)

  integer l,lmax
  real*8  a,b,c
  complex*16 mat(-lmax:lmax,-lmax:lmax)
  integer n,m
  complex*16 imag
  real*8  dlmat 

  imag=(0.0d0,1.0d0)
  mat=(0.0d0,0.0d0)
  do n=-l,l
     do m=-l,l
        call find_dlmat(l,n,m,b,dlmat) 
        tmp1=dble(n)*c
        tmp2=dble(m)*a
        mat(n,m)=cmplx(dlmat*cos(tmp1),-dlmat*sin(tmp2))
!        mat(n,m)=exp(-imag*n*c)*dlmat*exp(-imag*m*a)
     enddo
  enddo

 end subroutine find_rot_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine find_dlmat(l,n,m,b,dlmat)

  integer l,n,m
  real*8  dlmat,b,fact,b2
  integer k

  b2=b/2.0d0
  dlmat=0.0d0
  do k=0,2*l
 
     if (((l-n-k).ge.0).and.((l+m-k).ge.0).and.((k-m+n).ge.0).and.&
         ((l+n).ge.0).and.((l+m).ge.0).and.((l-n).ge.0).and.((l-m).ge.0)) then

         dlmat=dlmat+((-1.0d0)**(k-m+n))*&
             (cos(b2)**(2*l+m-n-2*k))*(sin(b2)**(2*k+n-m))*&
             sqrt(fact(l+n)*fact(l+m)*fact(l-n)*fact(l-m))/&
             (fact(l-n-k)*fact(l+m-k)*fact(k)*fact(k-m+n))     

     endif

  enddo

 end subroutine find_dlmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real*8 function fact(n)
   
   integer n,j

   if (n.eq.0) then
      fact=1
   else
      fact=1
      do j=1,n
         fact=fact*j
      end do
   end if

 end Function fact
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine apply_inversion_Ylm(l,lmax,mat)

   integer l
   complex*16 mat(-lmax:lmax,-lmax:lmax)

   mat=((-1.0d0)**l)*mat
   
 end subroutine apply_inversion_Ylm
