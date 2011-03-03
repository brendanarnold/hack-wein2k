      PROGRAM afminput

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'   

!     *****************************************************************
!     *****************************************************************
!
!     PARAMETERS
!
!     DON'T JUST CHANGE VALUES HERE, BUT ALSO WITHIN THE 
!     SUBROUTINES!!!
!     (e.g. use a query-replace command for the whole line)
!
!     *****************************************************************
!     *****************************************************************

!     NATO stands for maximum Number of ATOms  (per cell, ).


!     NCOM stands for maximum Number of LM (LM for a combination
!     of the two quantum numbers, default=48).

!     LOG: Set it to value 1, if you want a logfile to be created,
!     0 if not (default=1).

      PARAMETER(LOG=1)

!     *****************************************************************
!     *****************************************************************
 
      
      CHARACTER*77   FNAME
      CHARACTER*11   STATUS,FORM  
      CHARACTER*1    TMP1
      CHARACTER*2    TMP2,TMP21,TMP22
      CHARACTER*3    TMP3
      CHARACTER*4    TMP4
      CHARACTER*6    TMP6,TMP61
      CHARACTER*12   STR1,STR2
      character*4    lattic,lattic1

      INTEGER,allocatable ::        NDATAP(:), LMSEQ(:,:,:), &
                     NDATAL(:),NLM(:)
      INTEGER,allocatable ::        LMOUT(:,:),LMOUT1(:,:),CHATI(:,:)
      real*8,allocatable ::  DATA(:,:,:),clmfac(:,:),symout(:,:,:)
      real*8,allocatable :: pos1(:),pos2(:),pos3(:),locrot(:,:,:)
      integer,allocatable :: mult(:),iatomlist(:),iatomlist1(:),iatomlist2(:)

      LOGICAL        EQUAL


      dimension sym(3,3),rrot(3),rrott(3),rot(3,3),rotinv(3,3),rottrans(3,3)
      dimension br1(3,3),br1inv(3,3)
      COMPLEX*16 tmat(2*lmax+1,2*lmax+1),tmat1(-lmax:lmax,-lmax:lmax),tmat2(-lmax:lmax,-lmax:lmax)
      integer lsym0(3,3)
      REAL*8 TAU(3,48),tau0(3)
      COMMON /SYMM/  LSYM(3,3,48),IORD,INUM(48)                     
      REAL*8 TAU1(3,48)
      COMMON /SYMM1/  LSYM1(3,3,48),IORD1,INUM1(48)                     



!     *****************************************************************
!     *****************************************************************
!
!     BEGIN OF MAIN PROGRAM
!
!     *****************************************************************
!     *****************************************************************
      do i=1,3
      do j=1,3
      lsym0(i,j)=0
      enddo
      enddo

!     Opening of file stated in the argument.

      CALL GETARG(2,FNAME)
      IF (FNAME .EQ. '      ') CALL GETARG(1,FNAME)
      OPEN(1,FILE=FNAME,STATUS='OLD',ERR=20)


!     Opening of files listed in the file opened above.


 10   READ(1,*,END=40) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=30)
      IF (LOG .EQ. 1) WRITE(14,'("FILE OPENED SUCCESSFULLY: ",A20)') &
                                  FNAME 
      GOTO 10
 20   WRITE(*,*) ' ERROR IN OPENING ', FNAME, ' !!!!'
      STOP 'afminput.DEF'
 30   WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 40   CONTINUE

!     *****************************************************************
!
!     READING INFORMATION FROM THE *.STRUCT FILE 
!     (number of atoms and lines per data-block)
!     
!     *****************************************************************


      READ(20,9001) lattic,NATOMS
      read(20,'(6f10.5)') alat,blat,clat,alpha,beta,gamma
      READ(21,9001,iostat=isuper) lattic1,NATOMS1
      if(isuper.eq.0) print *, 'case.struct_supergroup present'
      if(isuper.ne.0) then
      print *, 'case.struct_supergroup NOT present!!!'
      print *, 'It is strongly recommended that you copy the (nonmagnetic) supergroup '
      print *, 'struct file to case.struct_supergroup (unless they are KLASSENGLEICH)'
      print *, 'Otherwise:'
      endif
      nato=natoms  
      allocate (        NDATAP(NATO), LMSEQ(NATO,NCOM,2), &
                  NDATAL(NATO),NLM(NATO))
      allocate (        LMOUT(NATO,NCOM*2),CHATI(NATO,3),symout(3,3,nato))
      allocate (        LMOUT1(NATO,NCOM*2),Clmfac(NATO,ncom*2))
      allocate (  DATA(2,NATO,NCOM),MULT(NATO))
      allocate (pos1(48*nato),pos2(48*nato),pos3(48*nato),iatomlist(48*nato))
      allocate (iatomlist1(nato),iatomlist2(nato),locrot(3,3,nato) )
!      IF (NATOMS .GT. NATO) STOP 'NATO too small'
      IF (LOG .EQ. 1) WRITE(14,'(/,"NUMBER OF ATOMS: ",I3,/)') NATOMS
      IF (LOG .EQ. 0) WRITE(14,*) 'SET THE LOG PARAMETER TO 1, IF', &
                      ' YOU WANT A LOG-FILE TO BE CREATED!!!'
!
      INDEX=0                                                           
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
!                                                                       
      nat=natoms
      DO 50 JATOM = 1,NAToms                                                
         INDEX=INDEX+1                                                  
         READ(20,1040) IATNR,POS1(index),POS2(index),             &
                       POS3(index),MULT(JATOM),ISPLIT            
            MULTND=MULT(JATOM)-1                                        
            DO 11 M=1,MULTND                                            
            INDEX=INDEX+1                                               
 1040 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1041 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
            READ(20,1041) IATNR,POS1(index),POS2(index),          &
                          POS3(index)                                   
 11         CONTINUE                                                    
         READ(20,1050) ANAME,NDATAP(JATOM),RNOT,RMTT,         &
                       ZZ                                        
         NDATAL(jatom)=NDATAP(jatom)/4
         IF (NDATAL(jatom)*4 .NE. NDATAP(jatom)) NDATAL(jatom)=NDATAL(jatom)+1
         READ(20,'(20x,3f10.5)') (locrot(1,j1,jatom), J1=1,3 ) 
         READ(20,'(20x,3f10.5)') (locrot(2,j1,jatom), J1=1,3 ) 
         READ(20,'(20x,3f10.5)') (locrot(3,j1,jatom), J1=1,3 ) 
 50   CONTINUE                                                          
!     READ SYMMETRY-OPERATIONS FROM TAPE20=POTE                         
      READ(20,1055) IORD                                                
      DO 12 J=1,IORD                                                    
 12   READ(20,1056) ( (LSYM(J1,J2,J),J2=1,3),TAU(J1,J), J1=1,3 ),INUM(J)
 1055 FORMAT(I4)                                                        
 1056 FORMAT(3(3I2,F10.7,/),I8)                                          
!
!     get the LM values for all atoms
      CALL GETDATA(1,16,NATOMS,DATA,LMSEQ,NDATAL,nlm)
!
! check if case.struct_supergroup exists
      if(isuper.eq.0) then
      DO 51 i=1,NATOMS1
 61      READ(21,9004,END=51) TMP3,TMP4
         IF (TMP3 .EQ. 'NPT') GOTO 71 
         GOTO 61
 71      continue
 51   CONTINUE
      read(21,*)
      read(21,*)
      read(21,*)
!     READ SYMMETRY-OPERATIONS FROM TAPE20=POTE                         
      READ(21,1055) IORD1                                                
      DO  J=1,IORD1                                                    
      READ(21,1056) ( (LSYM1(J1,J2,J),J2=1,3),TAU1(J1,J), J1=1,3 ),INUM1(J)
      enddo
      if(iord.eq.iord1) then
            print *, 'The super and subgroups are KLASSENGLEICH'
            print *, 'You must specify a translation vector which transforms the'
            print *, 'spin-up into the spin-dn atom (e.g.: 0.5,0.5,0.5 for AFM bcc Cr)!'
            read(*,*) tau0
            lsym0(1,1)=1
            lsym0(2,2)=1
            lsym0(3,3)=1
      else
            print *, 'The super and subgroups are TRANSLATIONENGLEICH'
            do 56 j=1,iord1
              do 55 i=1,iord
                do j1=1,3
                do j2=1,3
                if(lsym1(j1,j2,j).ne.lsym(j1,j2,i)) goto 55
                enddo
                enddo
                goto 56
 55           continue
                print *, 'Found a symmetry operation:'
                do j1=1,3
                tau0(j1)=tau1(j1,j)
                do j2=1,3
                lsym0(j1,j2)=lsym1(j1,j2,j)
                enddo
                write(6,'(3i4,f10.5)') (lsym0(j1,j2),j2=1,3),tau0(j1)
                enddo
                goto 57
 56         continue  
 57         continue
      endif
! no case.struct_supergroup present
      else
            print *, 'You must specify a symmetry operation (rotation + translation vector) '
            print *, 'which transforms the spin-up into the spin-dn atom (e.g. for AFM bcc Cr:)!'
            print *, ' 1 0 0  0.5'
            print *, ' 0 1 0  0.5'
            print *, ' 0 0 1  0.5'
            do j1=1,3
            read(5,*) (lsym0(j1,j2),j2=1,3),tau0(j1)
            enddo
      endif
!
!      Ensure that lsymm1 has unitary as first -- for simplicity
       do j=1,iord1
                if( (abs(lsym(1,1,j)-1).gt.1d-4) &
                .or.(abs(lsym(2,2,j)-1).gt.1d-4) &
                .or.(abs(lsym(3,3,j)-1).gt.1d-4)) goto 1001

                if( (abs(lsym(1,2,j)).gt.1d-4).or.(abs(lsym(1,3,j)).gt.1d-4) &
                .or.(abs(lsym(2,1,j)).gt.1d-4).or.(abs(lsym(2,3,j)).gt.1d-4) &
                .or.(abs(lsym(3,1,j)).gt.1d-4).or.(abs(lsym(3,2,j)).gt.1d-4)) goto 1001
                iunit=j
                goto 1002
1001            continue
       enddo
!      Re-order
1002   do j=iunit,2,-1
                do j1=1,3
                do j2=1,3
                        lsym(j1,j2,j)=lsym(j1,j2,j-1)
                enddo
                enddo
       enddo
        do j1=1,3
        do j2=1,3
                lsym(j1,j2,1)=0
        enddo
        lsym(j1,j1,1)=1
        enddo


!      OPEN(6,FILE='afminput.outputafminput',STATUS='unknown',FORM='formatted')
!find equivalent atoms due to afm sym oper.
!
      nchat=0
      index=0
      do jatom=1,nat
      write(14,*) ' '
      do m=1,mult(jatom)
      index=index+1
      if(m.ne.1) goto 15
         DO  J1=1,3                                                     
         rrott(J1)=lsym0(j1,1)*pos1(index)+lsym0(j1,2)*pos2(index)+lsym0(j1,3)*pos3(index)+TAU0(J1)
         if(rrott(j1).gt.0.99999999999999d0) rrott(j1)=rrott(j1)-1.d0
         if(rrott(j1).lt.-0.0000000000001d0) rrott(j1)=rrott(j1)+1.d0
         ENDDO
         index1=0
         do jatom1=1,nat
         do m1=1,mult(jatom1)
         index1=index1+1
!check rotated position
         rrot=rrott
         jrot0=0
 145     continue
         write(14,'(3(2f10.5,5x))') rrot(1),pos1(index1),rrot(2),pos2(index1),rrot(3),pos3(index1)
         if(abs(rrot(1)-pos1(index1)).lt.1.d-5.and.abs(rrot(2)-pos2(index1)).lt.1.d-5.and.abs(rrot(3)-pos3(index1)).lt.1.d-5) then       
!           iatomlist(index)=index1
           iatomlist1(jatom)=jatom1
           write(14,16) index,pos1(index),pos2(index),pos3(index),index1,rrot(1),rrot(2),rrot(3) 
 16        format('atom',i4,3f8.4,' transfered to atom',i4,3f8.4)
           if(jatom.le.jatom1) then
           nchat=nchat+1
           CHATI(nchat,1)=jatom
           CHATI(nchat,2)=jatom1
           endif
! find symop which transforms rotated atom into first position
              do j=1,iord
                DO  J1=1,3                                                     
                rrot(J1)=lsym(j1,1,j)*pos1(index1)+lsym(j1,2,j)*pos2(index1)+lsym(j1,3,j)*pos3(index1)+TAU(J1,j)
                if(rrot(j1).gt.0.99999999999999d0) rrot(j1)=rrot(j1)-1.d0
                if(rrot(j1).lt.-0.0000000000001d0) rrot(j1)=rrot(j1)+1.d0
                ENDDO
                irot0=0
 144            continue
         write(14,'(3(2f10.5,2x),3i3)') rrot(1),pos1(index1-m1+1),rrot(2),pos2(index1-m1+1),rrot(3),pos3(index1-m1+1),index1-m1+1,index1,m1
                if(abs(rrot(1)-pos1(index1-m1+1)).lt.1.d-5.and. &
                   abs(rrot(2)-pos2(index1-m1+1)).lt.1.d-5.and. &
                   abs(rrot(3)-pos3(index1-m1+1)).lt.1.d-5) then   
              iatomlist2(jatom)=j
              write(14,*) 'rotation to first atom by operation',iatomlist2(jatom)                
                goto 15
                else
                 call shiftpos(rrot,lattic,irot0)
                 if(irot0.ne.0) goto 144
                endif
              enddo
              write(14,*) 'rotation to first atom not found'
              stop 'rotation to first atom not found'
           goto 15
         else
          call shiftpos(rrot,lattic,jrot0)
          if(jrot0.ne.0) goto 145
         endif
         enddo
         enddo
         write(14,*) 'equivalent position not found',pos1(index),pos2(index),pos3(index),rrot(1),rrot(2),rrot(3)
         stop 'rrot not found'
 15   continue
      enddo
      enddo
!
      call latgen(lattic,alat,blat,clat,alpha,beta,gamma,br1)
      CALL INVERSSYMDEF(BR1,BR1inv)
     
! find rotations for all atomic positions
!
      do ichat=1,nchat
      jatom=CHATI(ichat,1)
!      do jatom=1,nat
        write(14,*) ' '
        write(14,*) 'symmetries for atom',jatom
        pi=acos(-1.d0)
        sym=0
        do i=1,3
        do j=1,3
         do k=1,3
         sym(i,j)=sym(i,j)+lsym0(i,k)*lsym(k,j,iatomlist2(jatom))
!         sym(i,j)=sym(i,j)+lsym(i,k,iatomlist2(jatom))*lsym0(k,j)
        enddo
        enddo
        write(14,'(3f10.4,10x,3i5)') (sym(i,j),j=1,3),(lsym(i,j,iatomlist2(jatom)),j=1,3)
        enddo
!! transfer to local rotation matrix (moved after orthogonal transf.)
!        rot(1:3,1:3)=transpose(locrot(1:3,1:3,jatom))
!!        rot(1:3,1:3)=(locrot(1:3,1:3,jatom))
!         call INVERSSYMDEF(rot,rotinv)
!        write(14,'(3f10.4,10x,3f10.4)') ((rot(i,j),j=1,3),(rotinv(i,j),j=1,3),i=1,3)
!         call transform(rottrans,rot,sym,rotinv)
!         sym=rottrans
!        do i=1,3
!        write(14,'(3f10.4,10x,3i5)') (sym(i,j),j=1,3),(lsym(i,j,iatomlist2(jatom)),j=1,3)
!        enddo
!transfer to orthogonal system
          rot=sym
          sym=0.d0
          do j=1,3
          do k=1,3
            do l=1,3
            do m=1,3
          sym(j,k)= sym(j,k) +  &
               BR1(j,l) * rot(l,m) * BR1inv(m,k)
            end do
            end do
          end do
          end do
! transfer to local rotation matrix
        rot(1:3,1:3)=transpose(locrot(1:3,1:3,jatom))
!        rot(1:3,1:3)=(locrot(1:3,1:3,jatom))
         call INVERSSYMDEF(rot,rotinv)
        write(14,'(3f10.4,10x,3f10.4)') ((rot(i,j),j=1,3),(rotinv(i,j),j=1,3),i=1,3)
         call transform(rottrans,rot,sym,rotinv)
         sym=rottrans
        do i=1,3
        write(14,'(3f10.4,10x,3i5)') (sym(i,j),j=1,3),(lsym(i,j,iatomlist2(jatom)),j=1,3)
        enddo
!
        if(br1(1,1).ne.1.d0.and.br1(2,2).ne.1.d0.and.br1(3,3).ne.1.d0) then
           write(14,*) 'symmetry operation in orthogonal coordinates'
           do i=1,3
            write(14,'(3f10.4)') (sym(i,j),j=1,3)
           enddo
        endif 
        symout(1:3,1:3,ichat)=sym(1:3,1:3)
        dd =sym(1,1)*(sym(2,2)*sym(3,3)- &
                         sym(3,2)*sym(2,3))- &
             sym(1,2)*(sym(2,1)*sym(3,3)- &
                         sym(3,1)*sym(2,3))+ &
             sym(1,3)*(sym(2,1)*sym(3,2)- &
                         sym(2,2)*sym(3,1))
!_____________________________________________________________
	if (dd.lt.0d0) then
        do i=1,3
        do j=1,3
        sym(i,j)=-sym(i,j)
        end do
        end do
        end if
 	call euler(0,sym,a11,b11,c11)
   write(14,3) a11*180/pi,b11*180/pi,c11*180/pi,dd
 3  format(2x,'euler angles: a,b,c ',3f6.1,' det:',f10.5)
        l0=1
	do  l=1,lmax
            call find_rot_mat(l,a11,b11,c11,tmat1,lmax)
        if (dd.lt.0) call apply_inversion_Ylm(l,lmax,tmat1)
           write(14,*) 'dd-matrix for l=',l,dd
           do m=-l,l
               do n=-l,l
               tmat(l+m+1,l+n+1)=tmat1(m,n)
               enddo 
            write(14,5)(tmat(l+m+1,l+n+1),n=-l,l)
           end do
 5         format(9(2f7.3,2x))
           DO 17 k=2,NLM(jatom)
             if(abs(LMSEQ(jatom,k,1)).eq.l) then
             n=LMSEQ(jatom,k,2)
              if(LMSEQ(jatom,k,2).eq.0) then
!m.eq.0
                jrtmat_pp=nint(dble(tmat(l+n+1,l+n+1)))
                jitmat_pp=nint(aimag(tmat(l+n+1,l+n+1)))
                write(14,*) 'test lm',LMSEQ(jatom,k,1),LMSEQ(jatom,k,2)
                write(14,'(4(2i3,2x))') jrtmat_pp,jitmat_pp
                if(jrtmat_pp.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                endif
              goto 17
              endif
             jrtmat_pp=nint(dble(tmat(l+n+1,l+n+1)))
             jitmat_pp=nint(aimag(tmat(l+n+1,l+n+1)))
             jrtmat_mm=nint(dble(tmat(l-n+1,l-n+1)))
             jitmat_mm=nint(aimag(tmat(l-n+1,l-n+1)))
             jrtmat_pm=nint(dble(tmat(l+n+1,l-n+1)))
             jitmat_pm=nint(aimag(tmat(l+n+1,l-n+1)))
             jrtmat_mp=nint(dble(tmat(l-n+1,l+n+1)))
             jitmat_mp=nint(aimag(tmat(l-n+1,l+n+1)))
             write(14,*) 'test lm',LMSEQ(jatom,k,1),LMSEQ(jatom,k,2)
             write(14,'(4(2i3,2x))') jrtmat_pp,jitmat_pp,jrtmat_mm,jitmat_mm,jrtmat_pm,jitmat_pm,jrtmat_mp,jitmat_mp
!nondiagonal dd matrix
!case +l, even m
              if(LMSEQ(jatom,k,1).gt.0.and.(n/2)*2.eq.n) then
                if(jrtmat_pp.eq.1.and.jrtmat_mm.eq.1) then
                   goto 17
                else if(jrtmat_pp.eq.-1.and.jrtmat_mm.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.1.and.jitmat_mm.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm listm1'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.-1.and.jitmat_mm.eq.1) then
                   clmfac(iCHAT,l0)=1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm listm2'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                else if(jrtmat_pm.eq.1.and.jrtmat_mp.eq.1) then
                   goto 17
                else if(jrtmat_pm.eq.-1.and.jrtmat_mp.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.1.and.jitmat_mp.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm list1'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.-1.and.jitmat_mp.eq.1) then
                   clmfac(iCHAT,l0)=1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm list2'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                endif
!case -l, even m
              else if(LMSEQ(jatom,k,1).lt.0.and.(n/2)*2.eq.n) then
                if(jrtmat_pp.eq.-1.and.jrtmat_mm.eq.-1) then
!back for fayalit!interchanged for k2v3o8 (-4 4)  
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jrtmat_pp.eq.1.and.jrtmat_mm.eq.1) then
!for feooh!back for fayalit!cuo1,2   muss +1 sein (identitaet!)
!                   clmfac(iCHAT,l0)=-1.d0
 !                  LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
  !                 LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
   !                LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
    !               LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
     !              l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.1.and.jitmat_mm.eq.-1) then
                   clmfac(iCHAT,l0)=1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm listm3'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.-1.and.jitmat_mm.eq.1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm listm4'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                else if(jrtmat_pm.eq.-1.and.jrtmat_mp.eq.-1) then
!back ! cuo 3/4 1 1 and -1 -1 rule interchanged
                   goto 17
                else if(jrtmat_pm.eq.1.and.jrtmat_mp.eq.1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.1.and.jitmat_mp.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm list3'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.-1.and.jitmat_mp.eq.1) then
                   clmfac(iCHAT,l0)=1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm list4'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                endif
!case l, odd m
              else if(LMSEQ(jatom,k,1).gt.0.and.(n/2)*2.ne.n) then
                if(jrtmat_pp.eq.1.and.jrtmat_mm.eq.1) then
!for Fayalit last atom
!                   clmfac(iCHAT,l0)=-1.d0  !op8
!                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
!                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
!                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
!                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
!                   l0=l0+2
                   goto 17
                else if(jrtmat_pp.eq.-1.and.jrtmat_mm.eq.-1) then
!ja for feooh 3.atom
                   clmfac(iCHAT,l0)=-1.d0  !op8
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.1.and.jitmat_mm.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0    
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm listm5'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.-1.and.jitmat_mm.eq.1) then
                   clmfac(iCHAT,l0)=1.d0   
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm listm6'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                else if(jrtmat_pm.eq.-1.and.jrtmat_mp.eq.-1) then
                   goto 17
                else if(jrtmat_pm.eq.1.and.jrtmat_mp.eq.1) then
!back!cuo 1,2 (i) interchanged
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.1.and.jitmat_mp.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm list5'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.-1.and.jitmat_mp.eq.1) then
                   clmfac(iCHAT,l0)=1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k+1,1)) stop 'lm list6'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k+1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k+1,2)
                   l0=l0+2
                   goto 17
                endif
! case -l, odd m
              else if(LMSEQ(jatom,k,1).lt.0.and.(n/2)*2.ne.n) then
                if(jrtmat_pp.eq.-1.and.jrtmat_mm.eq.-1) then
!fayalit                   clmfac(iCHAT,l0)=1.d0     !op8
                   clmfac(iCHAT,l0)=-1.d0    !op7
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jrtmat_pp.eq.1.and.jrtmat_mm.eq.1) then
!cuo1,2   muss +1 sein (identitaet!)
!                   clmfac(iCHAT,l0)=-1.d0
 !                  LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
  !                 LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
   !                LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
    !               LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
     !              l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.1.and.jitmat_mm.eq.-1) then
!                   clmfac(iCHAT,l0)=-1.d0   !op8
                   clmfac(iCHAT,l0)=1.d0   !op7
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm listm7'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pp.eq.-1.and.jitmat_mm.eq.1) then
!                   clmfac(iCHAT,l0)=1.d0   !op8
                   clmfac(iCHAT,l0)=-1.d0   !op7
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm listm8'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                else if(jrtmat_pm.eq.1.and.jrtmat_mp.eq.1) then
!back ! cuo 1/2 1 1 and -1 -1 rule interchanged
                   goto 17
                else if(jrtmat_pm.eq.-1.and.jrtmat_mp.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.1.and.jitmat_mp.eq.-1) then
                   clmfac(iCHAT,l0)=-1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm list7'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                else if(jitmat_pm.eq.-1.and.jitmat_mp.eq.1) then
                   clmfac(iCHAT,l0)=1.d0
                   LMOUT(iCHAT,l0)=LMSEQ(jatom,k,1)
                   LMOUT(iCHAT,l0+1)=LMSEQ(jatom,k,2)
                   if(LMSEQ(jatom,k,1).ne.-LMSEQ(jatom,k-1,1)) stop 'lm list8'
                   LMOUT1(iCHAT,l0)=LMSEQ(jatom,k-1,1)
                   LMOUT1(iCHAT,l0+1)=LMSEQ(jatom,k-1,2)
                   l0=l0+2
                   goto 17
                endif
              endif
             stop 'lm list0'
             endif
 17        continue

        enddo
            CHATI(iCHAT,3)=(L0-1)/2

      enddo
!
!     Output to unit 15 (*.inclmcopy)
      WRITE(15,'(i4,46x,''NUMBER of ATOMS to CHANGE'')') NCHAT
      DO 110 l=1,NCHAT
         WRITE(15,'(2i4,42x,''INTERCHANGE these ATOMS'')') CHATI(l,1), &
                    CHATI(l,2)
         write(15,'(3f15.11,5x,"SYMMETRY OPERATION",/,3f15.11,/3f15.11)') ((symout(i,j,l),j=1,3),i=1,3)
         WRITE(15,'(i4,46x,''NUMBER of LM to CHANGE SIGN'')') CHATI(l,3)
          do m=1,CHATI(l,3)*2,2
         WRITE(15,'((i3,i2,2x,i3,i2,f6.2))') LMOUT(l,m),LMOUT(l,m+1),LMOUT1(l,m),LMOUT1(l,m+1),clmfac(l,m)
          enddo
 110  CONTINUE


      do j1=1,3
      write(15,'(3i4,f10.5)') (lsym0(j1,j2),j2=1,3),tau0(j1)
      enddo
      write(15,*) 'The symmetry operation above is one of the operations of the NM-supergroup'
      write(15,*) 'missing in the AFM-subgroup (transfers spin-up into spin-dn atom)'
      stop

!     FORMAT-Expressions

 9000 FORMAT(///)
 9001 FORMAT(/,a4,24X,I2,/)
 9002 FORMAT(15X,I5)
 9003 FORMAT(////)
 9004 FORMAT(11X,A3,2X,A4)
 9005 FORMAT(/)
 9006 FORMAT(I6)
 9007 FORMAT(I2,1X,I2,1X)
 9008 FORMAT(2X,A6,5X,A6)
 9009 FORMAT(3X,3I5,F19.12)
 9010 FORMAT(A4,A4,A4,2f4.1)
 9011 FORMAT(A19,I3,A1,I3,3X,A2,I3,3X,A2,i3,2(2X,E16.9))
 9012 FORMAT(A13,1X,3A4,1X,3I5,2(2X,E16.9))
      END


!     *****************************************************************
!     *****************************************************************
!
!     SUBROUTINES 
!
!     *****************************************************************
!     *****************************************************************


!     *****************************************************************
!
!     SUBROUTINE GEDTATA
!
!     *****************************************************************


!     read-in of data from LM-data-blocks  from the *.clmup/*.clmdn 
!     files

      SUBROUTINE GETDATA(number,iunit,NATOMS,data,lmseq,ndatal,nlm)

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'   

!     Parameters have to be edited again.

!      PARAMETER(NATO=100)
!      PARAMETER(NCOM=48) 
      PARAMETER(LOG=1)

      CHARACTER*1    TMP1
      CHARACTER*2    TMP2,TMP21,TMP22
      CHARACTER*3    TMP3
      CHARACTER*4    TMP4
     
      dimension   DATA(2,NATOms,NCOM), LMSEQ(NATOms,NCOM,2) &
                     ,NDATAL(NATOms),NLM(NATOms)



!     Outer loop (i) runs through the number of Atoms

      DO 500 i=1,NATOMS

!     Scanning for LM-string. This should be less error-sensitive than 
!     counting lines.

 510      READ(iunit,9902,END=500) TMP2,TMP3
         IF (TMP2 .EQ. 'LM') GOTO 520
         GOTO 510
 520     READ(TMP3,9903) NLM(i)
         IF (LOG .EQ. 1)  &
            WRITE(14,'(/,"NUMBER OF LM IN ATOM ",I3,A9,I2,A3,I3)') &
                  i,' IN UNIT ',iunit,': ',NLM(i)
             if(nlm(i).gt.ncom) stop 'NCOM too small'

!     With the knowledge of the number of LM, the inner loop for each
!     atom can be started

!     First the scanning for the line holding information about 
!     the L and M values

         DO 530 j=1,NLM(i)
 540        READ(iunit,9904,END=530) TMP1,TMP3,TMP21
            IF (TMP1 .EQ. 'L') GOTO 550
            GOTO 540
 550        READ(TMP3,9903) LMSEQ(i,j,1)
            READ(TMP21,9903) LMSEQ(i,j,2)


!     Now read in of the first value in the last line of the data-
!     block.

            DO 560 k=1,NDATAL(i)
               READ(iunit,*,END=560)
 560        CONTINUE
            READ(iunit,9905) VALUE
            DATA(number,i,j)=VALUE
            IF (LOG .EQ. 1) WRITE(14,9906)  &
             'L=',LMSEQ(i,j,1),'M=', &
                 LMSEQ(i,j,2),DATA(number,i,j)            

 530     CONTINUE
 500  CONTINUE



!     FORMAT-EXPRESSIONS for GETDATA


 9900 FORMAT(///)
 9901 FORMAT(/)
 9902 FORMAT(13X,A2,A3)
 9903 FORMAT(I3)
 9904 FORMAT(14X,A1,A3,5X,A2)
 9905 FORMAT(3X,E19.12)
 9906 FORMAT(A2,I3,5X,A2,I3,5X,E19.12)
      END


subroutine euler(mod,rot_new,a,b,c)
  implicit real*8 (a-h,o-z)
  
  integer mod
  real tmp

  dimension rot_old(3,3),rot_new(3,3),z(3),zz(3)
  dimension y(3),yy(3),yyy(3),pom(3),x(3),xx(3)
  
  zero=1.0d-20
  pi=acos(-1d0)
  
  do i=1,3
     do j=1,3
        rot_old(i,j)=0
        if (i.eq.j) rot_old(i,i)=1
     end do
  end do
  
  do j=1,3
     y(j)=rot_old(j,2)
     yyy(j)=rot_new(j,2)
     z(j)=rot_old(j,3)
     zz(j)=rot_new(j,3)
  end do
  
  call vecprod(z,zz,yy)
  y_norm=dsqrt(dot(yy,yy))
  
  if (y_norm.lt.1d-10) then
     
     if (abs(dot(y,yyy)).gt.1d0) then
        aa=dot(y,yyy)/abs(dot(y,yyy))
        a=acos(aa)
     else
        a=acos(dot(y,yyy))
     end if
     
     if (dot(z,zz).gt.zero) then
        c=zero
        b=zero
        if (yyy(1).gt.zero) a=2*pi-a
     else
        c=a
        a=zero
        b=pi
        if (yyy(1).lt.zero) c=2*pi-c
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
     if (yy(1).gt.zero) a=2*pi-a
     call vecprod(yy,yyy,pom)
     if (dot(pom,zz).lt.-zero) c=2*pi-c
  end if

 if (mod.eq.1) then
!   space fixed axes (exchange a,c)
    tmp=a
    a=c
    c=tmp
 endif 


end subroutine euler

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
 
 subroutine determinant(sym,det)
  
  real*8   sym(3,3),det

  det=sym(1,1)*(sym(2,2)*sym(3,3)-sym(3,2)*sym(2,3))-&
      sym(1,2)*(sym(2,1)*sym(3,3)-sym(3,1)*sym(2,3))+&
      sym(1,3)*(sym(2,1)*sym(3,2)-sym(2,2)*sym(3,1))

 end subroutine determinant

 subroutine find_rot_mat(l,a,b,c,mat,lmax)

  integer l,lmax
  real*8  a,b,c
  complex*16 mat(-lmax:lmax,-lmax:lmax)
  integer n,m
  complex*16 imag
  real*8  dlmat 

  imag=(0.0d0,1.0d0)

  do n=-l,l
     do m=-l,l
        call find_dlmat(l,n,m,b,dlmat) 
        mat(n,m)=exp(-imag*n*c)*dlmat*exp(-imag*m*a)
!    if(l.eq.1) write(*,*) n,m,mat(n,m),dlmat,exp(-imag*n*c),exp(-imag*m*a) 
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
 
 subroutine apply_inversion_Ylm(l,lmax,mat)

   integer l
   complex*16 mat(-lmax:lmax,-lmax:lmax)

   mat=((-1.0d0)**l)*mat
   
 end subroutine apply_inversion_Ylm



        subroutine shiftpos (r,lattic,irot0)
	IMPLICIT REAL*8 (A-H,O-Z)
        character*4 lattic
        real*8 r(3)
        if(lattic(1:1).eq.'F') then
         if(irot0.eq.0) then
         r(1)=r(1)+0.5d0
         r(2)=r(2)+0.5d0
         irot0=1
         if(r(1).gt.0.99999999999999d0) r(1)=r(1)-1.d0
         if(r(1).lt.-0.0000000000001d0) r(1)=r(1)+1.d0
         if(r(2).gt.0.99999999999999d0) r(2)=r(2)-1.d0
         if(r(2).lt.-0.0000000000001d0) r(2)=r(2)+1.d0
         else if(irot0.eq.1) then
         r(1)=r(1)-0.5d0
         r(3)=r(3)+0.5d0
         irot0=2
         if(r(1).gt.0.99999999999999d0) r(1)=r(1)-1.d0
         if(r(1).lt.-0.0000000000001d0) r(1)=r(1)+1.d0
         if(r(3).gt.0.99999999999999d0) r(3)=r(3)-1.d0
         if(r(3).lt.-0.0000000000001d0) r(3)=r(3)+1.d0
         else if(irot0.eq.2) then
         r(2)=r(2)-0.5d0
         r(1)=r(1)+0.5d0
         irot0=3
         if(r(1).gt.0.99999999999999d0) r(1)=r(1)-1.d0
         if(r(1).lt.-0.0000000000001d0) r(1)=r(1)+1.d0
         if(r(2).gt.0.99999999999999d0) r(2)=r(3)-1.d0
         if(r(2).lt.-0.0000000000001d0) r(2)=r(3)+1.d0
         else
         irot0=0
         endif
        endif        
        if(lattic(1:1).eq.'B') then
         if(irot0.eq.0) then
         r(1)=r(1)+0.5d0
         r(2)=r(2)+0.5d0
         r(3)=r(3)+0.5d0
         irot0=1
         if(r(1).gt.0.99999999999999d0) r(1)=r(1)-1.d0
         if(r(1).lt.-0.0000000000001d0) r(1)=r(1)+1.d0
         if(r(2).gt.0.99999999999999d0) r(2)=r(2)-1.d0
         if(r(2).lt.-0.0000000000001d0) r(2)=r(2)+1.d0
         if(r(3).gt.0.99999999999999d0) r(3)=r(3)-1.d0
         if(r(3).lt.-0.0000000000001d0) r(3)=r(3)+1.d0
         else
         irot0=0
         endif
        endif
        return
        end


      subroutine transform(T,Pinv,A,P)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(3,3),P(3,3),A(3,3),Pinv(3,3)
      
      do i=1,3
         do j=1,3
            sum=0
            do k=1,3
               do l=1,3
                  sum=sum+Pinv(i,k)*A(k,l)*P(l,j)
               end do
            end do
            T(i,j)=sum
         end do
!       write(6,*)'Transf. matrix:',(T(i,k),k=1,3)
      end do
      return
      end
      SUBROUTINE INVERSSYMDEF(A,AINV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(3,3),AINV(3,3)
        det= a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
            +a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3) &
            -a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
      AINV(1,1) =(   A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / det
      AINV(2,1) =( - A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / det
      AINV(3,1) =(   A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / det
      AINV(1,2) =( - A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / det
      AINV(2,2) =(   A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / det
      AINV(3,2) =( - A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / det
      AINV(1,3) =(   A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / det
      AINV(2,3) =( - A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / det
      AINV(3,3) =(   A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / det
      RETURN
      END

      SUBROUTINE LATGEN(lattic,aa,bb,cc,alph,beta,gamma,br1)
!                                                                       
!     LATGEN GENERATES TWO BRAVAIS MATRICES, DEFINES THE VOLUME OF      
!     THE UNIT CELL AND CALLS ROTDEF                                    
!     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS      
!                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN      
!                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN      
!                 WFTAPE) INTO CARTESIAN SYSTEM                         
!     BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-      
!                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )       
!                 TO CARTESIAN SYSTEM                                   
!        USE param
!        USE struct
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      character*4 lattic
      LOGICAL          ORTHO           
!                           
!      COMMON /ORTH/   ORTHO                
!      COMMON /GENER/  BR1(3,3),BR2(3,3) 
      dimension  BR1(3,3),BR2(3,3) ,alpha(3),pia(3)
                                   
!---------------------------------------------------------------------  
!                      
      PI=ACOS(-1.0D0)                                                   
      SQRT3=SQRT(3.D0)
      ALPHA(1)=ALPH*PI/180.0D0                                             
      ALPHA(2)=beta*PI/180.0D0                                             
      ALPHA(3)=gamma*PI/180.0D0                                             
      PIA(1)=2.D0*PI/AA                                                 
      PIA(2)=2.D0*PI/BB                                                 
      PIA(3)=2.D0*PI/CC                                                 
      IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
      IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'F') GOTO 30                                    
      IF(LATTIC(1:1).EQ.'B') GOTO 40                                    
      IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
      IF(LATTIC(1:1).EQ.'R') GOTO 60                                    
!
!        Error: wrong lattice, stop execution
!
         GOTO 900
!                                                                       
!.....HEXAGONAL LATTICE                                                 
 10   CONTINUE                                                          
      BR1(1,1)=2.D0/SQRT3*PIA(1)                                        
      BR1(1,2)=1.D0/SQRT3*PIA(1)                                        
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)=2.D0/SQRT3*PIA(1)                                        
      BR2(1,2)=1.D0/SQRT3*PIA(1)                                        
      BR2(1,3)=0.0D0                                                    
      BR2(2,1)=0.0D0                                                    
      BR2(2,2)=PIA(2)                                                   
      BR2(2,3)=0.0D0                                                    
      BR2(3,1)=0.0D0                                                    
      BR2(3,2)=0.0D0                                                    
      BR2(3,3)=PIA(3)                                                   
!                                                                       
      RVFAC=2.D0/SQRT(3.D0)                                             
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....RHOMBOHEDRAL CASE                                                    
 60   BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
      BR1(2,1)=-1.0d0*PIA(2)                                                  
      BR1(2,2)=1.0d0*PIA(2)                                                    
      BR1(2,3)=0.0d0*PIA(2)                                                    
      BR1(3,1)=1.0d0*PIA(3)                                                    
      BR1(3,2)=1.0d0*PIA(3)                                                    
      BR1(3,3)=1.0d0*PIA(3)                                                    
!
      BR2(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR2(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR2(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
      BR2(2,1)=-1.0d0*PIA(2)                                                   
      BR2(2,2)=1.0d0*PIA(2)                                                    
      BR2(2,3)=0.0d0*PIA(2)                                                    
      BR2(3,1)=1.0d0*PIA(3)                                                    
      BR2(3,2)=1.0d0*PIA(3)                                                    
      BR2(3,3)=1.0d0*PIA(3)                                                    
      RVFAC=6.D0/SQRT(3.D0)
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....PRIMITIVE LATTICE                                                 
!                                                                       
  20  CONTINUE
      SINBC=SIN(ALPHA(1))
      COSAB=COS(ALPHA(3))
      COSAC=COS(ALPHA(2))
      COSBC=COS(ALPHA(1))
      WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
      BR2(1,1)= SINBC/WURZEL*PIA(1)
      BR2(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
      BR2(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
      BR2(2,1)= 0.0
      BR2(2,2)= PIA(2)/SINBC
      BR2(2,3)= -PIA(3)*COSBC/SINBC
      BR2(3,1)= 0.0
      BR2(3,2)= 0.0
      BR2(3,3)= PIA(3)
!
      BR1(1,1)= SINBC/WURZEL*PIA(1)
      BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
      BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
      BR1(2,1)= 0.0
      BR1(2,2)= PIA(2)/SINBC
      BR1(2,3)= -PIA(3)*COSBC/SINBC
      BR1(3,1)= 0.0
      BR1(3,2)= 0.0
      BR1(3,3)= PIA(3)
!
      RVFAC= 1.d0/WURZEL
      ORTHO=.TRUE.
      if(abs(alpha(1)-pi/2.d0).gt.0.0001) ortho=.false.
      if(abs(alpha(2)-pi/2.d0).gt.0.0001) ortho=.false.
      if(abs(alpha(3)-pi/2.d0).gt.0.0001) ortho=.false.
!
      GOTO 100
!                                                                       
!.....FC LATTICE                                                        
 30   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                    
!     definitions according to column, rows convention for BR2
!
      BR2(1,1)=-PIA(1)                                                  
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= PIA(1)                                                  
      BR2(2,1)= PIA(2)                                                  
      BR2(2,2)=-PIA(2)                                                  
      BR2(2,3)= PIA(2)                                                  
      BR2(3,1)= PIA(3)                                                  
      BR2(3,2)= PIA(3)                                                  
      BR2(3,3)=-PIA(3)                                                  
!                                                                       
      RVFAC=4.D0                                                        
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....BC LATTICE                                                        
 40   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= 0.0D0                                                    
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= PIA(1)                                                  
      BR2(2,1)= PIA(2)                                                  
      BR2(2,2)= 0.0D0                                                    
      BR2(2,3)= PIA(2)                                                  
      BR2(3,1)= PIA(3)                                                  
      BR2(3,2)= PIA(3)                                                  
      BR2(3,3)= 0.0D0
!
      RVFAC=2.D0
      ORTHO=.TRUE.                                             
      GOTO 100                                                    
!             
 50   CONTINUE                                                          
      IF(LATTIC(2:3).EQ.'XZ') GOTO 51                                    
      IF(LATTIC(2:3).EQ.'YZ') GOTO 52                                    
!.....CXY LATTICE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= PIA(1)                                                    
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= 0.0D0                                                  
      BR2(2,1)=-PIA(2)                                                  
      BR2(2,2)= PIA(2)                                                    
      BR2(2,3)= 0.0D0                                                 
      BR2(3,1)= 0.0D0                                                  
      BR2(3,2)= 0.0D0                                                 
      BR2(3,3)= PIA(3)                                                    
!                                                                       
      RVFAC=2.D0                                                        
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....CXZ CASE (CXZ LATTICE BUILD UP)                                     
 51   CONTINUE                                     
!.....CXZ ORTHOROMBIC CASE 
      IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
         BR1(1,1)=PIA(1)                                                   
         BR1(1,2)=0.0D0                                                    
         BR1(1,3)=0.0D0                                                    
         BR1(2,1)=0.0D0                                                    
         BR1(2,2)=PIA(2)                                                   
         BR1(2,3)=0.0D0                                                    
         BR1(3,1)=0.0D0                                                    
         BR1(3,2)=0.0D0                                                    
         BR1(3,3)=PIA(3)                                                   
!                                                                       
         BR2(1,1)= PIA(1)                                                   
         BR2(1,2)= 0.0                                                   
         BR2(1,3)= PIA(1)                                                      
         BR2(2,1)= 0.0                                                      
         BR2(2,2)= PIA(2)                                                     
         BR2(2,3)= 0.0                                                     
         BR2(3,1)=-PIA(3)                                                     
         BR2(3,2)= 0.0                                                     
         BR2(3,3)= PIA(3)                                                     
!                                                                       
         RVFAC=2.0                                                         
         ORTHO=.TRUE.                                             
         GOTO 100                                                          
      ELSE
!.....CXZ MONOCLINIC CASE 
         write(*,*) '  gamma not equal 90'
         SINAB=SIN(ALPHA(3))
         COSAB=COS(ALPHA(3))
!                                                                       
         BR1(1,1)= PIA(1)/SINAB 
         BR1(1,2)= -PIA(2)*COSAB/SINAB
         BR1(1,3)= 0.0                                                   
         BR1(2,1)= 0.0                                                      
         BR1(2,2)= PIA(2)                                                     
         BR1(2,3)= 0.0                                                     
         BR1(3,1)= 0.0                                                     
         BR1(3,2)= 0.0                                                     
         BR1(3,3)= PIA(3)                                                     
!                                                                       
         BR2(1,1)= PIA(1)/SINAB 
         BR2(1,2)= -PIA(2)*COSAB/SINAB
         BR2(1,3)= PIA(1)/SINAB 
         BR2(2,1)= 0.0                                                      
         BR2(2,2)= PIA(2)                                                     
         BR2(2,3)= 0.0                                                     
         BR2(3,1)=-PIA(3)                                                     
         BR2(3,2)= 0.0                                                     
         BR2(3,3)= PIA(3)                                                     
!                                                                       
         RVFAC=2.0/SINAB                                                   
         ORTHO=.FALSE.                                             
         GOTO 100                                                          
      ENDIF
!                                                                       
!.....CYZ CASE (CYZ LATTICE BUILD UP)                                     
 52   CONTINUE                                     
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= PIA(1)                                                      
      BR2(1,2)= 0.0                                                   
      BR2(1,3)= 0.0                                                      
      BR2(2,1)= 0.0                                                      
      BR2(2,2)= PIA(2)                                                     
      BR2(2,3)= PIA(2)                                                     
      BR2(3,1)= 0.0                                                     
      BR2(3,2)=-PIA(3)                                                     
      BR2(3,3)= PIA(3)                                                     
!                                                                       
      RVFAC=2.0                                                         
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....DEFINE VOLUME OF UNIT CELL                                        
 100  CONTINUE                                                          
      write(14,*)' '
      write(14,*)' BR1,  BR2'
      do i=1,3
      write(14,654)(br1(i,j),j=1,3),(br2(i,j),j=1,3)
654   format(3f10.5,3x,3f10.5)
      enddo
      write(14,*)' '
      VOL=AA*BB*CC/RVFAC                                                
!                                                                       
!.....DEFINE ROTATION MATRICES IN NONSYMMORPHIC CASE                    
!      CALL ROTDEF                                            
!
      RETURN
!
!        Error messages
!
  900 continue
! 900  CALL OUTERR('LATGEN','wrong lattice.')
      STOP 'LATGEN - Error'
!
!        End of 'LATGEN'
!
      END                                                               









