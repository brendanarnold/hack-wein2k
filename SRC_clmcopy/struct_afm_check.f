      PROGRAM struct_afm_check                                                     
!                                                                       
      use struct
      use symetr
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*11   STATUS,FORM                                        
      CHARACTER*77   FNAME
!
      integer,allocatable :: jatch(:,:),lmchm(:)   ,lsym(:,:,:)          
      real*8,allocatable :: symmat(:,:,:),taul(:,:)
!
      integer msym(3,3)
      real*8  taum(3),rrott(3),delta(3),rrot1(3)
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
      STOP 'struct_afm_check.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)

!.....READ TAPE20=POTE                                                  
!                                                                       
      call readstruct

! reading input parameters : changing atoms
      read (5,*) jmax
      allocate (jatch(jmax,2),symmat(3,3,jmax),lmchm(jmax),lsym(3,3,jmax),taul(3,jmax))
      do 20 j=1,jmax
        read(5,*) (jatch(j,jx),jx=1,2) 
        read(5,*) (symmat(1,i,j),i=1,3),(lsym(1,i,j),i=1,3),taul(1,l)
        read(5,*) (symmat(2,i,j),i=1,3),(lsym(2,i,j),i=1,3),taul(2,l)
        read(5,*) (symmat(3,i,j),i=1,3),(lsym(3,i,j),i=1,3),taul(3,j)
        read(5,*) lmchm (j)
           do iclm=1,lmchm (j)
!                    old-L,M   new-L,M  factor
           read(5,*) lmch
           enddo
 20   continue 
      do j2=1,3
      READ(5,*) (mSYM(J1,J2),J1=1,3),TAUm(J2)
      enddo
      do 32 j=1,jmax
        write(6,3000) j,(jatch(j,jx),jx=1,2) 
 3000 format(i2,' case: atom',i3,'   interchanged with',i3)
        if(jatch(j,1).eq.jatch(j,2)) goto 32
        index1=1
        index2=1
        do jatom=1,jatch(j,1)-1
          index1=index1+mult(jatom)
        enddo
        do jatom=1,jatch(j,2)-1
          index2=index2+mult(jatom)
        enddo
        do 31 ind=0,mult(jatch(j,1))-1
         do J1=1,3
          rrott(J1)=msym(j1,1)*pos(1,index1+ind)+msym(j1,2)*pos(2,index1+ind)+msym(j1,3)*pos(3,index1+ind)+TAUm(J1)
          if(rrott(j1).gt.0.99999999999999d0) rrott(j1)=rrott(j1)-1.d0
          if(rrott(j1).lt.-0.0000000000001d0) rrott(j1)=rrott(j1)+1.d0
         ENDDO
                DO  J1=1,3
                rrot1(J1)=lsym(j1,1,j)*rrott(1)+lsym(j1,2,j)*rrott(2)+lsym(j1,3,j)*rrott(3)+TAUL(J1,j)
                if(rrot1(j1).gt.0.99999999999999d0) rrot1(j1)=rrot1(j1)-1.d0
                if(rrot1(j1).lt.-0.0000000000001d0) rrot1(j1)=rrot1(j1)+1.d0
                delta(j1)=rrot1(j1)-pos(j1,index2+ind)
                if(delta(j1).gt.0.99d0) delta(j1)=delta(j1)-1.d0
                if(delta(j1).lt.-0.99d0) delta(j1)=delta(j1)+1.d0
                ENDDO

         do 30 j1=1,3
          if(abs(delta(j1)).lt.0.0000001d0) goto 30
          write(6,*) 'Position not stricktly AFM by symmetry. Symmetry enforced and positions changed'
          write(6,'(" from",3f11.8," to",3f11.8)') (pos(j11,index2+ind),j11=1,3),(rrot1(j11),j11=1,3) 
          pos(1:3,index2+ind)=rrot1(1:3)
          goto 31
 30      continue 
 31     continue
 32   continue
      call writestruct
!
      STOP 'struct_afm_check END'                                                  
      END
