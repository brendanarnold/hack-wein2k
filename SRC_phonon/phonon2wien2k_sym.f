      program phonon2wien2k
!     ########################################      
!
!     Program generates WIEN struct files + job from phonon-file.
!
!     ########################################
      implicit real*8(a-h,o-z)
      character*79    fname2,prob
      
      character*4     LATTIC, IREL,CFORM, UNIT
      real*8          ROTLOC(3,3)

      character*79    fname
      character*80    TITLE
      character       ch
      character,allocatable::    phNAME(:)
      character*2,allocatable::    ANAME(:)
      integer,allocatable::        IATNR(:),izz(:),inum(:),ityp(:),natdispl(:)
      dimension       cmat(3,3),pmat(3,3),tmat(3,3)
      real*8, allocatable:: pos(:,:),displ(:,:)
      logical ortho      
!     ### program input
      
      write(*,*) 'Program generates  WIEN struct files from phonon file.'
      write(*,*)
      write (*,*) 'Filename of phonon file: '
      read (*,17) fname
   17 format(A79)

      open (20,FILE=fname,STATUS='old',FORM='formatted',ERR=500)
 100  read(20,'(a80)',end=999) title
      if(title(3:7).ne.'ATOMS') goto 100
        read(20,*,err=999) nequ_at
        read(20,*,err=999) nat
        allocate (aname(nequ_at),phname(nequ_at),izz(nequ_at),iatnr(nequ_at))
        allocate (inum(nequ_at))
        read(20,'(a80)',err=999) title
        do i=1,nequ_at
          read(20,*) iatnr(i),izz(i),aname(i)
        enddo
        read(20,'(a80)',err=999) title
        do i=1,nequ_at
          read(20,*) iatnr(i),phname(i),aname(i)
        enddo
        read(20,'(a80)',err=999) title
        read(20,*) (inum(i),i=1,nequ_at)
        read(20,'(a80)',err=999) title
        read(20,*) lattic
        read(20,'(a80)',err=999) title
        do i=1,3
          read(20,*) (cmat(i,j),j=1,3)
        enddo
        read(20,'(a80)',err=999) title
        do i=1,3
          read(20,*) (pmat(i,j),j=1,3)
        enddo
        read(20,'(a80)',err=999) title
        do i=1,3
          read(20,*) (tmat(i,j),j=1,3)
        enddo
!
        ortho=.true.
        if(cmat(1,2).ne.0.d0) ortho=.false.
        if(cmat(1,3).ne.0.d0) ortho=.false.
        if(cmat(2,1).ne.0.d0) ortho=.false.
        if(cmat(2,3).ne.0.d0) ortho=.false.
        if(cmat(3,1).ne.0.d0) ortho=.false.
        if(cmat(3,2).ne.0.d0) ortho=.false.
        lattic='P'
        if(ortho) then
! check for FC lattice
          if(tmat(1,1).eq.tmat(1,3).and.tmat(1,2).eq.0.d0) then
           tmatnorm=0.5d0/tmat(1,1)
           if(tmat(2,1).eq.tmat(2,2).and.tmat(2,3).eq.0.d0) then
             if(tmat(3,2).eq.tmat(3,3).and.tmat(3,1).eq.0.d0) then
             print*, 'FC-lattic 1'
             Lattic='F'
             endif
           else if(tmat(2,2).eq.tmat(2,3).and.tmat(2,1).eq.0.d0) then
             if(tmat(3,1).eq.tmat(3,2).and.tmat(3,3).eq.0.d0) then
             print*, 'FC-lattic 2'
             Lattic='F'
             endif
           endif
          else if(tmat(1,1).eq.tmat(1,2).and.tmat(1,3).eq.0.d0) then
           tmatnorm=0.5d0/tmat(1,1)
           if(tmat(2,1).eq.tmat(2,3).and.tmat(2,2).eq.0.d0) then
             if(tmat(3,2).eq.tmat(3,3).and.tmat(3,1).eq.0.d0) then
             print*, 'FC-lattic 3'
             Lattic='F'
             endif
           else if(tmat(2,2).eq.tmat(2,3).and.tmat(2,1).eq.0.d0) then
             if(tmat(3,1).eq.tmat(3,3).and.tmat(3,2).eq.0.d0) then
             print*, 'FC-lattic 4'
             Lattic='F'
             endif
           endif
          else if(tmat(1,2).eq.tmat(1,3).and.tmat(1,1).eq.0.d0) then
           tmatnorm=0.5d0/tmat(1,2)
           if(tmat(2,1).eq.tmat(2,2).and.tmat(2,3).eq.0.d0) then
             if(tmat(3,1).eq.tmat(3,3).and.tmat(3,2).eq.0.d0) then
             print*, 'FC-lattic 5'
             Lattic='F'
             endif
           else if(tmat(2,1).eq.tmat(2,3).and.tmat(2,2).eq.0.d0) then
             if(tmat(3,1).eq.tmat(3,2).and.tmat(3,3).eq.0.d0) then
             print*, 'FC-lattic 6'
             Lattic='F'
             endif
           endif
! check for BCC lattic            
          else if(tmat(1,1).eq.tmat(1,2).and.tmat(1,1).eq.-tmat(1,3)) then
           tmatnorm=0.5d0/abs(tmat(1,1))
           if(tmat(2,1).eq.-tmat(2,2).and.tmat(2,3).eq.tmat(2,1)) then
             if(tmat(3,2).eq.tmat(3,3).and.-tmat(3,1).eq.tmat(3,3)) then
             print*, 'BC-lattic 1'
             Lattic='B'
             endif
           else if(-tmat(2,1).eq.tmat(2,2).and.tmat(2,2).eq.tmat(2,3)) then
             if(tmat(3,1).eq.-tmat(3,2).and.tmat(3,3).eq.tmat(3,1)) then
             print*, 'BC-lattic 2'
             Lattic='B'
             endif
           endif
          else if(tmat(1,1).eq.-tmat(1,2).and.tmat(1,3).eq.tmat(1,1)) then
           tmatnorm=0.5d0/abs(tmat(1,1))
           if(-tmat(2,1).eq.tmat(2,2).and.tmat(2,2).eq.tmat(2,3)) then
             if(tmat(3,2).eq.-tmat(3,3).and.tmat(3,1).eq.tmat(3,2)) then
             print*, 'BC-lattic 3'
             Lattic='B'
             endif
           else if(tmat(2,2).eq.-tmat(2,3).and.tmat(2,1).eq.tmat(2,2)) then
             if(-tmat(3,1).eq.tmat(3,3).and.tmat(3,2).eq.tmat(3,3)) then
             print*, 'BC-lattic 4'
             Lattic='B'
             endif
           endif
          else if(tmat(1,2).eq.tmat(1,3).and.-tmat(1,1).eq.tmat(1,2)) then
           tmatnorm=0.5d0/abs(tmat(1,1))
           if(tmat(2,1).eq.-tmat(2,2).and.tmat(2,3).eq.tmat(2,1)) then
             if(tmat(3,1).eq.-tmat(3,3).and.tmat(3,2).eq.tmat(3,1)) then
             print*, 'BC-lattic 5'
             Lattic='B'
             endif
           else if(tmat(2,1).eq.-tmat(2,3).and.tmat(2,2).eq.tmat(2,1)) then
             if(tmat(3,1).eq.-tmat(3,2).and.tmat(3,3).eq.tmat(3,1)) then
             print*, 'BC-lattic 6'
             Lattic='B'
             endif
           endif
          endif
        endif   
!
        pmat=0
        do i=1,3
          do j=1,3
            do k=1,3
            pmat(i,j)=pmat(i,j)+cmat(i,k)*tmat(k,j)
            enddo
          enddo
        enddo
      pi=acos(-1.d0)
!      aa=sqrt(pmat(1,1)**2+pmat(1,2)**2+pmat(1,3)**2)
!      bb=sqrt(pmat(2,1)**2+pmat(2,2)**2+pmat(2,3)**2)
!      cc=sqrt(pmat(3,1)**2+pmat(3,2)**2+pmat(3,3)**2)
      aa=sqrt(pmat(1,1)**2+pmat(2,1)**2+pmat(3,1)**2)
      bb=sqrt(pmat(1,2)**2+pmat(2,2)**2+pmat(3,2)**2)
      cc=sqrt(pmat(1,3)**2+pmat(2,3)**2+pmat(3,3)**2)
    if(lattic(1:1).eq.'P') then
      alpha=(pmat(1,2)*pmat(1,3)+pmat(2,2)*pmat(2,3)+pmat(3,2)*pmat(3,3))/bb/cc
      beta=(pmat(1,1)*pmat(1,3)+pmat(2,1)*pmat(2,3)+pmat(3,1)*pmat(3,3))/aa/cc
      gamma=(pmat(1,2)*pmat(1,1)+pmat(2,2)*pmat(2,1)+pmat(3,2)*pmat(3,1))/bb/aa
    else if(lattic(1:1).eq.'F') then
!      aa=tmat(1,1)*cmat(1,1)+tmat(1,2)*cmat(2,2)+tmat(1,3)*cmat(3,3)
!      bb=tmat(2,1)*cmat(1,1)+tmat(2,2)*cmat(2,2)+tmat(2,3)*cmat(3,3)
!      cc=tmat(3,1)*cmat(1,1)+tmat(3,2)*cmat(2,2)+tmat(3,3)*cmat(3,3)
      aa=aa*2.d0/sqrt(2.d0)
      bb=bb*2.d0/sqrt(2.d0)
      cc=cc*2.d0/sqrt(2.d0)
      alpha=0.d0
      beta=0.d0
      gamma=0.d0
    else if(lattic(1:1).eq.'B') then
      aa=aa*2.d0/sqrt(3.d0)
      bb=bb*2.d0/sqrt(3.d0)
      cc=cc*2.d0/sqrt(3.d0)
      alpha=0.d0
      beta=0.d0
      gamma=0.d0
    endif
      alpha=acos(alpha)*180.d0/pi
      beta=acos(beta)*180.d0/pi
      gamma=acos(gamma)*180.d0/pi
      aa=aa/0.529177
      bb=bb/0.529177
      cc=cc/0.529177
      write(*,*) aa,bb,cc,alpha,beta,gamma
!
        allocate (pos(3,nat),ityp(nat))
        read(20,'(a80)',err=999) title
        read(20,'(a80)',err=999) title
        do i=1,nat
          read(20,*) (pos(j,i),j=1,3),ityp(i)
        enddo
        read(20,'(a80)',err=999) title
        read(20,*) ndispl
        allocate (displ(3,ndispl),natdispl(0:ndispl))
        natdispl=-1
        read(20,'(a80)',err=999) title
        do i=1,ndispl
          read(20,*) (displ(j,i),j=1,3),ch,natdispl(i)
        enddo

      cform=''
      irel='RELA'
      unit='ang'
      jri=781
      rmt=2.0d0
      rotloc=0.d0
      rotloc(1,1)=1.d0
      rotloc(2,2)=1.d0
      rotloc(3,3)=1.d0
!
   do idispl=0,ndispl
      if(idispl.lt.10) write(fname2,'(a5,i1,a7)') 'case_',idispl,'.struct'   
      if(idispl.gt.9) write(fname2,'(a5,i2,a7)') 'case_',idispl,'.struct'   
      write(*,*) 'generating ',fname2(1:60)
      write(*,*) '     original positions               WIEN2k positions)'
      open (22,FILE=fname2,STATUS='unknown',FORM='formatted',ERR=500)
      write(22,*) 'Phonon-file by phonon2wien2k from ',fname(1:40)
      write(22,2010) lattic,NAT,CFORM,IREL,UNIT
      write(22,2020) AA,BB,CC,ALPHA,beta,gamma
      do i=1,nat
        if(natdispl(idispl).eq.i) then
           p1=POS(1,i)+displ(1,idispl)
           p2=POS(2,i)+displ(2,idispl)
           p3=POS(3,i)+displ(3,idispl)
           call pnorm(p1,p2,p3)
        else
          p1=pos(1,i)
          p2=pos(2,i)
          p3=pos(3,i)
        endif
        if(lattic(1:1).eq.'P')  then
          write(22,2030) -i,p1,p2,p3,1,8
        else if(lattic(1:1).eq.'F') then
          p1a=(p1*tmat(1,1)+p2*tmat(1,2)+p3*tmat(1,3))*tmatnorm
          p2a=(p1*tmat(2,1)+p2*tmat(2,2)+p3*tmat(2,3))*tmatnorm
          p3a=(p1*tmat(3,1)+p2*tmat(3,2)+p3*tmat(3,3))*tmatnorm
           call pnorm(p1a,p2a,p3a)
  write(*,'(i2,3f10.5,3x,3f10.5)') i,p1,p2,p3,p1a,p2a,p3a
          write(22,2030) -i,p1a,p2a,p3a,1,8
        else if(lattic(1:1).eq.'B') then
          p1a=(p1*tmat(1,1)+p2*tmat(1,2)+p3*tmat(1,3))*tmatnorm
          p2a=(p1*tmat(2,1)+p2*tmat(2,2)+p3*tmat(2,3))*tmatnorm
          p3a=(p1*tmat(3,1)+p2*tmat(3,2)+p3*tmat(3,3))*tmatnorm
           call pnorm(p1a,p2a,p3a)
  write(*,'(i2,3f10.5,3x,3f10.5)') i,p1,p2,p3,p1a,p2a,p3a
          write(22,2030) -i,p1a,p2a,p3a,1,8
        endif
         r0=0.000005d0
         if(izz(ityp(i)).le.71)   r0=0.00001d0
         if(izz(ityp(i)).le.36)   r0=0.00005d0
         if(izz(ityp(i)).le.18)   r0=0.0001d0
        prob=aname(ityp(i))
        if(prob(2:2).eq.' ') then
         write(22,3050) ANAME(ityp(i)),JRI,R0,RMT,dble(iZZ(ityp(i))) 	
        else
         write(22,2050) ANAME(ityp(i)),JRI,R0,RMT,dble(iZZ(ityp(i))) 	
        endif
        write(22,2051) ((ROTLOC(k,j),k=1,3),j=1,3)
      enddo
      
!     ### set symmetry operations to zero
      write(22,2060) 0
      write(22,'(3i2,f10.7)') 1, 0, 0, 0.0000000
      write(22,'(3i2,f10.7)') 0, 1, 0, 0.0000000
      write(22,'(3i2,f10.7)') 0, 0, 1, 0.0000000
      write(22,'(i8)')       1
      
      close (22)
   enddo

      stop


!     ### errors leading to termination of program
  500 prob='Error accessing file: ' // fname
      call faterr(prob)
      stop
  999 prob='Error reading file: ' // fname
      call faterr(prob)
      stop


 1000 FORMAT(A80)
 1010 FORMAT(A4,24X,I2,1X,A4,/,13X,A4,6X,A4)
 2010 FORMAT(A4,'LATTICE,NONEQUIV. ATOMS',I3,1X,A4,/,'MODE OF CALC=',A4,' unit=',A4)
 1020 FORMAT(6F10.6)
 2020 FORMAT(6F10.6)
 1030 FORMAT(5X,I3,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)
 1031 FORMAT(5X,I3,4X,F10.8,3X,F10.8,3X,F10.8)
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 1051 FORMAT(20X,3F10.8)
 2060 FORMAT(I4,6X,'NUMBER OF SYMMETRY OPERATIONS')
 2030 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,10X,'MULT=',I2,10X,'ISPLIT=',I2)
 2031 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)
 2050 FORMAT(A2,9X,'NPT=',I5,2X,'R0=',F10.8,1X,'RMT=',F10.4,3X,'Z:',F5.1)
 3050 FORMAT(A1,10X,'NPT=',I5,2X,'R0=',F10.8,1X,'RMT=',F10.4,3X,'Z:',F5.1)
 2051 FORMAT('LOCAL ROT MATRIX:',3X,3F10.7,/,20X,3F10.7,/,20X,3F10.7)
 
      end

      subroutine faterr(prob)
        character*80 prob
        write(*,*)
        write(*,*) 'Fatal Error occured:'
        write(*,*) prob
        write(*,*) 'Program terminated.'
        write(*,*)
        return
      end

      subroutine pnorm(p1,p2,p3)
      real*8 p1,p2,p3
           if(p1.ge.1.d0) p1=p1-1.d0
           if(p1.lt.0.d0) p1=p1+1.d0
           if(p2.ge.1.d0) p2=p2-1.d0
           if(p2.lt.0.d0) p2=p2+1.d0
           if(p3.ge.1.d0) p3=p3-1.d0
           if(p3.lt.0.d0) p3=p3+1.d0
      return
      end
