      program phonon2wien2k
!     ########################################      
!
!     Program generates WIEN struct files + job from phonon-file.
!
!     ########################################
      implicit real*8(a-h,o-z)
      logical ortho
      character*79    fname2,prob,fname3
      
      character*4     LATTIC, IREL,CFORM, UNIT
      real*8          ROTLOC(3,3),anorm(3)

      character*79    fname
      character*80    TITLE
      character       ch
      character,allocatable::    phNAME(:)
      character*2,allocatable::    ANAME(:)
      integer,allocatable:: IATNR(:),izz(:),inum(:),ityp(:),natdispl(:),i_case(:)
      integer        mult_case(48),iz(3,3,48)
      dimension   xmat(3,3),cmat(3,3),pmat(3,3),tmat(3,3),rvec(3),tau(3,48)
      real*8, allocatable:: pos(:,:),displ(:,:),force(:,:),pos_case(:,:),force_case(:,:)
      
!     ### program input

      write(*,*) 'Program generates  Phonon-Hellman-Feynman file from WIEN calculations'
      write(*,*)
      write (*,*) 'Filename of phonon file: '
      read (*,17) fname
      write (*,*) fname
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
      pi=acos(-1.d0)
      alpha=acos(alpha)*180.d0/pi
      beta=acos(beta)*180.d0/pi
      gamma=acos(gamma)*180.d0/pi
      aa=aa/0.529177
      bb=bb/0.529177
      cc=cc/0.529177
!      write(*,*) aa,bb,cc,alpha,beta,gamma
!
        allocate (pos(3,nat),ityp(nat),force(3,nat),pos_case(3,nat),i_case(nat),force_case(3,nat))
        read(20,'(a80)',err=999) title
        read(20,'(a80)',err=999) title
        do i=1,nat
          read(20,*) (pos(j,i),j=1,3),ityp(i)
        enddo
        read(20,'(a80)',err=999) title
        read(20,*) ndispl
        allocate (displ(3,ndispl),natdispl(ndispl))
        read(20,'(a80)',err=999) title
        do i=1,ndispl
          read(20,*) (displ(j,i),j=1,3),ch,natdispl(i)
        enddo

      fname2=fname
      fname3=fname
      do i=79,2,-1
      if(fname(i:i).eq.'.') then
       fname2(1:i)=fname(1:i)
       fname2(i+1:i+4)='hff'
       fname3(1:i)=fname(1:i)
       fname3(i+1:i+8)='hff_sym'
       exit
      endif
      enddo
      open (22,FILE=fname2,STATUS='unknown',FORM='formatted',ERR=510)
      open (23,FILE=fname3,STATUS='unknown',FORM='formatted',ERR=510)
      open (21,FILE='case.finM',STATUS='old',FORM='formatted',ERR=510)
      open (66,FILE='case.out',STATUS='unknown',FORM='formatted',ERR=510)

      do J=1,3
      write(22,'(3f11.6)')   (pmat(i,j),i=1,3)
      enddo
! form invers of tmat
!      if(lattic(1:1).eq.'F') tmat=tmat/sqrt(2.d0)
!      if(lattic(1:1).eq.'B') tmat=tmat/sqrt(3.d0)
!orig      call INVERSSYMDEF(tmat,pmat)
!      call INVERSSYMDEF(pmat,xmat)
!      if(lattic(1:1).eq.'F') tmat=tmat*sqrt(2.d0)
!      if(lattic(1:1).eq.'B') tmat=tmat*sqrt(3.d0)
       fmist=sqrt(3.d0)
       if(lattic(1:1).eq.'B') fmist=sqrt(2.d0)
        fmist=1.d0
              anorm(1)=sqrt(pmat(1,1)**2+pmat(2,1)**2+pmat(3,1)**2) /fmist     
              anorm(2)=sqrt(pmat(1,2)**2+pmat(2,2)**2+pmat(3,2)**2) /fmist     
              anorm(3)=sqrt(pmat(1,3)**2+pmat(2,3)**2+pmat(3,3)**2) /fmist     

      write(66,*) 'S- matrix, Norm'
      do i=1,3
      do j=1,3
      xmat(i,j)=pmat(i,j)/anorm(j)
      enddo
      write(66,'(3f10.6,3x,f10.6)')   (pmat(i,j),j=1,3),anorm(i)
      enddo
      write(66,*) 'normalized S-matrix'
      do i=1,3
      write(66,'(3f10.6,3x,f10.6)')   (xmat(i,j),j=1,3)
      enddo
!
      f1tot=0.d0
      f2tot=0.d0
      f3tot=0.d0
! write all cases
   DO idispl=1,ndispl
      write(22,'(7x,3f11.7,3x,3f12.7)') (POS(k,natdispl(idispl)),k=1,3),(displ(k,idispl),k=1,3)
! read case.finM for case idispl
      Read(21,*) nat_case
!      read(21,'(15x,i2,/)') (mult_case(i),i=1,nat_case)
      do i=1,nat
        read(21,'(5x,i3,4x,f10.8,3x,f10.8,3x,f10.8)') i_case(i),(pos_case(j,i),j=1,3)
      enddo
      do i=1,nat_case
         read(21,*) (force_case(k,i),k=1,3)
      enddo
      read(21,*) nsym_case
      do i=1,nsym_case
        read(21,'(3i2,f10.7,/,3i2,f10.7,/,3i2,f10.7,/i8)') ((iz(j,k,i),k=1,3),tau(j,i),j=1,3),inum0
!         write(*,*) 'inum',inum0
      enddo
!
      f1sum=0.d0
      f2sum=0.d0
      f3sum=0.d0
!
      do i=1,nat
        if(natdispl(idispl).eq.i) then
           p1=POS(1,i)+displ(1,idispl)
           p2=POS(2,i)+displ(2,idispl)
           p3=POS(3,i)+displ(3,idispl)
           if(p1.ge.1.d0) p1=p1-1.d0
           if(p1.lt.0.d0) p1=p1+1.d0
           if(p2.ge.1.d0) p2=p2-1.d0
           if(p2.lt.0.d0) p2=p2+1.d0
           if(p3.ge.1.d0) p3=p3-1.d0
           if(p3.lt.0.d0) p3=p3+1.d0
        else
          p1=pos(1,i)
          p2=pos(2,i)
          p3=pos(3,i)
        endif
        if(lattic(1:1).eq.'P')  then
        else if(lattic(1:1).eq.'F') then
          p1a=(p1*tmat(1,1)+p2*tmat(1,2)+p3*tmat(1,3))*tmatnorm
          p2a=(p1*tmat(2,1)+p2*tmat(2,2)+p3*tmat(2,3))*tmatnorm
          p3a=(p1*tmat(3,1)+p2*tmat(3,2)+p3*tmat(3,3))*tmatnorm
          call pnorm(p1a,p2a,p3a)
          p1=p1a
          p2=p2a
          p3=p3a
        else if(lattic(1:1).eq.'B') then
          p1a=(p1*tmat(1,1)+p2*tmat(1,2)+p3*tmat(1,3))*tmatnorm
          p2a=(p1*tmat(2,1)+p2*tmat(2,2)+p3*tmat(2,3))*tmatnorm
          p3a=(p1*tmat(3,1)+p2*tmat(3,2)+p3*tmat(3,3))*tmatnorm
          call pnorm(p1a,p2a,p3a)
          p1=p1a
          p2=p2a
          p3=p3a
        endif
! compare p1,2,3 with pos_case and rotate if necessary
        do j=1,nat
        write(66,'("ATOM",i3,3f10.4,3x,i3,3f10.4)') i,p1,p2,p3,j,pos_case(1,j),pos_case(2,j),pos_case(3,j)
          if(abs(p1-pos_case(1,j)).lt.0.00001d0.and.abs(p2-pos_case(2,j)).lt.0.00001d0.and.abs(p3-pos_case(3,j)).lt.0.00001d0) then
            index=iabs(i_case(j))
            do j1=j,1,-1
              if(i_case(j1).ne.i_case(j)) goto 50
            enddo
            j1=0
 50         j1=j1+1
! find rotation from j1 to j
            do k=1,nsym_case
              X=0.D0                                                      
              Y=0.D0                                                      
              Z=0.D0                                                      
              DO 26 J2=1,3                                                 
              Z=Z+IZ(3,J2,k)*POS_case(J2,j1)                      
              Y=Y+IZ(2,J2,k)*POS_case(J2,j1)                               
 26           X=X+IZ(1,J2,k)*POS_case(J2,j1)
              Z=Z+ tau(3,k) 
              Y=Y+ tau(2,k) 
              X=X+ tau(1,k) 
              irot0=0
 144          continue
              call pnorm(x,y,z)
              write(66,'("rot",i3," to",i3," op",i3,3f10.5)') j1,j,k,x,y,z
              if(abs(x-POS_case(1,j)).lt.1.d-5.and. &
                 abs(y-POS_case(2,j)).lt.1.d-5.and. &
                 abs(z-POS_case(3,j)).lt.1.d-5) then
               force(1,i)=0.d0
               force(2,i)=0.d0
               force(3,i)=0.d0
               DO 27 J2=1,3                                                 
               force(1,i)=force(1,i)+force_case(j2,iabs(i_case(j1)))*iz(1,j2,k)
               force(2,i)=force(2,i)+force_case(j2,iabs(i_case(j1)))*iz(2,j2,k)
               force(3,i)=force(3,i)+force_case(j2,iabs(i_case(j1)))*iz(3,j2,k)
 27            continue
               goto 54
              else
               rvec(1)=x
               rvec(2)=y
               rvec(3)=z
               call shiftpos(rvec,lattic,irot0)
               x=rvec(1)
               y=rvec(2)
               z=rvec(3)
               if(irot0.ne.0) goto 144
              endif
            enddo
            stop 'rotation not found'
 54         continue
            goto 55
          endif
        enddo
        stop 'pos_case not found'
 55     continue
! multiply force with inv bravais mat
              f1=0.d0
              f2=0.d0
              f3=0.d0
              if(lattic.eq.'F'.or.lattic.eq.'B') then
              f1=force(1,i)
              f2=force(2,i)
              f3=force(3,i)
              else
              DO 28 J2=1,3                                                 
!              f1=f1+force(j2,i)*pmat(1,j2)/anorm(1)
 !             f2=f2+force(j2,i)*pmat(2,j2)/anorm(2)
  !            f3=f3+force(j2,i)*pmat(3,j2)/anorm(3)
              f1=f1+force(j2,i)*xmat(1,j2)
              f2=f2+force(j2,i)*xmat(2,j2)
              f3=f3+force(j2,i)*xmat(3,j2)
 28           continue
              endif
!           do k=1,3
!             force(k,i)=force(k,i)*13.6d0/1000.d0/0.529177d0
!           enddo
          write(22,'(i3,1x,a1,2x,3f10.7,3x,3f10.4)') i,phname(ityp(i)),(POS(k,i),k=1,3),f1,f2,f3
      f1sum=f1sum+f1
      f2sum=f2sum+f2
      f3sum=f3sum+f3
!
      enddo
      f1tot=f1tot+abs(f1sum)
      f2tot=f2tot+abs(f2sum)
      f3tot=f3tot+abs(f3sum)
      write(66,'("Sum of forces (should be zero) for case",i4,":",3f10.4)') idispl,f1sum,f2sum,f3sum 
      write(*,'("Sum of forces (should be zero) for case",i4,":",3f10.4)') idispl,f1sum,f2sum,f3sum 
   enddo
      write(66,'("Sum of forces of all displacements:",3f10.4)') f1tot,f2tot,f3tot 
      write(*,'("Sum of forces of all displacements:",3f10.4)') f1tot,f2tot,f3tot 
      call sym_for(ndispl,nat,fname3)
      stop


!     ### errors leading to termination of program
  500 prob='Error accessing file: ' // fname
      call faterr(prob)
      stop
  510 prob='Error accessing file: ' // fname2
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
 2030 FORMAT('ATOM=',I3,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,10X,'MULT=',I2,10X,'ISPLIT=',I2)
 2031 FORMAT('ATOM=',I3,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)
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
         if(r(2).gt.0.99999999999999d0) r(2)=r(2)-1.d0
         if(r(2).lt.-0.0000000000001d0) r(2)=r(2)+1.d0
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
         call pnorm(r(1),r(2),r(3))
         else
         irot0=0
         endif
        endif
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

      subroutine sym_for(ndispl,nat,fname3)
      implicit real*8(a-h,o-z)
      real*8, allocatable::  pos1(:,:),displ(:,:),pos(:,:),force(:,:,:)
      character*1, allocatable:: ilabel(:)
      character*79    fname3
      real*8 a(3,3)

      allocate (pos1(3,ndispl),displ(3,ndispl),pos(3,nat),force(3,nat,ndispl))
      allocate(ilabel(nat))

      rewind 22
      read(22,*) a
      write(23,'(3f11.6)') a

      do idispl=1,ndispl
        read(22,*) pos1(1,idispl),pos1(2,idispl),pos1(3,idispl),displ(1,idispl),displ(2,idispl),displ(3,idispl)
        do i=1,nat
         read(22,'(i3,1x,a1,2x,3f10.7,3x,3f10.4)') i1,ilabel(i),pos(1,i),pos(2,i),pos(3,i),force(1,i,idispl),force(2,i,idispl),force(3,i,idispl)
        enddo
      enddo
! all of file 22 read
!
      ftotx=0.d0
      ftoty=0.d0
      ftotz=0.d0
      do idispl=1,ndispl-1,2
        ftx=0.d0
        fty=0.d0
        ftz=0.d0
        ftxp=0.d0
        ftyp=0.d0
        ftzp=0.d0
        ftxm=0.d0
        ftym=0.d0
        ftzm=0.d0
        if(pos1(1,idispl).eq.pos1(1,idispl+1).and. &
           pos1(2,idispl).eq.pos1(2,idispl+1).and. &
           pos1(3,idispl).eq.pos1(3,idispl+1)) then
           if(displ(1,idispl).eq.-displ(1,idispl+1).and. &
              displ(2,idispl).eq.-displ(2,idispl+1).and. &
              displ(3,idispl).eq.-displ(3,idispl+1)) then
              write(66,*) 'displacements',idispl,' and',idispl+1,' can be symmetrized'
              write(*,*) 'displacements',idispl,' and',idispl+1,' can be symmetrized'
              write(66,*) 'The following forces differ by more than 20 %:'
              write(*,*) 'The following forces differ by more than 20 %:'
              write(23,'(8x,3f10.7,3x,3f12.7)') pos1(1,idispl),pos1(2,idispl), &
                 pos1(3,idispl),displ(1,idispl),displ(2,idispl),displ(3,idispl)
              do i=1,nat
!
!check quality (forces differ by more than 20 percent)
              iwrit=0
              dfx=(abs(force(1,i,idispl))-abs(force(1,i,idispl+1)))/abs(force(1,i,idispl))
              dfy=(abs(force(2,i,idispl))-abs(force(2,i,idispl+1)))/abs(force(2,i,idispl))
              dfz=(abs(force(3,i,idispl))-abs(force(3,i,idispl+1)))/abs(force(3,i,idispl))
              if(abs(force(1,i,idispl)).gt.0.5d0 .and. dfx.gt.0.2d0) iwrit=1
              if(abs(force(2,i,idispl)).gt.0.5d0 .and. dfy.gt.0.2d0) iwrit=1
              if(abs(force(3,i,idispl)).gt.0.5d0 .and. dfz.gt.0.2d0) iwrit=1
              dfx=(abs(force(1,i,idispl))-abs(force(1,i,idispl+1)))/abs(force(1,i,idispl+1))
              dfy=(abs(force(2,i,idispl))-abs(force(2,i,idispl+1)))/abs(force(2,i,idispl+1))
              dfz=(abs(force(3,i,idispl))-abs(force(3,i,idispl+1)))/abs(force(3,i,idispl+1))
              if(abs(force(1,i,idispl+1)).gt.0.5d0 .and. dfx.gt.0.2d0) iwrit=1
              if(abs(force(2,i,idispl+1)).gt.0.5d0 .and. dfy.gt.0.2d0) iwrit=1
              if(abs(force(3,i,idispl+1)).gt.0.5d0 .and. dfz.gt.0.2d0) iwrit=1
!               
              if(iwrit.eq.1) then
                  write(66,'(" atom ",i3,3f10.3)') &
                        i,force(1,i,idispl),force(2,i,idispl),force(3,i,idispl)
                  write(66,'("      ",i3,3f10.3)') &
                        i,force(1,i,idispl+1),force(2,i,idispl+1),force(3,i,idispl+1)
                  write(*,'(" atom ",i3,3f10.3)') &
                        i,force(1,i,idispl),force(2,i,idispl),force(3,i,idispl)
                  write(*,'("      ",i3,3f10.3)') &
                        i,force(1,i,idispl+1),force(2,i,idispl+1),force(3,i,idispl+1)
              endif
!

              fx=(force(1,i,idispl)-force(1,i,idispl+1))/2.d0
              fy=(force(2,i,idispl)-force(2,i,idispl+1))/2.d0
              fz=(force(3,i,idispl)-force(3,i,idispl+1))/2.d0
              if(fx.gt.0.d0) ftxp=ftxp+fx  
              if(fy.gt.0.d0) ftyp=ftyp+fy  
              if(fz.gt.0.d0) ftzp=ftzp+fz  
              if(fx.lt.0.d0) ftxm=ftxm+fx  
              if(fy.lt.0.d0) ftym=ftym+fy  
              if(fz.lt.0.d0) ftzm=ftzm+fz  
              enddo
              renxp=(ftxp-ftxm)/2.d0/ftxp
              renxm=(ftxm-ftxp)/2.d0/ftxm
              renyp=(ftyp-ftym)/2.d0/ftyp
              renym=(ftym-ftyp)/2.d0/ftym
              renzp=(ftzp-ftzm)/2.d0/ftzp
              renzm=(ftzm-ftzp)/2.d0/ftzm
              do i=1,nat
              fx=(force(1,i,idispl)-force(1,i,idispl+1))/2.d0
              fy=(force(2,i,idispl)-force(2,i,idispl+1))/2.d0
              fz=(force(3,i,idispl)-force(3,i,idispl+1))/2.d0
              if(fx.gt.0.d0) fx=fx*renxp 
              if(fy.gt.0.d0) fy=fy*renyp  
              if(fz.gt.0.d0) fz=fz*renzp  
              if(fx.lt.0.d0) fx=fx*renxm   
              if(fy.lt.0.d0) fy=fy*renym  
              if(fz.lt.0.d0) fz=fz*renzm   
              ftx=ftx+fx  
              fty=fty+fy  
              ftz=ftz+fz  
              write(23,'(i3,1x,a1,2x,3f10.7,3x,3f10.4)') i,ilabel(i),pos(1,i),pos(2,i),pos(3,i),fx,fy,fz
              enddo
      write(66,'("Sum of forces (symmetrized) for case",i4,":",3f10.4)') idispl,ftx,fty,ftz 
!      write(*,'("Sum of forces (symmetrized) for case",i4,":",3f10.4)') idispl,ftx,fty,ftz 
           endif
        endif
        ftotx=ftotx+abs(ftx)
        ftoty=ftoty+abs(fty)
        ftotz=ftotz+abs(ftz)
      enddo
      write(66,'("Sum of forces of all symm.displ:   ",3f10.4)') ftotx,ftoty,ftotz 
      write(*,'("Sum of forces of all symm.displ:   ",3f10.4)') ftotx,ftoty,ftotz 
      return
      end
