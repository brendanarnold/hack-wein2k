      implicit real*8 (a-h,o-z)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
!     Main program
      include 'params.inc'
      dimension pos(3,48*natmax),pred(3,natmax),ptmp(3,natmax)
      dimension hess(3*natmax,3*natmax),hold(3*natmax,3*natmax)
      dimension hred(3*natmax,3*natmax)
      character*8 diagsval(3,natmax)
      dimension work(natmax*natmax),jpvt(natmax*3)
      dimension br1(3,3),br2(3,3),steps(3)
      integer   translat(4,neigmax,natmax*48),nuse(natmax)
      dimension distances(neigmax,natmax*48),htmp(3,3,natmax,natmax)
      dimension heigen(3*natmax,3*natmax)
      dimension W1(3*natmax),W2(3*natmax)
      dimension dmatrix(3,3,natmax),dmat2(3,3,natmax),tmp(3,3)
      character*80 fname,ename
      character*11 status,form
      character*67 ERRMSG
!
!     Name
      call gtfnam(fname,ename)
      CALL ERRFLG(Ename,'Error in PAIRHESS')
      OPEN (1,FILE=fname,STATUS='OLD',ERR=910)
8000  CONTINUE
      READ (1,*,END=8001,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
      GOTO 8000
8001  continue
!     Parameters from unit 5
!     Maximum distance in a.u.
      rmax=10.0D0
!     Scale term for exponential decay
      scl=2.0D0
!     Rescaling of Hessian estimate      
      rscale=0.25D0
!     How much to add to diagonal elements
      diagscale=1.0D0
!     Cutoff to ignore neighbors with < this weighting
      cutoff=0.05D0
      Simple=.true.
      read(10,*,err=100,end=100)r,sc,rs
      rmax=r
      scl=sc
      rscale=rs
      read(10,*,err=100,end=100)cut,diags,imode
      diagscale=diags
      cutoff=cut
      if(imode .ne. 0)Simple=.false.
100   close(unit=1)
      close(unit=10)
!     General framework
!
!     Initialize
      call init(pred,pos,br1,br2,N1,N2)
!
!     Some user stuff
      write(6,1000)rmax,scl,cutoff,rscale,diagscale
1000  format(' Pair-wise Hessian Approximation'/ &
      ' Version 1.0, February 2006'// &
      ' Spring Model'/ &
      ' Maximum neighbor distance ',F10.5,' a.u.'/ &
      ' Decay ',F10.5,' and cutoff ',F10.5/ &
      ' Hessian rescale ',F10.4,' Diagonal additive term ',F10.4)
        if(simple)then
                write(6,*)'Harmonic in dx**2+dy**2+dz**2 mode'
        else
                write(6,*)'Harmonic in (r-rmin)**2 mode'
        endif
      write(6,*)
      write(6,*)'Bravais Lattice to cartesian co-ordinates'
      write(6,1002)br1
1002  format('   ',3D16.8)
      write(6,*)' '

!
!     Find symmetry operations that relate reduced to larger set
      call findsymm(pred,pos,N1,N2,br1)
!
!     Find equivalents
      call findequivs(pred)
!     Test
      call daverage(pred)
      isfree(1:3,1:nat)=.false.
      depends(1:3,1:nat)=.false.
      diagsval(1:3,1:nat)='  Fixed '
      dmatrix(1:3,1:3,1:nat)=0.D0
      dmat2(1:3,1:3,1:nat)=0.D0
      icanopt=0
      do n=1,nat
         do i=1,3
            ptmp(1:3,1:nat)=pred(1:3,1:nat)
            disp=1.D-4/alat(i)
            ptmp(i,n)=ptmp(i,n)+disp
            call daverage(ptmp)
            do k=1,3
                t1=(ptmp(k,n)-pred(k,n))/disp
                if(abs(t1) .gt. 1D-8)then
                           dmatrix(k,i,n)=t1
                           dmat2(k,i,n)=1.D0
                endif
                t=(ptmp(k,n)-pred(k,n))*alat(k)
                if(abs(t) .gt. 1D-8)then
                          isfree(k,n)=.true.
                          diagsval(k,n)='  Free  '
                          icanopt=1
                endif
            enddo
         enddo
         do i=1,3
                t=0
                do k=1,3
                t=t+abs(dmatrix(k,i,n))
                enddo
                if(abs(t) .lt. 1D-8)depends(i,n)=.true.
         enddo
         write(6,60)n,((dmatrix(k,i,n),k=1,3),depends(i,n),isfree(i,n),i=1,3)
60       format('Connectivity for atom ',i3,/,(3D15.6,' Dependent ',L2,' Free ',L2))
      enddo
!     Apparently no free variables
      if(icanopt .eq. 0)goto 980
!
!     Steps to use, 1d-5 a.u.
      do i=1,3
!        1D-6 Angstroms
         steps(i)=1d-6/alat(i)
      enddo

!     Find Neighbours
      call findneigh(translat,nuse,distances,pred,br1,N1,N2,rmax)
!
!
!     Core, generate the Hessian
      call makehessb(translat,nuse,distances,pred,hess,br1,N1,N2,steps)

!
!     Fixup symmetry issues
	j=0
	do j1=1,nat
	  do i1=1,3
		j=j+1
		k=0
		do j2=1,nat
		   do i2=1,3
			k=k+1
			  htmp(i2,i1,j2,j1)=hess(k,j)
		   enddo
	        enddo
	  enddo
	enddo

	do j1=1,nat
	   do j2=1,nat
		t=0.D0
		do i1=1,3
		   do i2=1,3
			t=0
			do i3=1,3
			do i4=1,3
			t=t+htmp(i3,i4,j2,j1)*dmatrix(i2,i3,j2)*dmatrix(i1,i4,j1)
			enddo
			enddo
			tmp(i2,i1)=t
		   enddo
		enddo
		do i1=1,3
		   do i2=1,3
			htmp(i2,i1,j2,j1)=tmp(i2,i1)
		   enddo
		enddo
	   enddo
	enddo
	j=0
	do j1=1,nat
	  do i1=1,3
		j=j+1
		k=0
		do j2=1,nat
		   do i2=1,3
			k=k+1
			hess(k,j)=htmp(i2,i1,j2,j1)
		   enddo
	        enddo
	  enddo
	enddo
!
      do j=1,n1*3
         do k=1,n1*3
            t=(hess(k,j)+hess(j,k))*0.5D0
            hess(k,j)=t
            hess(j,k)=t
         enddo
      enddo

!     Trace, and trap zero diagonals
      trace=0.D0
      ntrace=0
      do j=1,N1*3
         if(hess(j,j) .gt. 1d-5)then
                trace=trace+hess(j,j)
                ntrace=ntrace+1
         endif
      enddo
!
      trace=trace/dble(ntrace)
      rscale=2.D0*rscale/trace
!     Rescale, and trap zeros on diagonal
      write(6,*)' '
      do j=1,N1*3
          do k=1,n1*3
           hess(k,j)=hess(k,j)*rscale
           if(abs(hess(k,j)) .lt. 1d-8)hess(k,j)=0.D0
           hold(k,j)=hess(k,j)
          enddo
          hess(j,j)=hess(j,j)+diagscale
          hold(j,j)=hess(j,j)
      enddo

      write(21,2101)
2101  format('PORT 2.0 0.35      # PORT/NEWT;  tolf, Initial Trust Radius')
      write(6,*)'Diagonal Elements of Approximate Hessian'
      write(6,*)'****************************************'
      do j=1,N1
         write(6,11)j,(hess(i,i),i=3*j-2,3*j),(diagsval(i,j),i=1,3)
11      format(' Hxx, Hyy, Hzz atom ',i3,3F10.5,' | ',3a) 
      X1=0
      X2=0
      X3=0
      if(isfree(1,j))X1=1.0
      if(isfree(2,j))X2=1.0
      if(isfree(3,j))X3=1.0
      write(21,2102)X1,X2,X3,1.0,'   #Atom ',j,' Generated by pairhess'
2102  format(f3.1,3F4.1,a,i2,a)
      enddo

!     Copy over free elements only
      if1=0
      i1=0
      do j1a=1,n1
        do j1b=1,3
           i1=i1+1
           if(isfree(j1b,j1a))then
            if1=if1+1
            i2=0
            if2=0
            do j2a=1,n1
              do j2b=1,3
                i2=i2+1  
                if(isfree(j2b,j2a))then
                        if2=if2+1
                        hred(if2,if1)=hess(i2,i1)
                        heigen(if2,if1)=-hess(i2,i1)
                endif
              enddo
            enddo
           endif
        enddo
     enddo
     do j=1,if1
        heigen(j,j)=heigen(j,j)+1
     enddo
!     Eigenvalues of reduced matrix
!      dimension W(3*natmax),WORK(10*natmax)
      LWORK=10*natmax
      call DSYEV('N', 'U', if1, heigen, 3*natmax, W1, WORK, LWORK, INFO )
      if(info .ne. 0)then
        write(*,*)'Death in dsyev, info=',info
        stop 'Chaos'
      endif
      do j=1,if1
        do k=1,if1
          heigen(j,k)=hred(j,k)
        enddo
      enddo
      call DSYEV('N', 'U', if1, heigen, 3*natmax, W2, WORK, LWORK, INFO )
      if(info .ne. 0)then
        write(*,*)'Death in dsyev, info=',info
        stop 'Chaos'
      endif
      write(6,*)'Eigenvalues'
      do j=1,if1
        write(6,703)j,W1(j),W2(j)
      enddo
703   format(i4,2F12.6)
      write(*,*)'Max min of Eigenvalues ',W2(1),w2(if1)
      write(*,*)'Max min of I-H Eigen   ',W1(1),w1(if1)
!     Calculate condition number for free elements only
      call matcno(hred,3*natmax,if1,condfree,iflag)
      write(*,*)'Free elements condition number ',condfree
      write(6,*)'Free elements condition number ',condfree
!     Calculate condition number
      call matcno(hold,3*natmax,3*N1,CONDNO,IFLAG)
        if( iflag .eq. -2) then
         write(*,*)'*** Code Error: Increase MX in matinv of matcon.f'
         write(6,*)'*** Code Error: Increase MX in matinv of matcon.f'
        else if(iflag .eq. -1)then
         write(*,*) '*** Code Error: N>LDA, and chaos may ensue'
         write(6,*) '*** Code Error: N>LDA, and chaos may ensue'
        else if(iflag .gt. 0)then
         write(*,*) '*** Warning: Hessian matrix appears to be singular'
         write(6,*) '*** Warning: Hessian matrix appears to be singular'
        else
         write(*,606) condno,ntrace,condno/sqrt(dble(ntrace))
         write(6,606) condno,ntrace,condno/sqrt(dble(ntrace))
606      format(' Condition number ',D12.5,' for ',i3,' variables; Cond/sqrt(N)=',D12.5)
        endif
!
!     Do Cholesky factorization without pivoting
      do j=1,N1*3
         write(22,10)(hess(k,j),k=1,j)
!         write(6,10)(hess(k,j),k=1,N1*3)
      enddo
      job=0
      call dchdc(hess,3*natmax,3*N1,work,jpvt,job,info)
      if(info .lt. 3*n1)goto 970
!
!     Output and make coffee
      open(unit=10,file='.minpair')
      do j=1,N1*3
         write(10,10)(hess(k,j),k=1,j)
!         write(6,10)(hess(k,j),k=1,N1*3)
      enddo
10    format(8D16.8)
!
!     Done
        CALL errclr(ename)
        inquire(unit=6,name=fname)
        i=len_trim(fname)
        write(*,*)'Check .minpair, the estimate, and output in ',fname(1:i)
        STOP ' PairHess END'
!
!        error handling
!
910 INFO = 1
!
!        pairhess.def couldnt be opened
!
        WRITE (ERRMSG,9000) FNAME
        CALL OUTERR('PairHess',ERRMSG)
        GOTO 999
920     INFO = 2
!
!        file FNAME couldn''t be opened
!
        WRITE (ERRMSG,9010) IUNIT
        CALL OUTERR('PairHess',ERRMSG)
        WRITE (ERRMSG,9020) FNAME
        CALL OUTERR('PairHess',ERRMSG)
        WRITE (ERRMSG,9030) STATUS, FORM
        CALL OUTERR('PairHess',ERRMSG)
        GOTO 999
960 INFO = 7
!        Error reading file 
        WRITE (ERRMSG,9050) FNAME
        CALL OUTERR('PairHess',ERRMSG)
        GOTO 999

970     Continue
        write(ERRMSG,9040)
        call outerr('PairHess',ERRMSG)
        write(*,*)'Error, info= ',info,' in Cholesky factorization'
        info=info+1
        write(*,*)'Error, column ',info,' of Hessian not definite'
        i1=(info-1)/3
        i2=info-i1*3
        if(i2.eq.1)then
                write(6,9041)'X',i1+1
                write(*,9041)'X',i1+1
        else if(i2.eq.2)then
                write(6,9041)'Y',i1+1
                write(*,9041)'Y',i1+1
        else if(i2.eq.3)then
                write(6,9041)'Z',i1+1
                write(*,9041)'Z',i1+1
        endif
9041    format('Problem with ',a,' coordinate, atom ',i3)
        write(6,9042)(hold(j,info),j=1,n1*3)
9042    format(8D14.5)
        GOTO 999
980     write(ERRMSG,9043)
9043    format('Error, structure does not appear to have variable positions!')
        call outerr('PairHess',ERRMSG)
        write(6,9043)
        goto 999

999     STOP 'PairHess - Error. Check file pairhess.error.'

9000    FORMAT('can''t open definition file ',A40)
9010    FORMAT('can''t open unit: ',I2)
9020    FORMAT('       filename: ',A50)
9030    FORMAT('         status: ',A,'  form: ',A)
9040    FORMAT('Error reading file: ',A47)
9050    FORMAT('Error, Hessian not positive definite')
      end
