!
      subroutine dfpmin(p,n,iter,fret,nat)
      implicit double precision (a-h,o-z)
      include 'param.inc'
      integer iter,n,nmax,itmax
      real*8 fret,p(n),func,eps,stpmx,tolx
      parameter (itmax=200,stpmx=100.,eps=3.d-8)
      external dfunc,func
      integer i,its,j,nat
      logical check
      CHARACTER*67       ERRMSG
      real*8 den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,test &
        ,stp 
       real*8,allocatable :: dg(:),g(:),hdg(:),hessin(:,:) &
        ,pnew(:),xi(:)
      common /ctoler/tolx
         nmax=nat*3
       allocate ( dg(nmax),g(nmax),hdg(nmax),hessin(nmax,nmax) &
        ,pnew(nmax),xi(nmax))
      write(6,2010)'  call func'
      fp=func(p)
      if (fp.eq.999.111d0) return
      write(6,2010)'  call dfunc'
      call dfunc(p,g)
      sum=0.d0
      do 10 i=1,n
        do 20 j=1,n
          hessin(i,j)=0.d0
 20     CONTINUE 
        hessin(i,i)=1.d0
        xi(i)=-g(i)
        sum=sum+p(i)**2
 10   CONTINUE 
!     stpmax=stpmx*dmax1(dsqrt(sum),dfloat(n))
      do 30 its=1,itmax
        iter=its
        write(6,2010)'  call maxstp ' 
2010    format(10x,'(dfpmin)',a) 
        call maxstp(p,xi,n,stpmax,nat)
        write(6,2011)stpmax
2011    format(10x,'(dfpmin)  stpmax = ',e20.10)
!
!.......check max. step due to geometry 
        stp = 0.d0
        do 40 i=1,n
           stp = dmax1(stp,dabs(xi(i)))
 40     CONTINUE
        write(6,'(10x,19h(dfpmin)  stp    = ,e20.10)')stp
        stp = stp * stpmax
        write(6,'(10x,19h(dfpmin)  stpmax = ,e20.10)')stp
        if( stp .le. tolx ) goto 900
!.......end check
!
        write(6,2010) '   call lnsrch'
        stpmax = dmin1(1.d0,stpmax)

        call lnsrch(n,p,   fp,  g,xi,pnew,fret,stpmax,check)
        if (fret.eq.999.111d0) return
!            lnsrch(n,xold,fold,g,p, x,   f,   stpmax,check)
        if(check)then
           write(*,*)'  !!!!!    check = .true.   !!!!! '
           write(*,*)'             ATTENTION  '
        endif
        write(6,2010)'  pnew'
        write(6,550)(pnew(i),i=1,n)
 550    format( 20x,3f20.10 )
        fp=fret
        do 50 i=1,n
          xi(i)=pnew(i)-p(i)
          p(i)=pnew(i)
 50     CONTINUE
        test=0.d0
        do 60 i=1,n
          test = test + xi(i)**2
 60     CONTINUE
        test = dsqrt(test)
        write(6,2010)'  test of xyz convergence : '
        write(6,560)test
 560    format(10x,'(dfpmin)      test =',f20.10)
        write(6,570)tolx
 570    format(10x,'(dfpmin)      tolx =',f20.10 )
        if(test.lt.tolx )return
        do 70 i=1,n
          dg(i)=g(i)
 70     CONTINUE
        write(6,2010)'  call dfunc'
        call dfunc(p,g)
        den=dmax1(fret,1.d0)
        do 80 i=1,n
          dg(i)=g(i)-dg(i)
 80     CONTINUE
        do 90 i=1,n
          hdg(i)=0.
          do 90 j=1,n
            hdg(i)=hdg(i)+hessin(i,j)*dg(j)
 90     CONTINUE
        fac=0.d0
        fae=0.d0
        sumdg=0.d0
        sumxi=0.d0
        do 100 i=1,n
          fac=fac+dg(i)*xi(i)
          fae=fae+dg(i)*hdg(i)
          sumdg=sumdg+dg(i)**2
          sumxi=sumxi+xi(i)**2
 100    CONTINUE
        if(fac**2.gt.eps*sumdg*sumxi) then
          fac=1.d0/fac
          fad=1.d0/fae
          do 110 i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
 110      CONTINUE
          do 120 i=1,n
            do 120 j=1,n
              hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j) &
                -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
 120      CONTINUE
        endif
!        write(6,*) 'hessin'
!        do i=1,n
!          write(6,'(20f6.2)') (hessin(i,j),j=1,n)
!        enddo
        do 130 i=1,n
          xi(i)=0.d0
          do 130 j=1,n
            xi(i)=xi(i)-hessin(i,j)*g(j)
 130    CONTINUE
!     calculate the skalar product between 
!     new search direction and gradient
        xig=0.d0
        do 140 i=1,n
          xig=xig+xi(i)*g(i)
 140    CONTINUE
        write(6,'(10x,''(dfpmin)  xi*g='',f20.10)') xig
        if(xig.ge.0.d0) then
          write(6,'(10x,''(dfpmin)  xi*g>0 -> initialize Hessian'')')
          do 150 i=1,n
            do 160 j=1,n
              hessin(i,j)=0.d0
 160        CONTINUE
            hessin(i,i)=1.d0
            xi(i)=-g(i)
 150      CONTINUE
        endif
 30   CONTINUE
      GOTO 910
!
!        error handling
!
  900 INFO = 1
!
!        muffin-tin radii too large
!
      CALL OUTERR('DFPMIN','  stop  after  maxstp  in mini-run ')
      CALL OUTERR('DFPMIN','  stp.le.tolx ')
      WRITE (ERRMSG,'(a12,f15.6)') '  tolx = ',tolx
      CALL OUTERR('DFPMIN',ERRMSG)
      WRITE (ERRMSG,'(a12,f15.6)') '  stp  = ',stp
      CALL OUTERR('DFPMIN',ERRMSG)
      CALL OUTERR('DFPMIN', &
                  '  possible reason: muffin-tin radii too large')
      WRITE (ERRMSG,'(a12,f15.6)') '  stpmax = ',stpmax
      CALL OUTERR('DFPMIN',ERRMSG)
      GOTO 999
  910 INFO = 2
!
!     
!
      CALL OUTERR('DFPMIN','  too many iteration ')
      GOTO 999
  999 STOP 'MINI - Error'
      end
