      subroutine nwtmin(minmod,pst,n,iter,fret,nat,damp,tolf1)
      implicit double precision (a-h,o-z)
!     implicit none
      include 'param.inc'
      logical tsdp
      integer iter,n,nmax,itmax,info
      real*8 fret,pst(n),func,eps,stpmx,tolx,damp(*)
      CHARACTER*67       ERRMSG
      character*4 minmod
      parameter (itmax=200,stpmx=100.,eps=3.d-8)
      common /ctoler/tolx
      external dfunc,func
      integer i,its,nat
      real*8 fp,stpmax,test,stp,dmp,fion,ff,fm,d_f 
      real*8,allocatable ::  g(:),gold(:),speed(:),p(:)
      nmax=nat*3
      allocate (  g(nmax),gold(nmax),speed(nmax),p(nmax))

      do 10 i=1,n
        gold(i)=0.d0
!        gold(nmax)=0.d0
        p(i)=pst(i)
        speed(i)=0.d0
 10   CONTINUE

      do 20 its=0,itmax
        iter=its

        write(6,2010)'  call func'
        fp=func(p)
        if (fp.eq.999.111d0) return
        write(6,2010)'  call dfunc'
        call dfunc(p,g)
        tsdp=.false.
           write(6,8384) (gold(i),i=1,n)
 8384      format('GOLD',10f10.5)
        do 30 i=1,n
!     force and gradient in Ry/bohr 
          ff=-g(i)
          fm=-gold(i)
!cc
!CCCC          if(its.gt.1) tsdp=.true.
!cc         
          if (tsdp) then
!     compute 'conjugate' correction
            if (ff*fm .ge. 0.0) then
              d_f = 1.5 * ff**2  / (fm**2+1e-12) * fm
              
!     don't allow for changes greater than 2.0 times the force
              if (abs(d_f) .gt. 2.0 * abs(ff)) d_f = 2.0 * ff

            else
!     decrease force magnitude to push somehow in between old and new 
!     coordinate
              if (abs(ff) .lt. abs(fm)) then
                d_f = 0.5 * ff**2  / fm**2 * fm
              else
!     allow at maximum shift back to old coordinate
                d_f = - ff - abs(ff/(ff - fm)) * fm
              endif
            endif
            fion=ff+d_f
            speed(i)=0.d0
            dmp=0.d0

          else
!     take just extrapolated force for moving if not tsdp
            fion=ff
            if(speed(i)*fion.ge.0.d0) then
              if(abs(fion).le.abs(3.d0*tolf1)) then
              dmp=1.-damp(i)
              else
              dmp=damp(i)
              endif
            else
              dmp=1.-damp(i)
            endif

           endif

          write(6,'(10x,''(nwtmin)  i,damp,speed,fion ='',i3,f6.2,3f10.4)') i,dmp,speed(i),fion
          if(minmod.eq.'NEW1')          dmp=damp(i)
          if(minmod.eq.'NEW1')          dmp=0.d0   
          speed(i)=dmp*speed(i)+fion
          gold(i)=g(i)
 30     CONTINUE     
        

        write(6,2010)'  call maxstp ' 
 2010   format(10x,'(nwtmin)',a) 
        call maxstp(p,speed,n,stpmax,nat)
        write(6,2011)stpmax
 2011   format(10x,'(nwtmin)  stpmax = ',e20.10)
!
!.......check max. step due to geometry 
        stp = 0.d0
        do 40 i=1,n
           stp = dmax1(stp,dabs(speed(i)))
 40     CONTINUE 
        write(6,'(10x,19h(nwtmin)  stp    = ,e20.10)')stp
        stp = stp * stpmax
        write(6,'(10x,19h(nwtmin)  stpmax = ,e20.10)')stp
!!blaha!        if( stp .le. tolx ) GOTO 900
!.......end check
!
        stpmax = dmin1(1.d0,stpmax)
        test=0.d0
        do 50 i=1,n
          speed(i)=stpmax*speed(i)
          test=test+speed(i)**2
 50     CONTINUE
        test = dsqrt(test)
        write(6,2010)'  test of xyz convergence : '
        write(6,560)test
 560    format(10x,'(nwtmin)  test =',f20.10)
        write(6,570)tolx
 570    format(10x,'(nwtmin)  tolx =',f20.10 )
 !!!blaha!       if(test.lt.tolx )return

        write(6,2010) '  move to new position'
        do 60 i=1,n
          p(i)=p(i)+speed(i)
 60     CONTINUE
        write(6,2010)'  new p'
        write(6,550)(p(i),i=1,n)
 550    format( 20x,3f12.4 )

 20   CONTINUE
      GOTO 910
!
!        error handling
!
  900 INFO = 1
!
!        muffin-tin radii too large
!
      CALL OUTERR('NWTMIN','  stop  after  maxstp  in mini-run ')
      CALL OUTERR('NWTMIN','  stp.le.tolx ')
      WRITE (ERRMSG,'(a10,f15.6)') '  tolx = ',tolx
      CALL OUTERR('NWTMIN',ERRMSG)
      WRITE (ERRMSG,'(a10,f15.6)') '  stp  = ',stp
      CALL OUTERR('NWTMIN',ERRMSG)
      CALL OUTERR('NWTMIN', &
                  '  possible reason: muffin-tin radii too large')
      WRITE (ERRMSG,'(a12,f15.6)') '  stpmax = ',stpmax
      CALL OUTERR('NWTMIN',ERRMSG)
      GOTO 999
  910 INFO = 2
!
!     
!
      CALL OUTERR('NWTMIN','  too many iteration ')
      GOTO 999
  999 STOP 'MINI - Error'
      end
