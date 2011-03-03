      
      subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check)
      implicit double precision (a-h,o-z)
      integer i,n
      logical check
      real*8 f,fold,stpmax,g(n),p(n),x(n),xold(n),func,alf,tolx &
        ,a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope &
        ,sum,test,tmplam
      parameter(alf=1.d-4)
      common /ctoler/tolx
      external func
      
!     write(6,2020)n
!2020  format(15x,'(lnsrch) n:',i5)
      write(6,2030)stpmax
2030  format(15x,'(lnsrch) stpmax:',f20.10)

      check=.false.
      sum=0.d0
      do 10 i=1,n
        sum=sum+p(i)*p(i)
 10   CONTINUE
      sum=dsqrt(sum)
!gv   if(sum.gt.stpmax) then
        do 20 i=1,n
          p(i)=p(i)*stpmax
!gv       p(i)=p(i)*stpmax/sum
 20   CONTINUE
!gv   endif
      slope=0.d0
      do 30 i=1,n
        slope=slope+g(i)*p(i)
 30   CONTINUE
      test=0.d0
      do 40 i=1,n
        test = test + p(i)**2
 40   CONTINUE
      test = dsqrt(test)
      if( test .eq. 0.d0 )then
         alamin = 1.d39
      else
         alamin=tolx/test
      endif
      write(6,'(15x,8h(lnsrch),9h  alamax=,f15.7)')1.d0
      write(6,'(15x,8h(lnsrch),9h  alamin=,f15.7)')alamin
      alam=1.d0

    1 continue
!
!gv-start  
      if(alam.lt.alamin)then
         do 50 i=1,n
            x(i)=xold(i)
 50      CONTINUE
         check=.true.
         return
      endif
!gv-end
!
      do 60 i=1,n
        x(i)=xold(i)+alam*p(i)
 60   CONTINUE
      write(6,2050)'  p=x'
      write(6,'(6(26x,3f16.8/))')(x(i),i=1,n)
      write(6,'(15x,8h(lnsrch),7h  alam=,f15.7)') alam
!     write(6,'(15x,a,2hx=)')'(lnsrch)'
!     write(6,'(3(26x,3f16.8/))')(x(i),i=1,n)
!     write(6,'(15x,8h(lnsrch),5hxold=)')
!     write(6,'(3(26x,3f16.8/))')(xold(i),i=1,n)
      write(6,2050)'  call func'
2050  format(15x,'(lnsrch)',a)
      f=func(x)
      if (f.eq.999.111d0) return
!gv   if(alam.lt.alamin) then
!gv     do i=1,n
!gv       x(i)=xold(i)
!gv     enddo
!gv     check=.true.
!gv     return
!gv   elseif(f.le.fold+alf*alam*slope) then
      if(f.le.fold+alf*alam*slope) then
        return
      else
        if(alam.eq.1.d0) then
          tmplam=-slope/(2.d0*(f-fold-slope))
        else
!ccccccc  fix for ifnc .gt. ipos   error
             write(6,*) 'LNSRCH: Energy not lower, return anyway'
             write(6,*) f,fold,alf,alam,slope
             return
!ccccccc
          rhs1=f-fold-alam*slope
          rhs2=f2-fold2-alam2*slope
          a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
          b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/ &
            (alam-alam2)
!         write(*,*)'rhs1 = ',rhs1
!         write(*,*)'rhs2 = ',rhs2
!         write(*,*)'a = ',a
!         write(*,*)'b = ',b
          if(a.eq.0.d0) then
            tmplam=-slope/(2.d0*b)
          else
            disc=b*b-3.d0*a*slope
            if( disc .lt. 0.d0 ) GOTO 900
            tmplam=(-b+dsqrt(disc))/(3.d0*a)
          write(*,*)' slope = ',slope
          endif
          if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
        endif
      endif
      alam2=alam
      f2=f
      fold2=fold
      alam=dmax1(tmplam,alamin)
!gv   alam=dmax1(tmplam,0.1d0 * alam)
      goto 1
!
!        error handling
!
  900 INFO = 1
!
!        muffin-tin radii too large
!
      CALL OUTERR('LNSRCH','  disc_lt_0 ')
      STOP 'MINI - Error'
      end
