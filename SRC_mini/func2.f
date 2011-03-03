      real*8 function func2(x)
      use mxpmgrid
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      dimension x(*)
!      common    /fmin/alat(3),xwert(mxpm,maxit),ewert(maxit), &
!                      fwert(mxpm,maxit),nparam,iter,ifnc,igrdc 
!      dimension xwert(mxpm,maxit),ewert(maxit),fwert(mxpm,maxit)
      common /ctoler/ tolx
      save ipos
      data eps/1.d-4/
!     Atomic positions have 1d-9 precision in case.struct
!     They are held, in angstroms in the code
!     A conservative estimate of the available precision is thus
!     1d-8*max(a,b,c)
      eps=1d-8*max(alat(1),alat(2),alat(3))
!
      write(6,*) 'in func,nparam',nparam
      do 10 ii=1,nparam
         ial = mod(ii-1,3) + 1
!     Note: this is dangerous, could create problems in drmng
!     It might be better to use a temporary array
         x(ii) = dmod( x(ii) + alat(ial), alat(ial) )
 10   CONTINUE
 1001 format(3f20.12 )
      ipos=iter+1
      do 20 ii=1,iter
        do 30 ip=1,nparam
          if ((abs(x(ip)-xwert(ip,ii))).gt.eps) goto 20
   30   continue
        ipos=ii
!       If we are here, we found it
        goto 400
   20 continue
!     Not available
900   do 40 ip=1,nparam
        xwert(ip,ipos)=x(ip)
   40 continue
!
      write(6,2020)
2020  format(15x,'(func)  call outp')
      func2=999.111d0
!      call outp
       return
      
400     func2=ewert(ipos)
        ifnc=ifnc+1
!       write(6,2000) ifnc,ipos
2000    format(15x,'(func)  ifnc:',i4,'  ipos:',i4 )
        return
!
      end
      subroutine dfunc2(x,g)
      use mxpmgrid
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      dimension x(*),g(*)
!      common    /fmin/alat(3),xwert(mxpm,maxit),ewert(maxit), &
!                      fwert(mxpm,maxit),nparam,iter,ifnc,igrdc 
!      dimension xwert(mxpm,maxit),ewert(maxit),fwert(mxpm,maxit)
      common /ctoler/ tolx
      save ipos
      data eps/1.d-4/
!     Atomic positions have 1d-9 precision in case.struct
!     They are held, in angstroms in the code
!     A conservative estimate of the available precision is thus
!     1d-8*max(a,b,c)
      eps=1d-8*max(alat(1),alat(2),alat(3))
!
      do 50 ii=1,nparam
         ial = mod(ii-1,3) + 1
!     Note: this is dangerous, could create problems in drmng
!     It might be better to use a temporary array
         x(ii) = dmod( x(ii) + alat(ial), alat(ial) )
   50 continue
!
      do 60 ii=1,iter
        do 70 ip=1,nparam
          if ((abs(x(ip)-xwert(ip,ii))).gt.eps) goto 60
   70   continue
        ipos=ii
   60 continue
!
      do 80 ip=1,nparam
        if ((abs(x(ip)-xwert(ip,ipos))).gt.eps) goto 910
        g(ip)=fwert(ip,ipos)
   80 continue
      igrdc=igrdc+1
!     write(6,2030) igrdc,ipos
2030  format(15x,'(dfunc)  igrdc:', i4,'  ipos:',i4)
      return
!        error handling
!
  910 INFO = 2
      CALL OUTERR('FUNC','  stop_fcn_grn_call' )
      write(6,5000)(x(ii),ii=1,nparam)
5000  format(' >>> (grd) x for grd- and fcn-call not equal !!!'/  &
       '    x = '/,100(3f20.10/) )
      CALL OUTERR('FUNC','>>> (grd) x for grd- and fcn-call not equal')
      GOTO 999
  999 STOP 'MINI - Error'
      end
