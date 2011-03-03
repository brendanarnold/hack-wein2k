      program eosfit
      implicit real*8 (a-h,o-z)
      CHARACTER*11      STATUS,FORM                                     
      CHARACTER*79      TITLE,NAME                                            
      CHARACTER*80      DEFFN, ERRFN ,fname                                   
      common voleos2(400), eeos2(400), volmurna(400), emurna(400), volbm(400), ebm(400)
      dimension  x(4),y(400),diag(100),fjac(100,4),yeos2(100)
      dimension ipvt(4),qtf(4),wa1(4),wa2(4),wa3(4),wa4(100)
      dimension vol(400), e(400) 
      dimension ymurna(400), ybm(400)
      external eos2,murna,birchmurna
!
      CALL GTFNAM(DEFFN,ERRFN)
      OPEN(1,FILE=DEFFN,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING eosfit.def !!!!'
      STOP 'eosfit.def'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'                                
 8001 CONTINUE
!
      ftol=1.d-16
      i=1
 1    read(55,*,end=2) vol(i),e(i)
!      vol(i)=vol(i)**3/4.d0
	voleos2(i)=vol(i)
	eeos2(i)=e(i)
	volmurna(i)=vol(i)
	emurna(i)=e(i)
	volbm(i)=vol(i)
	ebm(i)=e(i)
      i=i+1
      goto 1
 2    continue
      x(1)=e(1)
      x(2)=1.d0
      x(3)=-10.d0
      x(4)=100.d0
      m=i-1
      n=4
      maxfev=100000
      epsfcn=0.01d0
      mode=1
      factor=100.d0
      nprint=0
      ldfjac=100
      call lmdif(eos2,m,n,x,y,ftol,ftol,ftol,maxfev,epsfcn,diag,mode, &
      factor,nprint, &
       info,nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

      write(66,*) 'Equation of state: EOS2 (PRB52,8064)        info',info
      write(66,100) x
 100  format(' a,b,c,d',4f17.6)
      call eos2(i-1,4,x,yeos2,iflag)
!      do j=1,i-1
!      write(*,*) vol(j),eeos2(j),yeos2(j)
!      enddo
      v01=-x(3)/x(2)+sqrt((x(3)/x(2))**2-3.d0*x(4)/x(2))
      v02=-x(3)/x(2)-sqrt((x(3)/x(2))**2-3.d0*x(4)/x(2))
!      write(*,*) 'v0',v01**3,v02**3
      B01=-2.d0/9.d0*v01**(-4)*(x(2)+x(3)*v01**(-1))*14703.6
      B02=-2.d0/9.d0*v02**(-4)*(x(2)+x(3)*v02**(-1))*14703.6
      b0p1=(10.d0/3.d0*x(2)+11.d0/3.d0*x(3)*v01**(-1))/ &
           (x(2)+x(3)*v01**(-1))
      b0p2=(10.d0/3.d0*x(2)+11.d0/3.d0*x(3)*v02**(-1))/ &
           (x(2)+x(3)*v02**(-1))
      bdf=x(2)
      cdf=x(3)
!      write(*,*) 'B0',b01,b02
      if(abs(v01**3-vol(3)).lt.abs(v02**3-vol(3))) then
        write(66,101) v01**3,b01,b0p1
      else
        write(66,101) v02**3,b02,b0p2
      endif
!

!...murnaghan
      x(1)=emurna(3)
      x(2)=vol(3)
      x(3)=100.d0
      x(4)=5.d0
       n=3
      call lmdif(murna,m,n,x,y,ftol,ftol,ftol,maxfev,epsfcn,diag,mode, &
      factor,nprint, &
       info,nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!      write(*,*) 'Equation of state: Murnaghan (3)            info',info
!      write(*,101) (x(j),j=2,4),x(1)
       n=4
      call lmdif(murna,m,n,x,y,ftol,ftol,ftol,maxfev,epsfcn,diag,mode, &
      factor,nprint, &
       info,nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

      write(66,*)
      write(66,*) 'Equation of state: Murnaghan                info',info
      write(66,*) 'E=E0+[B*V/BP*(1/(BP-1)*(V0/V)**BP +1)-B*V0/(BP-1)]/14703.6'
      write(66,*) 'Pressure=B/BP*((V0/V)**BP -1)'
!      write(9,101) (x(j),j=2,4),x(1)
      write(66,101) (x(j),j=2,4),x(1)
 101  format(' V0,B(GPa),BP,E0',3f15.4,f18.6)
      call murna(i-1,4,x,ymurna,iflag)
      write(66,102)
 102  format(9x,'vol',7x,'energy',9x,'de(EOS2)',6x,'de(Murnaghan)  Pressure(GPa)')
      sigma1=0.d0
      sigma2=0.d0
      xmin=10000.
      xmax=-10000.
      do j=1,i-1
      pressure=x(3)/x(4)*( (x(2)/volmurna(j))**x(4) -1 )
      write(66,103) volmurna(j),emurna(j),yeos2(j),ymurna(j),pressure
      sigma1=sigma1+yeos2(j)**2
      sigma2=sigma2+ymurna(j)**2
 103  format(f13.4,f16.6,2f13.6,f15.3)
      if(xmax.lt.volmurna(j)) xmax=volmurna(j)
      if(xmin.gt.volmurna(j)) xmin=volmurna(j)
      enddo
      write(66,104) sqrt(sigma1/(i-1)),sqrt(sigma2/(i-1))
 104  format('                  Sigma:',4x,f14.6,f13.6)
!
!.....create plotfile
      xmin=xmin-(xmax-xmin)*0.1
      xmax=xmax+(xmax-xmin)*0.1
      dx=(xmax-xmin)/399
      do j=1,400
       emurna(j)=0.d0
      volmurna(j)=xmin+(j-1)*dx
      enddo
      call murna(400,4,x,ymurna,iflag)
!      write(9,*)
!      write(9,*)'Murnaghan equation of state'
!      write(9,*)'E=E0+[B*V/BP*(1/(BP-1)*(V0/V)**BP +1)-B*V0/(BP-1)]/14703.6'
!      write(9,*)'Pressure=B/BP*((V0/V)**BP -1)'
!      write(9,*)
!      write(9,*)'fitted equilibrium volume V0 = ',x(2),' au**3'
!      write(9,*)'fitted bulk modulus B = ',x(3),' GPa'
!      write(9,*)'fitted first derivative of B w.r.t. pressure (BP) = ',x(4), &
!	'(dimensionless)'
!      write(9,*)'fitted lowest energy E0 = ',x(1), ' Ry'
!      write(9,*)
!      write(9,125)
 125  format(3x,'Volume',10x,'Energy',15x,'Pressure')
      do j=1,400
      pressure=x(3)/x(4)*( (x(2)/volmurna(j))**x(4) -1 )
      write(9,126) volmurna(j),ymurna(j),pressure
 126  format(f12.6,4x,f16.6,5x,f10.4)
      enddo
!

!...birch-murnaghan
      x(1)=ebm(3)
      x(2)=vol(3)
      x(3)=100.d0
      x(4)=5.d0
       n=3
      call lmdif(birchmurna,m,n,x,y,ftol,ftol,ftol,maxfev,epsfcn,diag,mode, &
      factor,nprint, &
       info,nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!      write(*,*) 'Equation of state: Birch-Murnaghan (3)            info',info
!      write(*,101) (x(j),j=2,4),x(1)
       n=4
      call lmdif(birchmurna,m,n,x,y,ftol,ftol,ftol,maxfev,epsfcn,diag,mode, &
      factor,nprint, &
       info,nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

      write(66,*)
      write(66,*) 'Equation of state: Birch-Murnaghan                info',info
      write(66,*) 'E = E0 + 9/16*(B/14703.6)*V0*[(eta**2-1)**3*BP + (eta**2-1)**2*(6-4*eta**2)]'
      write(66,*) ' 	  --> eta = (V0/V)**(1/3)'
      write(66,*) 'Pressure = 3/2*B*(eta**7 - eta**5)*(1 + 3/4*(BP-4)*[eta**2 - 1])'
!      write(9,101) (x(j),j=2,4),x(1)
      write(66,101) (x(j),j=2,4),x(1)
      call birchmurna(i-1,4,x,ybm,iflag)
      write(66,112)
 112  format(9x,'vol',7x,'energy',9x,'de(Birch-Murnaghan)  Pressure(GPa)')
      sigma1=0.d0
      sigma2=0.d0
      xmin=10000.
      xmax=-10000.
      do j=1,i-1
      eta=(x(2)/volbm(j))**(1.d0/3.d0)
      pressure=1.5d0*x(3)*(eta**7.d0 - eta**5.d0) * &
 	(1.d0 + .75d0*(x(4)-4.d0)*(eta**2-1.d0))
      write(66,113) volbm(j),ebm(j),ybm(j),pressure
      sigma2=sigma2+ybm(j)**2
 113  format(f13.4,f16.6,f13.6,6x,f15.3)
      if(xmax.lt.volbm(j)) xmax=volbm(j)
      if(xmin.gt.volbm(j)) xmin=volbm(j)
      enddo
      write(66,114) sqrt(sigma2/(i-1))
 114  format('                  Sigma:',4x,f14.6)
      write(66,*)
!
!.....create plotfile
      xmin=xmin-(xmax-xmin)*0.1
      xmax=xmax+(xmax-xmin)*0.1
      dx=(xmax-xmin)/399
      do j=1,400
       ebm(j)=0.d0
      volbm(j)=xmin+(j-1)*dx
      enddo
      call birchmurna(400,4,x,ybm,iflag)
!      write(11,*)
!      write(11,*)'Birch-Murnaghan equation of state'
!      write(11,*)'E = E0 + 9/16*(B/14703.6)*V0*[(eta**2-1)**3*BP + (eta**2-1)**2*(6-4*eta**2)]'
!      write(11,*)'       --> eta = (V0/V)**(1/3)'
!      write(11,*)'Pressure = 3/2*B*(eta**7 - eta**5)*(1 + 3/4*(BP-4)*[eta**2 - 1])'
!      write(11,*)
!      write(11,*)'fitted equilibrium volume V0 = ',x(2),' au**3'
!      write(11,*)'fitted bulk modulus B = ',x(3),' GPa'
!      write(11,*)'fitted first derivative of B w.r.t. pressure (BP) = ',x(4), &
!	'(dimensionless)'
!      write(11,*)'fitted lowest energy E0 = ',x(1), ' Ry'
!      write(11,*)
!      write(11,125)
      do j=1,400
      eta=(x(2)/volbm(j))**(1.d0/3.d0)
      pressure=1.5d0*x(3)*(eta**7.d0 - eta**5.d0) * &
        (1.d0 + .75d0*(x(4)-4.d0)*(eta**2-1.d0))
      write(11,126) volbm(j),ybm(j),pressure
      enddo
!


!.....estimate of B and BP using exp. Volume:
5      write(*,*)
       write(*,*) 'Enter V0-exp. for an estimate of B and BP using EOS2 ' &
        ,'and the experimental'
       write(*,*) 'equillibrium volume. Make sure there is a `.` in the ' &
	,'number. (Fill out `0` if'
       write(*,*) 'you are interested the (Birch-)Murnaghan equation of ' &
	,'state only.) :'
      read(*,'(f20.5)') v01
      if(v01.eq.0.d0) stop
      B01=-2.d0/9.d0*v01**(-4.d0/3.d0)*(bdf+cdf*v01**(-1.d0/3.d0)) &
          *14703.6
      b0p1=(10.d0/3.d0*bdf+11.d0/3.d0*cdf*v01**(-1.d0/3.d0))/ &
           (bdf+cdf*v01**(-1.d0/3.d0))
        write(66,*) 'Estimate of B and BP using EOS2 and V0-exp:'
        write(66,105) v01,b01,b0p1
 105    format('V-exp,B(GPa),BP',3f15.4,f18.6)

      stop
      end
!
      subroutine eos2(m,n,x,fvec,iflag)
      implicit real*8 (a-h,o-z)
      common voleos2(400), eeos2(400), volmurna(400), emurna(400), volbm(400), ebm(400)
      dimension fvec(*),x(*)
      fsum=0.d0
      do i=1,m
      fvec(i)=x(1)+x(2)*voleos2(i)**(-1.d0/3.d0)+x(3)*voleos2(i)**(-2.d0/3.d0)+ &
              x(4)*voleos2(i)**(-1.d0)  -  eeos2(i)
      fsum=fsum+fvec(i)**2
!      write(*,*) voleos2(i),eeos2(i),fvec(i)
      enddo
!      write(*,*) m,(x(i),i=1,4),fsum
      return
      end

      subroutine murna(m,n,x,fvec,iflag)
      implicit real*8 (a-h,o-z)
      common voleos2(400), eeos2(400), volmurna(400), emurna(400), volbm(400), ebm(400)
      dimension fvec(*),x(*)
      if (n.eq.3) then
       xh=x(4)
       x(4)=5.0d0
      endif
      fsum=0.d0
      do i=1,m
      fvec(i)=x(1)+x(3)/14703.6d0*volmurna(i)/(x(4)*(x(4)-1.d0))*(x(4)* &
              (1.d0-x(2)/volmurna(i))+(x(2)/volmurna(i))**x(4)-1.d0)  - emurna(i)
!      fvec(i)=x(1)+x(2)*volmurna(i)**(-1.d0/3.d0)+x(3)*volmurna(i)**(-2.d0/3.d0)+
!     *        x(4)*volmurna(i)**(-1.d0)  -  e(i)
      fsum=fsum+fvec(i)**2
!      write(66,*) volmurna(i),emurna(i),fvec(i)
      enddo
!      write(*,*) m,(x(i),i=1,4),fsum
      if (n.eq.3) then
       x(4)=xh
      endif
      return
      end

      subroutine birchmurna(m,n,x,fvec,iflag)
      implicit real*8 (a-h,o-z)
      common voleos2(400), eeos2(400), volmurna(400), emurna(400), volbm(400), ebm(400)
      dimension fvec(*),x(*)
      if (n.eq.3) then
       xh=x(4)
       x(4)=5.0d0
      endif
      fsum=0.d0
      do i=1,m
      eta=(x(2)/volbm(i))**(1.d0/3.d0)
      fvec(i)=x(1)+((9.d0*(x(3)/14703.6d0)*x(2))/16.d0)* &
	( (eta**2.d0-1.d0)**3.d0*x(4) + &
	  (eta**2.d0 - 1.d0)**2.d0*(6.d0-4.d0*eta**2.d0) ) - ebm(i)
      fsum=fsum+fvec(i)**2
!      write(*,*) volbm(i),ebm(i),fvec(i)
      enddo
!      write(*,*) m,(x(i),i=1,4),fsum
      if (n.eq.3) then
       x(4)=xh
      endif
      return
      end

