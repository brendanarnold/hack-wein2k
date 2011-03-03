      program eosfit
      implicit real*8 (a-h,o-z)
      CHARACTER*11      STATUS,FORM                                     
      CHARACTER*79      TITLE,NAME                                            
      CHARACTER*80      DEFFN, ERRFN ,fname                                   
      parameter (ndata_points=4000)
      parameter (nfit_param=28)
      parameter (liv=82+nfit_param*4)
      parameter (lv=105+nfit_param*(ndata_points+2*nfit_param+21)+2*ndata_points)

      common alat(ndata_points,6), e(ndata_points),ndimfit
      common /other/ xc(nfit_param)
      dimension  x(nfit_param),y(ndata_points),bounds(2,nfit_param)
      integer iv(liv),nni(nfit_param)
      dimension v(lv),xmin(nfit_param),xmax(nfit_param)

      dimension yeos(ndata_points)
!      dimension ipvt(nfit_param),qtf(nfit_param),wa1(nfit_param),wa2(nfit_param),wa3(nfit_param),wa4(400)
      external eos,eosparabol,eosother
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
      x=0.1d0
      ftol=1.d-15
      i=1
      write(*,*) 'Enter dimension of fit (number of variable lattice parameters, 1-6):'
      read(5,*) ndimfit
      if(ndimfit.eq.1) then
        nfit=3
      else if(ndimfit.eq.2) then  
        nfit= 6
      else if(ndimfit.eq.3) then  
        nfit= 10
      x(6)=00.d0
      x(9)=00.d0
      x(10)=00.d0

      else if(ndimfit.eq.4) then  
        nfit=15 !! 14
      else if(ndimfit.eq.5) then  
        nfit=21
      else if(ndimfit.eq.6) then  
        nfit=28
      else  
        write(*,*) ndimfit,' case not supported, must be 1-6'
        stop 'ndimfit'
      endif  
      write(*,*) ndimfit,' fitcase', nfit,' parameter'
      emin=+1.d0
      emax=-99999999999.d0
! 1    read(10,*,end=2) e(i),(alat(i,j),j=1,ndimfit)
 1    read(10,*,end=2) e(i)
      read(11,*,end=2) (alat(i,j),j=1,ndimfit)
      if(e(i).lt.emin) ind_emin=i
      emin=min(e(i),emin)
      emax=max(e(i),emax)
      i=i+1
      goto 1
 2    continue
      print*, 'lowest data point:', e(ind_emin),alat(ind_emin,1:ndimfit)
      bounds(1,:)=-d1mach(2)
      bounds(2,:)=d1mach(2)
      x(1)=emin
      bounds(1,1)=emin+(emin-emax)
      bounds(2,1)=emax
      x(3)=alat(ind_emin,1)
      bounds(1,3)=alat(ind_emin,1)-alat(ind_emin,1)/25.d0    ! bounds +- 4 %
      bounds(2,3)=alat(ind_emin,1)+alat(ind_emin,1)/25.d0
    if(ndimfit.ge.2) then  
      x(5)=alat(ind_emin,2)
      bounds(1,5)=alat(ind_emin,2)-alat(ind_emin,2)/25.d0
      bounds(2,5)=alat(ind_emin,2)+alat(ind_emin,2)/25.d0
    endif
    if(ndimfit.ge.3) then  
      x(8)=alat(ind_emin,3)
      bounds(1,8)=alat(ind_emin,3)-alat(ind_emin,3)/25.d0
      bounds(2,8)=alat(ind_emin,3)+alat(ind_emin,3)/25.d0
    endif
    if(ndimfit.ge.4) then  
      x(12)=alat(ind_emin,4)
      bounds(1,12)=alat(ind_emin,4)-alat(ind_emin,4)/25.d0
      bounds(2,12)=alat(ind_emin,4)+alat(ind_emin,4)/25.d0
    endif
    if(ndimfit.ge.5) then  
      x(17)=alat(ind_emin,5)
      bounds(1,17)=alat(ind_emin,5)-alat(ind_emin,5)/25.d0
      bounds(2,17)=alat(ind_emin,5)+alat(ind_emin,5)/25.d0
    endif
    if(ndimfit.ge.6) then  
      x(23)=alat(ind_emin,6)
      bounds(1,23)=alat(ind_emin,6)-alat(ind_emin,6)/25.d0
      bounds(2,23)=alat(ind_emin,6)+alat(ind_emin,6)/25.d0
     endif
      m=i-1
      iv=0
      v=0.d0
  goto 1234
!purely parabolic fit, no crossterms, 3dim
      n=7
      x(7)=alat(ind_emin,3)
      bounds(1,7)=alat(ind_emin,3)-alat(ind_emin,3)/25.d0
      bounds(2,7)=alat(ind_emin,3)+alat(ind_emin,3)/25.d0
      call dn2fb (m,n,x,bounds,eosparabol,iv,liv,lv,v,ui,ur,uf)
       nf=1
      call eosparabol(m,n,x,nf,yeos,ui,ur,uf)
      write(*,100) (x(i),i=1,nfit)
      write(*,102)
      sigma1=0.d0
      sigma2=0.d0
      do j=1,m
      write(*,103) (alat(j,j1),j1=1,ndimfit),e(j),yeos(j)
      sigma1=sigma1+yeos(j)**2
      enddo
      write(*,104) sqrt(sigma1/(m-1))

      iv=0
      v=0.d0
!purely other fit,  crossterms, 3dim
      n=3
      xc=x
      x(1:3)=0.d0             
      bounds(1,1:3)=-d1mach(2)
      bounds(2,1:3)=d1mach(2)
      call dn2fb (m,n,x,bounds,eosother,iv,liv,lv,v,ui,ur,uf)
       nf=1
      call eosother(m,n,x,nf,yeos,ui,ur,uf)
      write(*,100) (x(i),i=1,nfit)
      write(*,102)
      sigma1=0.d0
      sigma2=0.d0
      do j=1,m
      write(*,103) (alat(j,j1),j1=1,ndimfit),e(j),yeos(j)
      sigma1=sigma1+yeos(j)**2
      enddo
      write(*,104) sqrt(sigma1/(m-1))

1234  continue
!full fit
      n=nfit
      iv=0
      v=0.d0
      bounds(1,7)=-d1mach(2)
      bounds(2,7)=d1mach(2)
      call dn2fb (m,n,x,bounds,eos,iv,liv,lv,v,ui,ur,uf)

      write(*,*) 'Parabolic equation of state:         info',iv(1)
      write(*,*) 'E = x1 + x2(a-x3)^2'
      if(ndimfit.ge.2) write(*,*) '       + x4(b-x5)^2 + x6(a-x3)(b-x5)'
      if(ndimfit.ge.3) write(*,*) '       + x7(c-x8)^2 + x9(a-x3)(c-x8) + x10(b-x5)(c-x8)'
      if(ndimfit.ge.4) write(*,*) '       + x11(d-x12)^2 + x13(a-x3)(d-x12) + x14(b-x5)(d-x12) + x15(c-x8)(d-x12)'
      if(ndimfit.ge.5) write(*,*) '       + x16(e-x17)^2 + x18(a-x3)(e-x17) + x19(b-x5)(e-x17) + x20(c-x8)(e-x17)',&
         ' + x21(d-x12)(e-x17)'
      if(ndimfit.eq.6) write(*,*) '       + x22(f-x23)^2 + x24(a-x3)(f-x23) + x25(b-x5)(f-x23) + x26(c-x8)(f-x23)',&
         ' + x27(d-x12)(f-x23) + x28(e-x17)(f-x23)'
      write(*,100) (x(i),i=1,nfit)
      write(66,*) 'Parabolic equation of state:         info',iv(1)
      write(66,*) 'E = x1 + x2(a-x3)^2'
      if(ndimfit.ge.2) write(66,*) '       + x4(b-x5)^2 + x6(a-x3)(b-x5)'
      if(ndimfit.ge.3) write(66,*) '       + x7(c-x8)^2 + x9(a-x3)(c-x8) + x10(b-x5)(c-x8)'
      if(ndimfit.ge.4) write(66,*) '       + x11(d-x12)^2 + x13(a-x3)(d-x12) + x14(b-x5)(d-x12) + x15(c-x8)(d-x12)'
      if(ndimfit.ge.5) write(66,*) '       + x16(e-x17)^2 + x18(a-x3)(e-x17) + x19(b-x5)(e-x17) + x20(c-x8)(e-x17)',&
         ' + x21(d-x12)(e-x17)'
      if(ndimfit.eq.6) write(66,*) '       + x22(f-x23)^2 + x24(a-x3)(f-x23) + x25(b-x5)(f-x23) + x26(c-x8)(f-x23)',&
         ' + x27(d-x12)(f-x23) + x28(e-x17)(f-x23)'
      write(66,100) (x(i),i=1,nfit)
 100  format('Fitparameter are',(/,4f18.6))
       nf=1
      call eos(m,nfit,x,nf,yeos,ui,ur,uf)
!

      write(*,102)
      write(66,102)
 102  format(9x,'lattic parameters',7x,'energy',9x,'de(EOS)')
      sigma1=0.d0
      sigma2=0.d0
      do j=1,m
      write(*,103) (alat(j,j1),j1=1,ndimfit),e(j),yeos(j)
      write(66,103) (alat(j,j1),j1=1,ndimfit),e(j),yeos(j)
      sigma1=sigma1+yeos(j)**2
 103  format(8f15.6)
      enddo
      write(*,104) sqrt(sigma1/(m-1))
      write(66,104) sqrt(sigma1/(m-1))
 104  format('                  Sigma:',4x,2f14.6)
!
!.....create plotfile
      write(*,*) 'Optionally create data points from fit function'
      write(*,*) 'Enter number of datapoints for your',ndimfit,'dimensional Energy surface'
      write(*,*) 'NI=0 terminates; NI=1 will use 1 specific value in I-th component and allows to generate 2D-cuts'
      nni=1
      read(*,*,end=999) nni(1:ndimfit)
       if(nni(1).eq.0) stop
      do i=1,ndimfit
       write(*,*) 'Enter X-min, X-max in',i,'-th direction'
       read(*,*) xmin(i),xmax(i)
      enddo
!
       i=0
         do j1=1,nni(1)
         do j2=1,nni(2)
         do j3=1,nni(3)
         do j4=1,nni(4)
         do j5=1,nni(5)
         do j6=1,nni(6)
          i=i+1
                           alat(i,1)=xmin(1)+(xmax(1)-xmin(1))/max(nni(1)-1,1)*(j1-1)
          if(ndimfit.ge.2) alat(i,2)=xmin(2)+(xmax(2)-xmin(2))/max(nni(2)-1,1)*(j2-1)
          if(ndimfit.ge.3) alat(i,3)=xmin(3)+(xmax(3)-xmin(3))/max(nni(3)-1,1)*(j3-1)
          if(ndimfit.ge.4) alat(i,4)=xmin(4)+(xmax(4)-xmin(4))/max(nni(4)-1,1)*(j4-1)
          if(ndimfit.ge.5) alat(i,5)=xmin(5)+(xmax(5)-xmin(5))/max(nni(5)-1,1)*(j5-1)
          if(ndimfit.ge.6) alat(i,6)=xmin(6)+(xmax(6)-xmin(6))/max(nni(6)-1,1)*(j6-1)
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      m=i
      e=0.d0
      call eos(m,nfit,x,nf,yeos,ui,ur,uf)
      do i=1,ndimfit
        write(66,'("Number of points",i4," range",2f10.4)') nni(i),xmin(i),xmax(i)
      enddo
      do j=1,m
        write(66,103) (alat(j,j1),j1=1,ndimfit),yeos(j)
      enddo
! check if it is 2D and which dimensions are used for printing on plotfile
      j=0
      do i=1,ndimfit
       if(nni(i).ne.1) then
         j=j+1
         nni(j)=nni(i)
         xmin(j)=xmin(i)
         xmax(j)=xmax(i)
      endif
      enddo
      write(12,*) nni(1:j),(xmax(i)-xmin(i),i=1,j)
      write(12,111) (yeos(i)-int(x(1)),i=1,m)
 111  FORMAT(5E16.8)                                                    

!      xmin=xmin-(xmax-xmin)*0.1
!      xmax=xmax+(xmax-xmin)*0.1
!      dx=(xmax-xmin)/399
!      do j=1,400
!       emurna(j)=0.d0
!      volmurna(j)=xmin+(j-1)*dx
!      enddo
!      call murna(400,4,x,ymurna,iflag)

 999  stop
      end
!
      subroutine eos(m,n,x,nf,fvec,ui,ur,uf)
      implicit real*8 (a-h,o-z)
      parameter (ndata_points=4000)
      common alat(ndata_points,6), e(ndata_points),ndimfit
      dimension fvec(*),x(*)
      data ic/0/
      ic=ic+1
      if(ic.eq.1) print*,m,(e(i),i=1,m)
      fsum=0.d0
      do i=1,m
      fvec(i)=x(1)+x(2)*(alat(i,1)-x(3))**2 - e(i)
      enddo
      if(ndimfit.ge.2) then  
      do i=1,m
      fvec(i)=fvec(i)+x(4)*(alat(i,2)-x(5))**2  + &
                      x(6)*(alat(i,1)-x(3))*(alat(i,2)-x(5))
      enddo
      endif 
      if(ndimfit.ge.3) then  
      do i=1,m
      fvec(i)=fvec(i)+x(7)*(alat(i,3)-x(8))**2  + &
               x(9)*(alat(i,1)-x(3))*(alat(i,3)-x(8))+ &
              x(10)*(alat(i,2)-x(5))*(alat(i,3)-x(8))
      enddo
      endif
      if(ndimfit.ge.4) then  
      do i=1,m
      fvec(i)=fvec(i)+x(11)*(alat(i,4)-x(12))**2 + &
              x(13)*(alat(i,1)-x(3))*(alat(i,4)-x(12))+ &
              x(14)*(alat(i,2)-x(5))*(alat(i,4)-x(12)) + &
              x(15)*(alat(i,3)-x(8))*(alat(i,4)-x(12))
      enddo
      endif
      if(ndimfit.ge.5) then  
      do i=1,m
      fvec(i)=fvec(i)+x(16)*(alat(i,5)-x(17))**2 + &
              x(18)*(alat(i,1)-x(3))*(alat(i,5)-x(17))+ &
              x(19)*(alat(i,2)-x(5))*(alat(i,5)-x(17))+ &
              x(20)*(alat(i,3)-x(8))*(alat(i,5)-x(17)) + &
              x(21)*(alat(i,4)-x(12))*(alat(i,5)-x(17))
      enddo
      endif
      if(ndimfit.eq.6) then  
      do i=1,m
      fvec(i)=fvec(i)+x(22)*(alat(i,6)-x(23))**2 + &
              x(24)*(alat(i,1)-x(3))*(alat(i,6)-x(23))+ &
              x(25)*(alat(i,2)-x(5))*(alat(i,6)-x(23))+ &
              x(26)*(alat(i,3)-x(8))*(alat(i,6)-x(23)) + &
              x(27)*(alat(i,4)-x(12))*(alat(i,6)-x(23))+ &
              x(28)*(alat(i,5)-x(17))*(alat(i,6)-x(23))
      enddo
      endif  
      do i=1,m
      fsum=fsum+fvec(i)**2
      enddo
!      if(mod(ic,1000000).eq.5) write(*,*) ic,(x(i),i=1,n),fsum
!      write(*,*) ic,(x(i),i=1,n),fsum
      return
      end
      subroutine eosparabol(m,n,x,nf,fvec,ui,ur,uf)
      implicit real*8 (a-h,o-z)
      parameter (ndata_points=4000)
      common alat(ndata_points,6), e(ndata_points),ndimfit
      dimension fvec(*),x(*)
      data ic/0/
      ic=ic+1
      if(ic.eq.1) print*,m,(e(i),i=1,m)
      fsum=0.d0
      do i=1,m
      fvec(i)=x(1)+x(2)*(alat(i,1)-x(3))**2 - e(i)
      enddo
      if(ndimfit.ge.2) then  
      do i=1,m
      fvec(i)=fvec(i)+x(4)*(alat(i,2)-x(5))**2 ! + &
                   !   x(6)*(alat(i,1)-x(3))*(alat(i,2)-x(5))
      enddo
      endif 
      if(ndimfit.ge.3) then  
      do i=1,m
      fvec(i)=fvec(i)+x(6)*(alat(i,3)-x(7))**2  
      !  fvec(i)=fvec(i)+x(7)*(alat(i,3)-x(8))**2  + &
        !       x(9)*(alat(i,1)-x(3))*(alat(i,3)-x(8))+ &
        !      x(10)*(alat(i,2)-x(5))*(alat(i,3)-x(8))
      enddo
      endif
      if(ndimfit.ge.4) then  
      do i=1,m
      fvec(i)=fvec(i)+x(11)*(alat(i,4)-x(12))**2 + &
              x(13)*(alat(i,1)-x(3))*(alat(i,4)-x(12))+ &
              x(14)*(alat(i,2)-x(5))*(alat(i,4)-x(12)) + &
              x(15)*(alat(i,3)-x(8))*(alat(i,4)-x(12))
      enddo
      endif
      if(ndimfit.ge.5) then  
      do i=1,m
      fvec(i)=fvec(i)+x(16)*(alat(i,5)-x(17))**2 + &
              x(18)*(alat(i,1)-x(3))*(alat(i,5)-x(17))+ &
              x(19)*(alat(i,2)-x(5))*(alat(i,5)-x(17))+ &
              x(20)*(alat(i,3)-x(8))*(alat(i,5)-x(17)) + &
              x(21)*(alat(i,4)-x(12))*(alat(i,5)-x(17))
      enddo
      endif
      if(ndimfit.eq.6) then  
      do i=1,m
      fvec(i)=fvec(i)+x(22)*(alat(i,6)-x(23))**2 + &
              x(24)*(alat(i,1)-x(3))*(alat(i,6)-x(23))+ &
              x(25)*(alat(i,2)-x(5))*(alat(i,6)-x(23))+ &
              x(26)*(alat(i,3)-x(8))*(alat(i,6)-x(23)) + &
              x(27)*(alat(i,4)-x(12))*(alat(i,6)-x(23))+ &
              x(28)*(alat(i,5)-x(17))*(alat(i,6)-x(23))
      enddo
      endif  
      do i=1,m
      fsum=fsum+fvec(i)**2
      enddo
!      if(mod(ic,1000000).eq.5) write(*,*) ic,(x(i),i=1,n),fsum
!      write(*,*) ic,(x(i),i=1,n),fsum
      return
      end
      subroutine eosother(m,n,x,nf,fvec,ui,ur,uf)
      implicit real*8 (a-h,o-z)
      parameter (ndata_points=4000)
      common alat(ndata_points,6), e(ndata_points),ndimfit
      common /other/ xc(28)
      dimension fvec(*),x(*)
      data ic/0/
      ic=ic+1
      if(ic.eq.1) print*,m,(e(i),i=1,m)
      fsum=0.d0
      do i=1,m
      fvec(i)=xc(1)+xc(2)*(alat(i,1)-xc(3))**2 - e(i)
      enddo
      if(ndimfit.ge.2) then  
      do i=1,m
      fvec(i)=fvec(i)+xc(4)*(alat(i,2)-xc(5))**2  + &
                      x(1)*(alat(i,1)-xc(3))*(alat(i,2)-xc(5))
      enddo
      endif 
      if(ndimfit.ge.3) then  
      do i=1,m
      fvec(i)=fvec(i)+xc(6)*(alat(i,3)-xc(7))**2 +& 
               x(2)*(alat(i,1)-xc(3))*(alat(i,3)-xc(7))+ &
              x(3)*(alat(i,2)-xc(5))*(alat(i,3)-xc(7))
      enddo
      endif
      if(ndimfit.ge.4) then  
      do i=1,m
      fvec(i)=fvec(i)+x(11)*(alat(i,4)-x(12))**2 + &
              x(13)*(alat(i,1)-x(3))*(alat(i,4)-x(12))+ &
              x(14)*(alat(i,2)-x(5))*(alat(i,4)-x(12)) + &
              x(15)*(alat(i,3)-x(8))*(alat(i,4)-x(12))
      enddo
      endif
      if(ndimfit.ge.5) then  
      do i=1,m
      fvec(i)=fvec(i)+x(16)*(alat(i,5)-x(17))**2 + &
              x(18)*(alat(i,1)-x(3))*(alat(i,5)-x(17))+ &
              x(19)*(alat(i,2)-x(5))*(alat(i,5)-x(17))+ &
              x(20)*(alat(i,3)-x(8))*(alat(i,5)-x(17)) + &
              x(21)*(alat(i,4)-x(12))*(alat(i,5)-x(17))
      enddo
      endif
      if(ndimfit.eq.6) then  
      do i=1,m
      fvec(i)=fvec(i)+x(22)*(alat(i,6)-x(23))**2 + &
              x(24)*(alat(i,1)-x(3))*(alat(i,6)-x(23))+ &
              x(25)*(alat(i,2)-x(5))*(alat(i,6)-x(23))+ &
              x(26)*(alat(i,3)-x(8))*(alat(i,6)-x(23)) + &
              x(27)*(alat(i,4)-x(12))*(alat(i,6)-x(23))+ &
              x(28)*(alat(i,5)-x(17))*(alat(i,6)-x(23))
      enddo
      endif  
      do i=1,m
      fsum=fsum+fvec(i)**2
      enddo
      if(mod(ic,1000000).eq.1) write(*,*) ic,(xc(i),i=1,7)
      write(*,*) ic,(x(i),i=1,n),fsum
      return
      end
