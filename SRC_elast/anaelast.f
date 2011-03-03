!******* This program calculates the cubic elastic tensor
!******* from the results of the calculations initiated
!******* by init_elast and elast_setup
!*******
!*******   Author: CHARPIN Thomas
!*******           Laboratoire des Geomateriaux de l'IPGP
!*******           4,pl Jussieu
!*******           F-75252 Paris Cedex 05
!*******           France
!*******           charpin@ipgp.jussieu.fr   
!*******
!*******


      PROGRAM anaelast
      CHARACTER*1 lat,bul
      INTEGER noa,i,neos,ntetra,nrhomb
      CHARACTER*79 row,title
      REAL*8 a0,a,cr0,eeos(100),etetra(100),erhomb(100)
      REAL*8 vol0,vol(100),testrain(100),rstrain(100)
      REAL*8 p,kt,res,resr,c11,c12,c44,conv

      PARAMETER (conv=14710.5164057d0)



!******* We first read unstrained state parameters...

      OPEN(unit=10,file='init.struct',status='old',err=8000)
      READ(10,1000) title
      READ(10,1010) lat,noa
      READ(10,1000) row
      READ(10,1020) a0
      CLOSE(10)

      
      vol0=a0**3.d0
      IF ( lat.eq.'F') then
         vol0=vol0/4.d0
         cr0=a0*(3.d0**(0.5d0))
      ELSEIF(lat.eq.'B') then
         vol0=vol0/2.d0
         cr0=a0*(3.d0**(0.5d0))/2.d0
      ELSEIF(lat.eq.'P') then
         cr0=a0*(3.d0**(0.5d0))
      ELSE
         GOTO 8000
      ENDIF


      WRITE(*,*) '***************************************'      
      WRITE(*,*) '***************************************'
      WRITE(*,*) ' We are calculating elastic tensor for'
      WRITE(*,"(9X,A1,A13)") lat,' Cubic phase:'
      WRITE(*,"(3X,A79)") title
      WRITE(*,*) '         At volume:'
      WRITE(*,"(F9.2,A19)") vol0,' bohr^3 per formula'
      WRITE(*,*) '***************************************'
      WRITE(*,*) '***************************************'
      WRITE(*,*) ' '
      WRITE(*,*) ' '

!******* ...then we calculate strains for all the calculated 
!******* structures using eos.lat, rhomb.lat,tetra.lat...

      OPEN(unit=20,file='eos.lat',status='old',err=8000)
      READ(20,1050) neos
      IF (neos.le.4) then
         WRITE(*,*) 'Not enough points for Birch-Murnaghan fit'
         GOTO 8000
      ENDIF
      DO i=1,neos
         READ(20,1060) a
         vol(i)=a**3.d0
      ENDDO
      CLOSE(20)

      IF ( lat.eq.'F') then
         DO i=1,neos
            vol(i)=vol(i)/4.d0
         ENDDO
      ELSEIF(lat.eq.'B') then
         DO i=1,neos
            vol(i)=vol(i)/2.d0
         ENDDO         
      ENDIF

      OPEN(unit=30,file='tetra.lat',status='old',err=8000)
      READ(30,1050) ntetra
      DO i=1,ntetra
         READ(30,1070) a
         testrain(i)=(a-a0)/a0
!         WRITE(*,"(F12.6)") testrain(i)
      ENDDO
      CLOSE(30)
      
      OPEN(unit=40,file='rhomb.lat',status='old',err=8000)
      READ(40,1050) nrhomb
      DO i=1,nrhomb
         READ(40,1080) a
         rstrain(i)=(a-cr0)/cr0
!         WRITE(*,"(F12.6)") rstrain(i)
      ENDDO
      CLOSE(40)

!******* ...then we read energy using eos.ene, rhomb.ene, tetra.ene

      OPEN(unit=20,file='eos.ene',status='old',err=8000)
      READ(20,1000) row
      DO i=1,neos
         READ(20,1090) eeos(i)
      ENDDO
      CLOSE(20)

      OPEN(unit=30,file='tetra.ene',status='old',err=8000)
      READ(30,1000) row
      DO i=1,ntetra
         READ(30,1100) etetra(i)
!         WRITE(*,"(F12.6)") etetra(i)
      ENDDO
      CLOSE(30)

      OPEN(unit=40,file='rhomb.ene',status='old',err=8000)
      READ(40,1000) row
      DO i=1,nrhomb
         READ(40,1100) erhomb(i)
!         WRITE(*,"(F12.6)") erhomb(i)
      ENDDO
      CLOSE(20)

!****** ...then we fit...
!****** (for 'eos' we use third order Birch-Murnaghan fit.
!******  for others we use polynomial fit)
!******

      CALL fiteos(vol,eeos,vol0,kt,p,neos)
      WRITE(*,2000)
      READ(*,*)
      CALL fittetra(testrain,etetra,vol0,res,ntetra)
      WRITE(*,2000)
      READ(*,*) 
      CALL fitrhomb(rstrain,erhomb,vol0,resr,nrhomb)
      WRITE(*,2000)
      READ(*,*)

!******
!****** ...Now we shall the solve system:
!******  kt=(c11+2*c12)/3
!******  res=c11-c12
!******  resr=c11+2*c12+4*c44
!******      =3*kt+4c44
!******


      c44=(resr-3.d0*kt)/4.d0
      WRITE(*,"(A5,F12.6,a9,f12.3,a4)") 'c44=',c44,' a.u.  or',c44*conv,' GPa'
      c11=kt+res*(2.d0/3.d0)
      WRITE(*,"(A5,F12.6,a9,f12.3,a4)") 'c11=',c11,' a.u.  or',c11*conv,' GPa' 
      c12=c11-res
      WRITE(*,"(A5,F12.6,a9,f12.3,a4)") 'c12=',c12,' a.u.  or',c12*conv,' GPa'

!********* ....system solved...lets flush the results 
      
      OPEN(unit=20,file='elast.output')
      WRITE(20,1110) lat,title
      WRITE(20,1120) vol0
      WRITE(20,1140) p,p*conv
      WRITE(20,*) 'in atomic units: '
      WRITE(20,1130) c11,c12,c44
      WRITE(20,*) ' '
      WRITE(20,*) 'in GPa: '
      WRITE(20,1130) c11*conv,c12*conv,c44*conv

      CLOSE(20)



 1000 FORMAT(A79)
 1010 FORMAT(A1,26X,I3)
 1020 FORMAT(F10.6)
 1050 FORMAT(4X,I2)
 1060 FORMAT(40X,F9.5)
 1070 FORMAT(42X,F9.5)
 1080 FORMAT(60X,F9.5)
 1090 FORMAT(53X,F20.6)
 1100 FORMAT(55X,F20.6)
 1110 FORMAT('Elastic parameters for cubic ',A1,' phase',/,A79)
 1120 FORMAT('At volume:',F18.5,' bohr^3 per formula:')
 1130 FORMAT('c11=',F16.6,' c12=',F16.6,' c44=',F16.6)
 1140 FORMAT('At calculated pressure: ',F8.6,' a.u. or ',F12.6,' GPa',/)
 2000 FORMAT('Hit return to continue',A1)

      GOTO 8001
 8000 write(*,*) 'Error in anaelast'

      STOP

 8001 write(*,*) 'Analyze done.....'

      END

      





      SUBROUTINE fitrhomb(rstrain,erhomb,vol0,res,n)
      REAL*8 rstrain(100),erhomb(100),vol0,res,yp(100)
      REAL*8 emin,emax,pas,eo,val,eps
      INTEGER n,ndeg,ierr,i
      REAL*8,allocatable :: x(:),y(:),w(:),r(:),a(:)
      allocate( x(n),y(n),w(n),r(n),a(6*n+3))


      DO i=1,n
         x(i)=rstrain(i)
         y(i)=erhomb(i)
         w(i)=-1.d0
      ENDDO

      eps=1.d-16

      CALL dpolft(n,x,y,w,n-1,ndeg,eps,r,ierr,a)
      CALL dp1vlu(ndeg,2,0.d0,val,yp,a)

      res=3.d0*yp(2)/vol0

      emin=rstrain(1)
      emax=emin
      DO i=1,n
         IF (rstrain(i).ge.emax) emax=rstrain(i)
         IF (rstrain(i).le.emin) emin=rstrain(i)
      ENDDO      

      OPEN(unit=20,file='rhomb.fit')   
      pas=(emax-emin)/100.d0
      DO i=1,101
         eo=emin+(i-1)*pas
         CALL dp1vlu(ndeg,0,eo,val,yp,a)
         write(20,800) eo,val
      ENDDO
      CLOSE(20)


      
      OPEN(unit=40,file='rhomb.strain')
      DO i=1,n
         WRITE(40,800) rstrain(i),erhomb(i)
      ENDDO
      CLOSE(40)


      OPEN(unit=30,file='rhomb.output')
      WRITE(30,"(a)") 'Polynomial fit for rhombohedral strain done'
      WRITE(*,*) '******************************************'
      WRITE(*,"(a)") 'Polynomial fit or rhombohedral strain done'
      WRITE(30,910) eps,ndeg
      WRITE(*,910) eps,ndeg
      WRITE(30,900) vol0
      WRITE(*,900) vol0
      WRITE(30,920) res,res*14710.5164057d0
      WRITE(*,920) res,res*14710.5164057d0
      WRITE(*,*) '****************************************'
      WRITE(*,*) ' '
      WRITE(30,*) ' '
      WRITE(30,930)
      DO i=1,n
         WRITE(30,940) rstrain(i),erhomb(i),r(i)-erhomb(i)
      ENDDO

      CLOSE(30)

      
      deallocate( x,y,w,r,a)

 800  FORMAT(F9.6,1X,F15.6)
 900  FORMAT(/,'At volume= ',F9.2,' bohr^3')
 910  FORMAT('A RMS of ',E12.6,' was achieved using a polynome of degree : ',I2)
 920  FORMAT('C11+2C12+4C44 is: ',F8.6,' a.u or ',F9.3,' GPa')
 930  FORMAT(9x,'Strain',7x,'energy',10x,'dE')
 940  FORMAT(F14.8,F14.6,E14.6)

      RETURN
      END      






      SUBROUTINE fittetra(testrain,etetra,vol0,res,n)
      REAL*8 testrain(100),etetra(100),vol0,res,yp(100)
      REAL*8 emin,emax,pas,eo,val,eps
      INTEGER n,ndeg,ierr,i
      REAL*8,allocatable :: x(:),y(:),w(:),r(:),a(:)
      allocate( x(n),y(n),w(n),r(n),a(6*n+3))


      DO i=1,n
         x(i)=testrain(i)
         y(i)=etetra(i)
         w(i)=-1.d0
      ENDDO

      eps=1.d-16

      CALL dpolft(n,x,y,w,n-1,ndeg,eps,r,ierr,a)
      CALL dp1vlu(ndeg,2,0.d0,val,yp,a)

      res=yp(2)/vol0/6.d0

      emin=testrain(1)
      emax=emin
      DO i=1,n
         IF (testrain(i).ge.emax) emax=testrain(i)
         IF (testrain(i).le.emin) emin=testrain(i)
      ENDDO      

      OPEN(unit=20,file='tetra.fit')   
      pas=(emax-emin)/100.d0
      DO i=1,101
         eo=emin+(i-1)*pas
         CALL dp1vlu(ndeg,0,eo,val,yp,a)
         write(20,800) eo,val
      ENDDO
      CLOSE(20)


      
      OPEN(unit=40,file='tetra.strain')
      DO i=1,n
         WRITE(40,800) testrain(i),etetra(i)
      ENDDO
      CLOSE(40)


      OPEN(unit=30,file='tetra.output')
      WRITE(30,"(a)") 'Polynomial fit or tetragonal strain done'
      WRITE(*,*) '****************************************'
      WRITE(*,"(a)") 'Polynomial fit or tetragonal strain done'
      WRITE(30,910) eps,ndeg
      WRITE(*,910) eps,ndeg
      WRITE(30,900) vol0
      WRITE(*,900) vol0
      WRITE(30,920) res,res*14710.5164057d0
      WRITE(*,920) res,res*14710.5164057d0
      WRITE(*,*) '****************************************'
      WRITE(*,*) ' '
      WRITE(30,*) ' '
      WRITE(30,930)
      DO i=1,n
         WRITE(30,940) testrain(i),etetra(i),r(i)-etetra(i)
      ENDDO

      CLOSE(30)
      

      deallocate( x,y,w,r,a)

      

 800  FORMAT(F9.6,1X,F15.6)
 900  FORMAT(/,'At volume= ',F9.2,' bohr^3')
 910  FORMAT('A RMS of ',E12.6,' was achieved using a polynome of degree : ',I2)
 920  FORMAT('C11-C12 is: ',F8.6,' a.u or ',F9.3,' GPa')
 930  FORMAT(9x,'Strain',7x,'energy',10x,'dE')
 940  FORMAT(F14.8,F14.6,E14.6)

      RETURN
      END



      
      SUBROUTINE fiteos(volu,ener,vol0,kt,p,m)
      implicit real*8 (a-h,o-z)
      REAL*8 volu(100),ener(100),p,kt,vol0,sigma
      INTEGER i,j,h,m,n
      common vol(400), e(400)
      dimension  x(4),y(400),diag(100),fjac(100,4)
      dimension ipvt(4),qtf(4),wa1(4),wa2(4),wa3(4),wa4(100)
      external birch,birchone,bulk

      n=4
      ftol=1.d-16
      maxfev=100000
      epsfcn=0.001d0
      mode=1
      factor=100.d0
      nprint=0
      ldfjac=100

      DO i=1,m
         vol(i)=volu(i)
         e(i)=ener(i)
      ENDDO

      x(1)=e(1)
      x(2)=1.d0
      x(3)=-10.d0
      x(4)=100.d0
      call lmdif(birch,m,n,x,y,ftol,ftol,ftol,maxfev,epsfcn,diag,mode, &
      factor,nprint, &
       info,nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

      p=pressure(n,x,vol0)
      kt=bulk(n,x,vol0)

      vmin=vol(1)
      vmax=vmin
      DO i=1,m
         IF (vol(i).ge.vmax) vmax=vol(i)
         IF (vol(i).le.vmin) vmin=vol(i)
      ENDDO
!      write(*,*) vmin,vmax

      

      OPEN(unit=20,file='eos.fit')   
      pas=(vmax-vmin)/100.d0
      DO i=1,101
         vo=vmin+(i-1)*pas
         write(20,800) vo,birchone(n,x,vo)
      ENDDO
      CLOSE(20)
      
      OPEN(unit=40,file='eos.strain')
      DO i=1,m
         WRITE(40,800) vol(i),e(i)
      ENDDO
      CLOSE(40)


      OPEN(unit=30,file='eos.output')
      WRITE(30,"(a)") 'Birch Murnaghan fit done'
      WRITE(*,*) '*****************************'
      WRITE(*,"(a)") 'Birch Murnaghan fit done'
      WRITE(30,900) vol0
      WRITE(*,900) vol0
      WRITE(30,910) p,p*14710.5164057d0
      WRITE(*,910) p,p*14710.5164057d0
      WRITE(30,920) kt,kt*14710.5164057d0
      WRITE(*,920) kt,kt*14710.5164057d0
      WRITE(*,*) '*****************************'
      WRITE(*,*) ' '
      WRITE(30,*) ' '
      WRITE(30,930)
      sigma=0.d0
      DO i=1,m
         WRITE(30,940) vol(i),e(i),y(i)
         sigma=sigma+(y(i)**2.d0)
      ENDDO
      WRITE(30,950) sigma
      CLOSE(30)

 800  FORMAT(F15.5,1X,F15.6)
 900  FORMAT(/,'At volume= ',F9.2,' bohr^3')
 910  FORMAT('Pressure is: ',F8.6,' a.u. or ',F7.3,' GPa')
 920  FORMAT('Bulk modulus is: ',F8.6,' a.u or ',F9.3,' GPa=(C11+2C12)/3 ')
 930  FORMAT(9x,'vol',7x,'energy',10x,'dE')
 940  FORMAT(F14.5,2F14.6)
 950  FORMAT('Sigma:',22x,E14.6)


      RETURN
      END


      SUBROUTINE birch(m,n,x,fvec,iflag)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER m,n,iflag
      REAL*8 fvec(m),x(n)
      COMMON vol(400),e(400)
      

      DO i=1,m
         v=vol(i)**(-1.d0/3.d0)
        fvec(i)=x(1)+x(2)*v**(2.d0)+x(3)*v**(4.d0)+x(4)*vol(i)**(-2.d0)
        fvec(i)=fvec(i)-e(i)
      ENDDO        
        
      RETURN
      END


      FUNCTION pressure(n,x,vo)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      REAL*8 x(n)

      v=vo**(-1.d0/3.d0)
      pressure=(2.d0/3.d0)*x(2)*v**(5.d0)+(4.d0/3.d0)*x(3)*v**(7.d0)+2.d0*x(4)*vo**(-3.d0)
      
!      pressure=pressure*14710.5164057d0

      return
      end



      FUNCTION bulk(n,x,vo)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      REAL*8 x(n)

      v=vo**(-1.d0/3.d0)
      bulk=(10.d0/9.d0)*x(2)*v**(5.d0)+(28.d0/9.d0)*x(3)*v**(7.d0)+6.d0 &
           *x(4)*vo**(-3.d0)

      return
      end



      FUNCTION birchone(n,x,vo)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      REAL*8 x(n)
      
      v=vo**(-1.d0/3.d0)
      birchone=x(1)+x(2)*v**(2.d0)+x(3)*v**(4.d0)+x(4)*vo**(-2.d0)
      
        
      RETURN
      END

