      subroutine fcore(ja,nr,r,rho)
!     input:  non-spherical potential <- tape19
!     output: [yu91]-force due to core-density -> tape73

      implicit real*8 (a-h,o-z)
      INCLUDE 'param.inc'

      integer    ir,ik
      integer    nr
      integer    lmmax,ja
      logical    forout
      real*8     r(nrad),v1mp(nrad),dv1mp(nrad)
      real*8     value(nrad),rho(nrad)
      real*8     sq4pi3,fint,cint
      real*8     fcor(3),f1(3),f2(3)
      PI = 4.0D+0*ATAN(1.0D+0)
      sq3=sqrt(3.d0)
      sq4pi=sqrt(4.d0*pi)
      sq4pi3=sqrt(4.d0*pi/3.d0)
      forout=.false.
      fcor(1)=0.d0
      fcor(2)=0.d0
      fcor(3)=0.d0

!.... read total nonspherical potential of tape19=vns                   
!     norm of vlmp=v(lmp)*1           , if lmp <> 00+
!     norm of vlmp=v(00+)*r*sqrt(4*pi), if lmp =  00+
!     norm of rho =c(00+)*r*r*sqrt(4*pi)
      
      read(19,'(3x)')
      read(19,'(15x,i3//)') lmmax
      DO 5000 lm1=1,lmmax
        read(19,'(15x,i3,5x,i2/)') lll,mm
        ll =iabs(lll)                                                   
        read(19,'(3x,4e19.12)') ( v1mp(ir), ir=1,nr)                   
        read(19,'(/)')                                                     

!     force-calculation

        if (ll.eq.1) then
        forout=.true.
!     lll,mm --> ik (component of force)
!     +1,1 --> 1; -1,1 --> 2; +1,0 --> 3
              ik=mod(lll+mm+1,3)+1

!     normalize  v1mp()=v1mp()/r()

          dx=log(r(nr)/r(1))/(nr-1)
          v1mp(1)=v1mp(1)/r(1)
          DO 5001 ir=2,nr
            v1mp(ir)=v1mp(ir)/r(ir)
 5001     CONTINUE
          
!     calculate derivative of v/r
          call dfradial(r,v1mp,dv1mp,nr)

!     calculate radial integral
          DO 5002 ir=1,nr
            value(ir)=rho(ir)*v1mp(ir)
 5002     CONTINUE
          call charge(r,dx,value,1,nr,fint)
          f1(ik)=-fint/sq4pi*sq3

          DO 5003 ir=1,nr
            value(ir)=rho(ir)*r(ir)*dv1mp(ir)
 5003     CONTINUE
          call charge(r,dx,value,1,nr,fint)
          f2(ik)=-fint/sq4pi/sq3

          DO 5004 ir=1,nr
            value(ir)=rho(ir)
 5004     CONTINUE
          call charge(r,dx,value,1,nr,cint)
!     ik. component of total core-force
          fcor(ik)=f1(ik)+f2(ik)

        endif
        
 5000 CONTINUE
      read(19,'(///)')                                                     
      if (forout) then
      write(6,'(2i3,a7,4e15.7)') ja,1,'> COR' &
        ,sqrt(fcor(1)**2+fcor(2)**2+fcor(3)**2) &
        ,(fcor(ik),ik=1,3)
      write(21,700) 
      write(21,710) ja,ja,sqrt(fcor(1)**2+fcor(2)**2+fcor(3)**2)*1000 &
        ,fcor(1)*1000,fcor(2)*1000,fcor(3)*1000
      endif
 700  FORMAT (7x,'CORE-FORCE in mRy/a.u. = |F|',6x,'Fx',13x, &
              'Fy',13x,'Fz')
 710  FORMAT (':FCO',i3.3,':',1x,i3,'.ATOM',4f15.3,/)
!        write(83,'(6x,a7,4e15.7)') 
!     &    '    1',sqrt(f1(1)**2+f1(2)**2+f1(3)**2),(f1(ik),ik=1,3),
!     &    '+   2',sqrt(f2(1)**2+f2(2)**2+f2(3)**2),(f2(ik),ik=1,3)
!        write(83,'(2i3,a7,4e15.7)') ja,1,'> COR'
!     &    ,sqrt(fcor(1)**2+fcor(2)**2+fcor(3)**2)
!     &    ,(fcor(ik),ik=1,3)
      return
      end
