      subroutine nose0(x0,xplus1,xnos0,xnospl,n,fwert, &
                      fnose,temp,delta,weight,mult)
      implicit double precision (a-h,o-z)
      include 'param.inc'
      COMMON /STRUK/  AA,BB,CC,ALPHA(3),PIA(3),VOL
      DIMENSION  x0(n),xplus1(n),fwert(n)
      DIMENSION  weight(*),mult(*)
!
!     initialize nose 
!
      fnose=fnose*1.0d12/(4.1341373358254d16)
      boltzk=3.166678911d-6
      API=4.D0*DATAN(1.D0)
      ekin=0.0d0
      xnos0=0.0d0
      free=-6.0d0
      do 10 i=1,n,3
         help=dble(i+2)/3.0d0
         iatom=nint(help)
         do 20 j=1,mult(iatom)
            free=free+3.0d0
  20     CONTINUE
  10  CONTINUE
      if (free.lt.0.1d0) free=1.0d0
      QNOSE=2.0d0/((2.0d0*API*FNOSE)**2)*2.0d0*FREE*boltzk*TEMP
      VNOS=delta/qnose*(ekin-FREE*boltzk*TEMP/2.0d0)
      facnos=1.0d0/(1.0d0+VNOS*delta/2.0d0)
      write(6,*) 'fnose,delta,qnose,free',fnose,delta,qnose,free
!
!     determine anner
!
      do 30 i=1,100
         xnospl=(delta**2)/qnose* &
                (2.0d0*ekin*(facnos**2)-FREE*boltzk*TEMP)
         VNOS=(xnospl)/2.0d0/delta
         facnos=1.0d0/(1.0d0+VNOS*delta/2.0d0)
         write(6,*) 'nose:i,vnos,facnos',i,vnos,facnos 
   30 CONTINUE
      ekin=ekin*(facnos**2)
      VNOS=(xnospl)/2.0d0/delta
      anner=(VNOS*delta)/2.0d0
      write(6,*) 'nose:ekin,free,anner',ekin,free,anner
!
!     force and gradient in Ry/bohr
!     determine new position
!
      iatom=1
      ihelp=0
      do 40 i=1,n
         ihelp=ihelp+1
         xplus1(i)=x0(i)/(1.0d0+anner)+ &
                   fwert(i)/2.0d0* &
                   (delta**2)/weight(iatom)/(1.0d0+anner)
!         write(6,'(10x,''(nwtmin)  i,damp,speed,force =''
!     &      ,i3,f6.2,2f10.4)') i,damp,speed,fwert(i)
         if (ihelp.eq.3) then
            ihelp=0
            iatom=iatom+1
         endif
 40   CONTINUE
      write(6,*)'  old position',x0(1),x0(2),x0(3)
      write(6,*)'  new position'
      call ininos(x0,xplus1,delta,ekin,n,weight,mult)
      temp=2.0d0*ekin/boltzk/free
      write(6,*)(xplus1(i),i=1,n)
      write(6,*)'fwert',(fwert(i),i=1,n)
      write(6,*)'new ekin,temp,delta',ekin,temp,delta
      return
      end
