      program u4
      
      implicit real*8 (k,x,y,z)
      character*3 katom
      character*3 matrix
      character*80 titel,zitel1
      integer atome,ergeb,zsave1,zsave2,faktor
!     a,b,c Gitterkonstanten      
      real*8 a,b,c,alpha,beta,gamma,rbon,dist0
      real*8 theta,phi,rot,magnf,persp,hight,rmi,rma
      real*8 x(0:1000),y(0:1000),z(0:1000),dista(0:1000)
      dimension ikern(0:1000),rad(0:1000),mult(0:1000)
      integer m,n,iclc,ibft,nshx,nshy,icon,initsv,nbco,kernz
      dimension nuca(1000),nuce(1000),rmin(1000),rmax(1000),rbond(1000)
      dimension nbcol(1000),kernz(0:1000)
      call o_file
!     oeffnet die noetigen files (fehlerkontrolle)!                
!      open(unit=10,file='case.struct',form='formatted',status='old')
!      open (unit=20,file='case.plt',form='formatted',status='unknown')
      write(*,*) ' Please enter the KEY-atom number:'
      read(*,*) iatom
      write(katom,'(i3)') iatom
      write(*,*) 'katom:',katom
      atome=1
      read(10,100) titel1,matrix
      read(10,*)
      read(10,*) a,b,c,alpha,beta,gamma
 1    read(10,90,end=999) titel
      if(.not.(titel(12:16).eq.'EQUIV'.and.titel(7:9).eq.katom)) goto 1
      read(titel(36:66),91) x(1),y(1),z(1)
      rad(1)=a/10.d0
      kernz(1)=iatom
 91   format(3f10.5)
      write(*,*) x(1),y(1),z(1),kernz(1)
      read(10,*)
      read(10,*)
 2    read(10,90,end=3) titel
      if(titel(2:2).eq.' ') goto 3
      atome=atome+1
      rad(atome)=a/10.d0
      read(titel,92) kernz(atome),x(atome),y(atome),z(atome), &
        dista(atome)
 92   format(6x,i3,14x,3f8.4,3x,f9.5)
      write(*,92) kernz(atome),x(atome),y(atome),z(atome), &
        dista(atome)
      goto 2
 3    write (*,*)matrix
      int=atome
      if(gamma.eq.0.)then
      gamma=90.
      endif
  90  format(a80)
 100  format(a80,/,a3)     
      close(10,status='keep')
      
      call bmatr(gamma,matrix,int,a,b,c,x,y,z)
      
      ergeb=int
      
      write (20,190) titel1
      write(20,200) ergeb,a,b,c
      do 20 j=1,int
!      ikern(j)=idint(kernz(j))
      write(20,210) x(j),y(j),z(j),rad(j),kernz(j) 
 20   continue     
      theta=75.
      phi=-70.
      rot=0.
      magnf=1.
      persp=16.875
      hight=-0.5
      write(20,220) theta,phi,rot,magnf,persp,hight
      m=19
      n=1
      iclc=1
      ibft=0
      nshx=0
      nshy=0
      initsv=0

      icon=1
      write(20,230) m,n,iclc,ibft,nshx,nshy,icon,initsv
      do 400 i=1,icon
      rbond(i)=a/100.d0
      nbcol(i)=4
      write(20,240) kernz(1),kernz(2),dista(2)-0.01,dista(2)+0.01, &
      rbond(i),nbcol(i)
  400 continue   
      close(20,status='keep')
      
 
 190  format(a72)
 200  format(i5,4f15.9)     
 210  format(4f15.9,i4)
 220  format(6f12.5)   
 230  format(8i5)     
 240  format(2i5,3d15.7,i5)
      stop      
 999  write(*,*) 'end of file reached'
      end

      subroutine bmatr(gamma,matrix,int,a,b,c,x,y,z)
      
      implicit real*8 (k,x,y,z)
      character*3 matrix
!     a,b,c Gitterkonstanten      
      real*8 a,b,c,gamma,pi
      
      dimension x(0:1000),y(0:1000),z(0:1000)
      integer int
      if (matrix(1:1).eq.'P') then
        pi=acos(-1d0)
        gamma=gamma*pi/180.d0
        do 75 i=1,int
        
        x1=(a*sin(gamma))   
        x2=a*cos(gamma)
        x3=0.
        y1=0.
        y2=b
        y3=0.
        z1=0.
        z2=0.
        z3=c
        xx=x1*x(i)+y1*y(i)+z1*z(i)
        yy=x2*x(i)+y2*y(i)+z2*z(i)
        zz=x3*x(i)+y3*y(i)+z3*z(i)
!        xx=x1*x(i)+x2*y(i)+x3*z(i)
!        yy=y1*x(i)+y2*y(i)+y3*z(i)
!        zz=z1*x(i)+z2*y(i)+z3*z(i)
        x(i)=xx
        y(i)=yy
        z(i)=zz
   75  continue
        endif
      
      
      if (matrix(1:1).eq.'R') then
        do 80 i=1,int
        x1=(a/sqrt(3.d0)/2.d0)   
        x2=-(a/2)
        x3=c/3
        y1=(a/sqrt(3.d0)/2.d0)
        y2=a/2
        y3=c/3
        z1=-(a/sqrt(3.d0))
        z2=0
        z3=c/3
        xx=x1*x(i)+y1*y(i)+z1*z(i)
        yy=x2*x(i)+y2*y(i)+z2*z(i)
        zz=x3*x(i)+y3*y(i)+z3*z(i)
        x(i)=xx
        y(i)=yy
        z(i)=zz
!        x(i)=x1*x(i)+x2*y(i)+x3*z(i)
!        y(i)=y1*x(i)+y2*y(i)+y3*z(i)
!        z(i)=z1*x(i)+z2*y(i)+z3*z(i)
   80  continue
        endif
        if (matrix(1:1).eq.'H') then
        do 85 i=1,int
        write(*,*) x(i),y(i),z(i)
        x1=dsqrt(3.d0)*a/2
        x2=-(a/2)
        x3=0
        y1=0
        y2=a
        y3=0
        z1=0
        z2=0
        z3=c
        xx=x1*x(i)+y1*y(i)+z1*z(i)
        yy=x2*x(i)+y2*y(i)+z2*z(i)
        zz=x3*x(i)+y3*y(i)+z3*z(i)
!        xx=x1*x(i)+x2*y(i)+x3*z(i)
!        yy=y1*x(i)+y2*y(i)+y3*z(i)
!        zz=z1*x(i)+z2*y(i)+z3*z(i)
        x(i)=xx
        y(i)=yy
        z(i)=zz
   85   continue 
        endif

        if (matrix(1:1).eq.'F'.or.matrix(1:1).eq.'B'.or. &
              matrix(1:1).eq.'C') then
        do 700 i=1,int
        x1=a
        x2=0
        x3=0
        y1=0
        y2=b
        y3=0
        z1=0
        z2=0
        z3=c
        xx=x1*x(i)+y1*y(i)+z1*z(i)
        yy=x2*x(i)+y2*y(i)+z2*z(i)
        zz=x3*x(i)+y3*y(i)+z3*z(i)
!        xx=x1*x(i)+x2*y(i)+x3*z(i)
!        yy=y1*x(i)+y2*y(i)+y3*z(i)
!        zz=z1*x(i)+z2*y(i)+z3*z(i)
        x(i)=xx
        y(i)=yy
        z(i)=zz
  700    continue 
        endif

      return
      end
      subroutine o_file
      character*80 namstr,namlat
      character*10 test
      write(*,*)'Filename without extension: '
      read(*,99900)namstr
      do 1000 i=80,1,-1
1000  if(namstr(i:i).ne.' ') goto 1010
1010  if(i.gt.76)goto 9992
      namlat=namstr      
      namstr(i+1:i+9)='.outputnn'
      namlat(i+1:i+4)='.plt'
      open(unit=10,file=namstr,form='formatted',status='old',err=9993)
      open(unit=20,file=namlat,form='formatted',status='new',err=9994)
      return
9992  write(*,*)'Filename too long' 
      stop      
9993  write(*,*)'struct-File does not exist'      
      stop
9994  write(*,*)'File ',namlat(1:i+4),' already exists.'
      write(*,*)' Overwrite (y/n)?'
      read(*,99901)test
      if(test.eq.'y' .or. test.eq.'Y') then
       open(unit=20,file=namlat,form='formatted')
       return
      else
       stop
      endif    
99900 format(a80) 
99901 format(a1)     
      end
      
