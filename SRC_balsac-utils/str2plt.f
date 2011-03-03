      program u4
      
      implicit real*8 (k,x,y,z)
      character*3 matrix
      character*10 titel
      integer atome,ergeb,zsave1,zsave2,faktor
!     a,b,c Gitterkonstanten      
      real*8 a,b,c,alpha,beta,gamma,rbon,dist0
      real*8 theta,phi,rot,magnf,persp,hight,rmi,rma
      real*8 x(0:1000),y(0:1000),z(0:1000)
      dimension ikern(0:1000),rad(0:1000),mult(0:1000),kernz(0:1000)
      integer m,n,iclc,ibft,nshx,nshy,icon,initsv,nbco
      dimension nuca(1000),nuce(1000),rmin(1000),rmax(1000),rbond(1000)
      dimension nbcol(1000)
      call o_file
!     oeffnet die noetigen files (fehlerkontrolle)!                
!      open(unit=10,file='case.struct',form='formatted',status='old')
!      open (unit=20,file='case.plt',form='formatted',status='unknown')
      
      read(10,90) titel
      read(10,100) matrix,atome
      read(10,*)
      read(10,110) a,b,c,alpha,beta,gamma
      faktor=1
      if (matrix(1:1).eq.'F') faktor=4
      if (matrix(1:1).eq.'B') faktor=2
      if (matrix.eq.'CXY') faktor=2
      if (matrix.eq.'CYZ') faktor=2
      write (*,*)matrix
      
      if(gamma.eq.0.)then
      gamma=90.
      endif
      int=0
      do 10 i=1,atome
      int=int+1
      read(10,120) x(int),y(int),z(int) 
      read(10,130) mult(i)         
      
      if (matrix.eq.'CXY') then
        x(int+1)=x(int)+0.5d0
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)
        if (x(int+1).gt.1.d0) then
           x(int+1)=x(int+1)-1.d0
        endif
        if (y(int+1).gt.1.d0) then
           y(int+1)=y(int+1)-1.d0
        endif
        int=int+1
      endif
      
      if (matrix.eq.'CYZ') then
      write (*,*)matrix
        x(int+1)=x(int)
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)+0.5d0
        if (z(int+1).gt.1.d0) then
           z(int+1)=z(int+1)-1.d0
        endif
        if (y(int+1).gt.1.d0) then
           y(int+1)=y(int+1)-1.d0
        endif
        int=int+1
      endif
      
      if (matrix(1:1).eq.'B') then
        x(int+1)=x(int)+0.5d0
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)+0.5d0
        if (x(int+1).gt.1.d0) then
           x(int+1)=x(int+1)-1.d0
        endif
        if (y(int+1).gt.1.d0) then
           y(int+1)=y(int+1)-1.d0
        endif
        if (z(int+1).gt.1.d0) then
           z(int+1)=z(int+1)-1.d0
        endif
        int=int+1
      endif            
      
      if (matrix(1:1).eq.'F') then
        x(int+1)=x(int)+0.5d0
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)
        if (x(int+1).gt.1.d0) then
           x(int+1)=x(int+1)-1.d0
        endif
        if (y(int+1).gt.1.d0) then
           y(int+1)=y(int+1)-1.d0
        endif
        if (z(int+1).gt.1.d0) then
           z(int+1)=z(int+1)-1.d0
        endif
        x(int+2)=x(int)
        y(int+2)=y(int)+0.5d0
        z(int+2)=z(int)+0.5d0
        if (x(int+2).gt.1.d0) then
           x(int+2)=x(int+2)-1.d0
        endif
        if (y(int+2).gt.1.d0) then
           y(int+2)=y(int+2)-1.d0
        endif
        if (z(int+2).gt.1.d0) then
           z(int+2)=z(int+2)-1.d0
        endif
        x(int+3)=x(int)+0.5d0
        y(int+3)=y(int)
        z(int+3)=z(int)+0.5d0
        if (x(int+3).gt.1.d0) then
           x(int+3)=x(int+3)-1.d0
        endif
        if (y(int+3).gt.1.d0) then
           y(int+3)=y(int+3)-1.d0
        endif
        if (z(int+3).gt.1.d0) then
           z(int+3)=z(int+3)-1.d0
        endif
        int=int+3
      endif            
            
         
      if (mult(i).gt.1) then  
        do 30 ii=1,(mult(i)-1)
          int=int+1
          read(10,120) x(int),y(int),z(int)
  
      if (matrix.eq.'CXY') then
        x(int+1)=x(int)+0.5d0
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)
        if (x(int+1).gt.1) then
           x(int+1)=x(int+1)-1
        endif
        if (y(int+1).gt.1) then
           y(int+1)=y(int+1)-1
        endif
        int=int+1
      endif
      
      if (matrix.eq.'CYZ') then
        x(int+1)=x(int)
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)+0.5d0
        if (z(int+1).gt.1) then
           z(int+1)=z(int+1)-1
        endif
        if (y(int+1).gt.1) then
           y(int+1)=y(int+1)-1
        endif
        int=int+1
      endif
      
      if (matrix(1:1).eq.'B') then
        x(int+1)=x(int)+0.5d0
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)+0.5d0
        if (x(int+1).gt.1) then
           x(int+1)=x(int+1)-1
        endif
        if (y(int+1).gt.1) then
           y(int+1)=y(int+1)-1
        endif
        if (z(int+1).gt.1) then
           z(int+1)=z(int+1)-1
        endif
        int=int+1
      endif            
      
      if (matrix(1:1).eq.'F') then
        x(int+1)=x(int)+0.5d0
        y(int+1)=y(int)+0.5d0
        z(int+1)=z(int)
        if (x(int+1).gt.1) then
           x(int+1)=x(int+1)-1
        endif
        if (y(int+1).gt.1) then
           y(int+1)=y(int+1)-1
        endif
        if (z(int+1).gt.1) then
           z(int+1)=z(int+1)-1
        endif
        x(int+2)=x(int)
        y(int+2)=y(int)+0.5d0
        z(int+2)=z(int)+0.5d0
        if (x(int+2).gt.1) then
           x(int+2)=x(int+2)-1
        endif
        if (y(int+2).gt.1) then
           y(int+2)=y(int+2)-1
        endif
        if (z(int+2).gt.1) then
           z(int+2)=z(int+2)-1
        endif
        x(int+3)=x(int)+0.5d0
        y(int+3)=y(int)
        z(int+3)=z(int)+0.5d0
        if (x(int+3).gt.1) then
           x(int+3)=x(int+3)-1
        endif
        if (y(int+3).gt.1) then
           y(int+3)=y(int+3)-1
        endif
        if (z(int+3).gt.1) then
           z(int+3)=z(int+3)-1
        endif
        int=int+3
      endif            
  
  30    continue     
      goto 70 
      endif
  70    read(10,140) akernz
        do 40 j=int,int-(mult(i)*faktor-1),-1
        kernz(j)=nint(akernz)      
        rad(j)=a/10.d0
  
  40    continue    
      read(10,*)
      read(10,*)
      read(10,*)
  10  continue
  
      do 25 j=1,int
      
      if (x(j).le.1.) then
        x(j)=x(j)+1
          if (x(j).le.1.) then
          int=int+1
            x(int)=x(j)
            y(int)=y(j)
            z(int)=z(j)
            kernz(int)=kernz(j)
            rad(int)=a/10.d0
            if(y(j).le.1.) then
            y(j)=y(j)+1
              if(y(j).le.1.) then
              int=int+1
              x(int)=x(j)
              y(int)=y(j)
              z(int)=z(j)
              kernz(int)=kernz(j)  
              rad(int)=a/10.d0  
                if (z(j).le.1.) then
                z(j)=z(j)+1
                  if(z(j).le.1.) then
                  int=int+1
                  x(int)=x(j)
                  y(int)=y(j)
                  z(int)=z(j)
                  kernz(int)=kernz(j)
                  rad(int)=a/10.d0
                  endif
                  z(j)=z(j)-1
                endif   
              endif
            y(j)=y(j)-1
            endif
              if (z(j).le.1.) then
              z(j)=z(j)+1
                if(z(j).le.1.) then
                int=int+1
                x(int)=x(j)
                y(int)=y(j)
                z(int)=z(j)
                kernz(int)=kernz(j)
                rad(int)=a/10.d0
                endif
              z(j)=z(j)-1
              endif   
          endif
        x(j)=x(j)-1
      endif
      
      if (y(j).le.1.) then
        y(j)=y(j)+1
          if(y(j).le.1.) then
          int=int+1
          x(int)=x(j)
          y(int)=y(j)
          z(int)=z(j)
          kernz(int)=kernz(j)
          rad(int)=a/10.d0
          endif
        y(j)=y(j)-1
        endif   
          
      if (z(j).le.1.) then
        z(j)=z(j)+1
        if(z(j).le.1.) then
          int=int+1
          x(int)=x(j)
          y(int)=y(j)
          z(int)=z(j)
          kernz(int)=kernz(j)  
          rad(int)=a/10.d0
            if (y(j).le.1.) then
            y(j)=y(j)+1
              if(y(j).le.1.) then
                int=int+1
                x(int)=x(j)
                y(int)=y(j)
                z(int)=z(j)
                kernz(int)=kernz(j)
                rad(int)=a/10.d0
                endif
            y(j)=y(j)-1
            endif   
        endif
        z(j)=z(j)-1
      endif   
  
        if (x(j).ge.1.) then
        x(j)=x(j)-1
          if (x(j).le.1..and.x(j).ge.0.d0) then
          int=int+1
            x(int)=x(j)
            y(int)=y(j)
            z(int)=z(j)
            kernz(int)=kernz(j)
            rad(int)=a/10.d0
            if(y(j).ge.1.) then
            y(j)=y(j)-1
              if(y(j).le.1..and.y(j).ge.0.d0) then
              int=int+1
              x(int)=x(j)
              y(int)=y(j)
              z(int)=z(j)
              kernz(int)=kernz(j)  
              rad(int)=a/10.d0  
                if (z(j).ge.1.) then
                z(j)=z(j)-1
                  if(z(j).le.1..and.z(j).ge.0.d0) then
                  int=int+1
                  x(int)=x(j)
                  y(int)=y(j)
                  z(int)=z(j)
                  kernz(int)=kernz(j)
                  rad(int)=a/10.d0
                  endif
                  z(j)=z(j)+1
                endif   
              endif
            y(j)=y(j)+1
            endif
              if (z(j).ge.1.) then
              z(j)=z(j)-1
                if(z(j).le.1..and.z(j).ge.0.d0) then
                int=int+1
                x(int)=x(j)
                y(int)=y(j)
                z(int)=z(j)
                kernz(int)=kernz(j)
                rad(int)=a/10.d0
                endif
              z(j)=z(j)+1
              endif   
          endif
        x(j)=x(j)+1
      endif
      
      if (y(j).ge.1.) then
        y(j)=y(j)-1
          if(y(j).le.1..and.y(j).ge.0.d0) then
          int=int+1
          x(int)=x(j)
          y(int)=y(j)
          z(int)=z(j)
          kernz(int)=kernz(j)
          rad(int)=a/10.d0
          endif
        y(j)=y(j)+1
        endif   
          
      if (z(j).ge.1.) then
        z(j)=z(j)-1
        if(z(j).le.1..and.z(j).ge.0.d0) then
          int=int+1
          x(int)=x(j)
          y(int)=y(j)
          z(int)=z(j)
          kernz(int)=kernz(j)  
          rad(int)=a/10.d0
            if (y(j).ge.1.) then
            y(j)=y(j)-1
              if(y(j).le.1..and.y(j).ge.0.d0) then
                int=int+1
                x(int)=x(j)
                y(int)=y(j)
                z(int)=z(j)
                kernz(int)=kernz(j)
                rad(int)=a/10.d0
                endif
            y(j)=y(j)+1
            endif   
        endif
        z(j)=z(j)+1
      endif   
      

  25  continue    
  90  format(1x,a10)
 100  format(a3,24x,i3)     
   
 110  format(6f10.5)
 120  format(12x,f9.6,4x,f9.6,4x,f9.6)
 130  format(15x,i2)     
 140  format(56x,f4.1)     
      close(10,status='keep')
      
      call bmatr(gamma,matrix,int,a,b,c,x,y,z)
      
      icon=0
      do 350 i=1,int
      if(x(i).eq.0..and.y(i).eq.0..and.z(i).eq.0.) then
         icon=3
         do 360 ii=1,icon
         nuca(ii)=kernz(i)
         nuce(ii)=kernz(i)
         
360      continue
      rmin(1)=a-0.01
      rmax(1)=a+0.01
      rmin(2)=b-0.01
      rmax(2)=b+0.01
      rmin(3)=c-0.01
      rmax(3)=c+0.01
      endif
350   continue         
      dist0=9999.9
      do 500 j=1,int-1
        do 510 jj=j+1,int
          xx=x(j)-x(jj)
          yy=y(j)-y(jj)
          zz=z(j)-z(jj)
          dist=sqrt(xx*xx+yy*yy+zz*zz)
          if(dist.lt.dist0) then
             dist0=dist
             zsave1=kernz(j)
             zsave2=kernz(jj)
          endif
 510    continue
 500  continue            
      icon=icon+1
      
      ergeb=int
      
      write (20,190) titel
      write(20,200) ergeb,a,b,c
      do 20 j=1,int
      ikern(j)=idint(kernz(j))
      write(20,210) x(j),y(j),z(j),rad(j),ikern(j) 
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
      write(20,230) m,n,iclc,ibft,nshx,nshy,icon,initsv
      do 400 i=1,icon-1

      
      rbond(i)=a/100.d0
      nbcol(i)=4
      write(20,240) nuca(i),nuce(i),rmin(i),rmax(i),rbond(i), &
      nbcol(i)
  400 continue   
      rbon=a/50.d0
      nbco=5
      rmi=dist0-0.01
      rma=dist0+0.01
      write(20,240) zsave1,zsave2,rmi,rma,rbon,nbco
      close(20,status='keep')
      
 
 190  format(a60)
 200  format(i5,4f15.9)     
 210  format(4f15.9,i4)
 220  format(6f12.5)   
 230  format(8i5)     
 240  format(2i5,3d15.7,i5)
      stop      
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
      namstr(i+1:i+7)='.struct'
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
      
