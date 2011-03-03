      program strlat
!     felddimensionierung, hoechstzahl behandelbarer atome
      implicit real*8 (a-h,o-z)
      parameter (nato=100)
      character*4 typ
      dimension x(4,nato),klz(nato)
!     x.....Atomkoordinaten, Kugelradius
!     klz...Kernladungszahl      
      call o_file
!     oeffnet die noetigen files (fehlerkontrolle)!                
      call latsta
!     schreibt standarddaten in lat-file. erweiterbar zur nutzung von
!     zusaetzlichen optionen
      call matrix(nat,typ,a,b,c)
!     schreibt bravais-matrix, uebergibt nat (number of atoms)
      call koord(nat,nato,na,x,klz)
      call write(na,x,klz,typ,a,b,c)
      call option
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
      namlat(i+1:i+4)='.lat'
      open(unit=1,file=namstr,form='formatted',status='old',err=9993)
      open(unit=2,file=namlat,form='formatted',status='new',err=9994)
      return
9992  write(*,*)'Filename too long' 
      stop      
9993  write(*,*)'struct-File does not exist'      
      stop
9994  write(*,*)'File ',namlat(1:i+4),' already exists.'
      write(*,*)' Overwrite (y/n)?'
      read(*,99901)test
      if(test.eq.'y' .or. test.eq.'Y') then
       open(unit=2,file=namlat,form='formatted')
       return
      else
       stop
      endif    
99900 format(a80) 
99901 format(a1)     
      end

      subroutine latsta
      write(2,99900)50,0,' ',1,12
      write(2,99901)10,1.0,'wurzite'
      return
99900 format(i6,i2,a20,25i2)
99901 format(i5,f15.9,a)
      end
      
      subroutine matrix(nat,typ,a,b,c)
      implicit real*8 (a-h,o-z)      
      dimension bra(3,3)
      character*4 typ
      do 10 i=1,3
       do 10 j=1,3
10      bra(j,i)=0d0
      read(1,99900)
      read(1,99901) typ,nat
      read(1,99902)
      read(1,99903) a,b,c,alpha,beta,gamma
      if(gamma.eq.0d0) gamma=90d0
      gamma=gamma*acos(0d0)/90d0
!     umrechnung ins bogenmass      
      if(typ(1:3).eq.'CXY') then
       bra(1,1)=a/2d0
       bra(2,1)=-b/2d0
       bra(1,2)=a/2d0
       bra(2,2)=b/2d0
       bra(3,3)=c       
      else if (typ(1:3).eq.'CYZ') then
       bra(1,1)=a
       bra(2,2)=-b/2d0
       bra(3,2)=c/2
       bra(2,3)=b/2d0
       bra(3,3)=c/2d0
      else if(typ(1:3).eq.'CXZ') then
       bra(1,1)=a*sin(gamma)/2d0
       bra(2,1)=a*cos(gamma)/2d0
       bra(3,1)=-c/2d0
       bra(2,2)=b
       bra(1,3)=a*sin(gamma)/2d0
       bra(2,3)=a*cos(gamma)/2d0
       bra(3,3)=c/2d0
      else if(typ(1:1).eq.'P') then
       bra(1,1)=a*sin(gamma)
       bra(2,1)=a*cos(gamma)
       bra(2,2)=b
       bra(3,3)=c
      else if(typ(1:1).eq.'F') then
       bra(1,1)=a/2d0
       bra(2,1)=b/2d0
       bra(1,2)=a/2d0
       bra(3,2)=c/2d0
       bra(2,3)=b/2d0
       bra(3,3)=c/2d0
      else if(typ(1:1).eq.'B') then
       bra(1,1)=a/2d0
       bra(2,1)=-b/2d0
       bra(3,1)=c/2d0
       bra(1,2)=a/2d0
       bra(2,2)=b/2d0
       bra(3,2)=-c/2d0
       bra(1,3)=-a/2d0
       bra(2,3)=b/2d0
       bra(3,3)=c/2
      else if(typ(1:1).eq.'R') then
       bra(1,1)=a/sqrt(3d0)/2d0
       bra(2,1)=-a/2d0
       bra(3,1)=c/3d0
       bra(1,2)=a/sqrt(3d0)/2d0
       bra(2,2)=a/2d0
       bra(3,2)=c/3d0
       bra(1,3)=-a/sqrt(3d0)
       bra(3,3)=c/3d0
      else if(typ(1:1).eq.'H') then
       bra(1,1)=a*sqrt(3d0)/2d0
       bra(2,1)=-a/2d0
       bra(2,2)=a
       bra(3,3)=c
      else 
       goto 9990
      end if
      write(2,99904)bra
      return
9990  write(*,*)'Unknown lattice type'
      stop      
99900 format(a80)      
99901 format(a4,23x,i3)      
99902 format(13x,a4)
99903 format(6F10.7)
99904 format(2(3F15.9/),3f15.9)
      end      
      subroutine koord(nat,nato,na,x,klz)
      implicit real*8 (a-h,o-z)
      dimension x(4,*),klz(*)
      na=0
      do 10 i=1,nat
       na=na+1
       if(na.gt.nato) goto 9990
       read(1,99990)x(1,na),x(2,na),x(3,na)
       read(1,99991)mult
       do 11 im=2,mult
        na=na+1
        if(na.gt.nato) goto 9990
11      read(1,99990)x(1,na),x(2,na),x(3,na)        
       read(1,99992)rmt,z
       do 12 im=na-mult+1,na
        x(4,im)=rmt/3d0
12      klz(im)=z
       read(1,99993)
       read(1,99993)
10     read(1,99993)
      return      
9990  write(*,*)'Too many atoms. Increase parameter nato in source ', &
      'code.'
      stop
99990 format(5x,3x,4x,f10.7,3x,f10.7,3x,f10.7)
99991 format(15x,i2,17x,i2)
99992 format(10x,5x,5x,5x,10x,5x,f10.5,5x,f5.2)
99993 format(20x,3f10.8)
      end

      subroutine write(na,x,klz,typ,a,b,c)
      real*8 x,a,b,c
      character*4 typ
      dimension x(4,*),klz(*)
      if((typ(1:1).eq.'F').or.(typ(1:1).eq.'B')) then
       write(2,99900)na,0
      else
       write(2,99900)-na,0 
      end if
      if((typ(1:1).eq.'F').or.(typ(1:1).eq.'B')) then
         do 10 i=1,na
10       write(2,99901)x(1,i)*a,x(2,i)*b,x(3,i)*c,x(4,i),klz(i)
      else
         do 11 i=1,na
11       write(2,99901)x(1,i),x(2,i),x(3,i),x(4,i),klz(i)
      end if
      return
99900 format(2i5)
99901 format(4f15.9,i10)
      end
      
      subroutine option
      implicit real*8 (a-h,o-z)
      write(2,99900)0,0,1,0
!     Millersche Indizes der zu waehlenden Gitterebene
      write(2,99901)4,4,5,1,0d0,0d0,0d0,0d0,0d0      
!     
      nrclx=0
      if(nrclx.gt.0) goto 1000 
!     sonst sind weitere Zeilen einzufuegen
      write(2,99902)nrclx
!     nrclx: Anzahl der zu restrukturierenden Ebenen (<12)
      write(2,99905)80d0,-30d0,0d0,1d0,11.25d0,-0.5d0
!     zwei Winkel zur Angabe der Blickrichtung
!     
      icon=0
      initsv=0
      if(icon.ne.0) goto 1000
      if(initsv.ne.0) goto 1000
!     sonst sind weitere Zeilen einzufuegen      
      write(2,99906)19,1,1,1,0,0,icon,initsv      
!     Grafik-Option, Farbcode, Code zur Atombenennung, ,
!     horizontaler, vertikaler Shift des Plotes in Pixeleinheiten, ,
!     Initialisierungsflag
      close(2)
      write(*,*)'Translation has been finished.'
      return
1000  write(*,*)'Option not aviable'
      stop      
99900 format(4i5)      
99901 format(4i5,5f10.3)
99902 format(12i5)
99903 format(i5,7f10.6)
99904 format(4f10.6,i5)
99905 format(6f12.5)
99906 format(8i5)
99907 format(2i5,3d15.7,i5)
      end      

