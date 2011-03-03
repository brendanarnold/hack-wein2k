      subroutine inldau(ifr,ifw,nat,ieq)
!to this program serves to read input files in FLAPW WIEN program
! then to construct the input files 
! when s-o coupling and magnetization are present
       implicit real*8(a-h,o-z)
        dimension a(80),b(100,80)
        integer ieq(100)
        character*40 ch(2)
                 ch(1)=' LDA+U calculation, spin down'
                 ch(2)=' LDA+U calculation, spin up  '
        do j=1,3
        read(ifr,1980,end=12)(a(i),i=1,62)
        write(ifw,1980)(a(i),i=1,62)
        enddo
        do 80 jatom=1,nat
         read(ifr,*)l
         if(l.gt.0)read(ifr,1980)(a(i),i=1,71)
        do 81 jeq=1,ieq(jatom)
         write(ifw,100)l
         if(l.gt.0)write(ifw,1980)(a(i),i=1,71)
 81     continue
 80     continue
        read(ifr,1980)(a(i),i=1,62)
        write(ifw,1980)(a(i),i=1,62)
        do j=1,2
        ifr=ifr+1
        ifw=ifw+1
        jatso=0
        notU=0
        do jatom=1,nat
        read(ifr,*)jat
        if(jat.gt.0)then
        read(ifr,*)l
        l2=(2*l+1)**2
          do k=1,l2
          read(ifr,1980)(b(k,i),i=1,36)
          enddo
         endif
            do jeq=1,ieq(jatom)
            jatso=jatso+1
            if(jat.le.0)then
            write(ifw,103)notU
            else
            if(j.eq.1)then
            write(ifw,101)jatso
            else
            write(ifw,102)jatso
            endif
            write(ifw,104)l
              do k=1,l2
              write(ifw,1980)(b(k,i),i=1,36)
              enddo
            endif
            enddo
          enddo
        enddo
12      continue
 1980   format(80a1)
100     format(i3)
101       format(7x,i6,'  atom  LDA+U calculation, spin down')
102       format(7x,i6,'  atom  LDA+U calculation, spin up')
103       format(i4,'   LDA+U not considered')
104       format(i4,'  L for LDA+U calculation')
        return
        end
