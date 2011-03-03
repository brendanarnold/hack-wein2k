      subroutine in1ch(ifr,ifw,nat,ieq)
! this program serves to read input files in FLAPW WIEN program
! then to construct the input files 
! when s-o coupling and magnetization are present
       implicit real*8(a-h,o-z)
       character a,b
       character*79 aa
        dimension a(79),b(100,60)
        integer ieq(*)
        read(ifr,1980)(a(i),i=1,27)
        write(ifw,1980)(a(i),i=1,27)
        read(ifr,1980)(a(i),i=1,52)
        write(ifw,1980)(a(i),i=1,52)
        do 80 jatom=1,nat
        read(ifr,1981) aa
        do i=1,79
        if(aa(i:i).ne.' ') then
           do j=i+1,79
           if(aa(j:j).eq.' ') then
             do k=j+1,79
             if(aa(k:k).ne.' ') then
              a(1)=aa(k:k)
              read(a,'(i1)') n
              goto 100
             endif
             enddo
           endif
           enddo
        endif
        enddo
 100    continue
         do j=1,n
         read(ifr,1980)(b(j,i),i=1,29)
         enddo
        do 81 jeq=1,ieq(jatom)
        write(ifw,1981) aa
         do j=1,n
         write(ifw,1980)(b(j,i),i=1,29)
         enddo
 81     continue
 80     continue
 82     read(ifr,1980,end=888,err=888)(a(i),i=1,79)
        write(ifw,1980)(a(i),i=1,79)
        goto 82
 888    continue
 1980   format(80a1)
 1981   format(a79)
        return
        end
