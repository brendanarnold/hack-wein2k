 subroutine scan_octave(filein,nato)

   use struct, only : tobohr,lattyp,sgnum,HRtransf,hex2rho,rho2hex

   integer       nato
   character*80  filein,line
   real*8        orig(3),x(3,nato)
   real*8        a,b,c,alfa,beta,gamma
   real*8        aa,bb,cc
   character*50  sgname,aname(nato)
   integer       n,nat
   real*8        pi

      pi=acos(-1.0d0)
      tobohr=1.0d0
      orig=0.0d0

      open(unit=1,file=filein,status='old')
      call jumpto('a')
      read(1,*) a,b,c
      call jumpto('alpha')
      read(1,*) alfa,beta,gamma
      call jumpto('spgn')
      read(1,*)
      read(1,'(a)') sgname
      call jumpto('nat')
      read(1,*) nat
      call jumpto('pos')
      do i=1,nat
         read(1,*) x(1:3,i)
      enddo   
      call jumpto('aname')
      do i=1,nat      
         read(1,*)
         read(1,'(a)') aname(i)
      enddo
      close (1)

!      nat=1

     call test_sgname(sgname,.false.)

     if (sgnum.eq.167.or.sgnum.eq.166.or.sgnum.eq.161.or.sgnum.eq.160.or.&
         sgnum.eq.155.or.sgnum.eq.148.or.sgnum.eq.146) then
        HRtransf=.true.
        small=1.0d-6
        if ((abs(alfa-90.0d0).lt.small).and.&  
            (abs(beta-90.0d0).lt.small).and.&  
            (abs(gamma-120.0d0).lt.small)) then
!           do nothing
        else
           if (lattyp.ne.'R  ') then
              write(*,'(2a)') 'lattyp should be R and is: ',trim(lattyp)
              stop
           endif   
           alfa=alfa*pi/180.0d0
           aa=a*2.0d0*cos((pi-alfa)/2.0d0)
           bb=aa
           cc=3.0d0*sqrt(a**2-(aa**2)/3.0d0)
           a=aa
           b=bb
           c=cc
           alfa=90.0d0
           beta=90.0d0
           gamma=120.0d0
           do i=1,nat
              x(1:3,i)=matmul(x(1:3,i),rho2hex)
              do j=1,3
                 if (x(j,i).lt.0.0d0) x(j,i)=x(j,i)+1.0d0
                 if (x(j,i).gt.1.0d0) x(j,i)=x(j,i)-1.0d0
              enddo
           enddo

        endif
     endif

        open(unit=1,status='scratch')        
!     open(unit=1,file='rolask')        
     write(1,'(3f10.6)') orig(1:3)
     write(1,'(6f12.6)') a,b,c,alfa,beta,gamma
     sgname=''''//trim(sgname)//''''
     write(1,'(a)') sgname
     do i=1,nat
        call insert_string(aname(i),'     ')
        aname(i)=''''//trim(aname(i))//''''         
        write(1,'(a)') aname(i)
        write(1,'(3f16.7)') (x(j,i),j=1,3)  
     enddo
     rewind 1


   end subroutine scan_octave

    subroutine jumpto(var)

       character var*(*)
       integer iunit

       character h*1,line*100,var1*50 

       iunit=1
       rewind(iunit)

1      continue

       read(iunit,'(a)',end=2) line 
!           write(0,*) line
       if (line(1:7).eq.'# name:') then
          var1=adjustl(line(8:100))
!          write(*,'(a,a)') 'var1=',var1
!          write(*,'(a,a)') ' var=',var
          if (trim(var1).eq.trim(var)) then
3            continue
             read(iunit,'(a)') line
!             write(*,'(a)') line(1:7) 
             if (line(1:7).eq.'# type:') then
                var1=adjustl(line(8:100))
!                write(*,'(a,a)') 'var1=',var1 
                if (trim(var1).eq.'cell') then
                   read(iunit,*)
                   read(iunit,*)
                   read(iunit,*)
                   goto 3
                elseif (trim(var1).eq.'matrix') then
                   read(iunit,*)
                   read(iunit,*)
                   goto 2
                elseif (trim(var1).eq.'scalar') then
                   goto 2
                elseif (trim(var1).eq.'string') then
                   read(iunit,*)
                   goto 2
                endif
             endif 
          endif
       endif

       goto 1

2     continue

     end subroutine jumpto
