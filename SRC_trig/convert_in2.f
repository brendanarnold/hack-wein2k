      character*80 a
      character*196 b
      read(*,'(a)') a
      write(*,'(a)') a
      read(*,'(a)') a
      write(*,'(a)') a
      read(*,'(a)') a
      write(*,'(a)') a
      ifirst=0
 1     read(*,'(a)',end=10,err=10) b
      if (ifirst.eq.0) then
        if(b(2:2).ne.'0') then
         write(*,*) 'WIEN2k_03'
         goto 10
        endif
        ifirst=1
      endif  
      do i=196,1,-1
      if(b(i:i).eq.'X') goto 2
      if(b(i:i).ne.' ') exit
      enddo
      if(b(2:2).eq.'0') then
          write(*,'(196(1x,a2,a2))') (b(j:j+1),b(j+2:j+3),j=1,i-3,4)
      else
          write(*,'(196(a3,a2))') (b(j:j+2),b(j+3:j+4),j=1,i-4,5)
      endif
      goto 1
2     write(*,'(a)') b(1:80)
      read(*,'(a)') a
      write(*,'(a)') a
10    continue
      end

 
