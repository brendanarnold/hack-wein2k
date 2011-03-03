subroutine xyz

  USE struct
  USE atom_list   
  
  INTEGER       :: i,j
  CHARACTER     :: atom*10
  
  WRITE(2,*) nlist
  WRITE(2,'(a)') title     
  
  do i=1,nlist
     iat=list_iat(i)
     atom=aname(iat)
     WRITE(2,'(a2,3f15.8)') &
          atom(1:3),(list_rpos(j,i),j=1,3)
  enddo

end subroutine xyz
