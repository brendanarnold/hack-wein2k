subroutine xtl
  
  USE struct
  USE atom_list   
  
  INTEGER       :: i,j,ind
  CHARACTER     :: atom*10
  
  WRITE(2,'(2a)') 'TITLE     ',title     
  WRITE(2,'(a)') 'DIMENSION 3'     
  WRITE(2,'(a)') 'CELL'     
  WRITE(2,'(6(f10.5,x))') aa,bb,cc,alpha(1),alpha(2),alpha(3)
  WRITE(2,'(a)') 'SYMMETRY NUMBER  1'
  WRITE(2,'(a)') 'SYMMETRY LABEL P1'     
  WRITE(2,'(a)') 'ATOMS'     
  WRITE(2,'(a)') 'NAME         X           Y           Z'
  
  do i=1,nlist
     iat=list_iat(i)
     atom=aname(iat)
     WRITE(2,'(a,4x,3f12.5)') &
          atom(1:3),(list_pos(j,i),j=1,3)
  enddo
  
  WRITE(2,'(a)') 'EOF'     
  
end subroutine xtl
