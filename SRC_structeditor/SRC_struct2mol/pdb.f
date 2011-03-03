subroutine pdb

   USE struct
   USE atom_list   

   INTEGER       :: i,j,ind
   CHARACTER     :: atom*10

      WRITE(2,'(2a)') 'TITLE     ',title     
      WRITE(2,'(a,3f9.3,3f7.2,a)') 'CRYST1', &
           aa,bb,cc,alpha(1),alpha(2),alpha(3),' P1'
   
      do i=1,nlist
         iat=list_iat(i)
         atom=aname(iat)
         WRITE(2,'(a,i5,2x,a,14x,3f8.3)') &
              'ATOM  ',i,atom(1:3),(list_rpos(j,i),j=1,3)
      enddo

    end subroutine pdb
