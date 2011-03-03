 subroutine addatom2list(pos,iat,ind,movex,movey,movez,molf)
   
   USE atom_list
   USE struct, only: aname
   
   REAL*8        :: pos(3),movex,movey,movez
   integer       :: iat,ind
   CHARACTER     :: molf*3
   real*8        :: zero

   zero=1.0d-5
   nlist=nlist+1
   
   if (nlist.gt.natlist) then
      write(6,*) 'nlist.gt.natlist'
      stop
   endif
   
   list_pos(1,nlist)=pos(1)+movex
   list_pos(2,nlist)=pos(2)+movey
   list_pos(3,nlist)=pos(3)+movez

   if (molf.ne.'pov') then
   if (molf.ne.'dx ') then
      do k=1,3
         if (list_pos(k,nlist).lt.zero) list_pos(k,nlist)=list_pos(k,nlist)+0.001 
         if (abs(list_pos(k,nlist)-1.0d0).lt.zero) list_pos(k,nlist)=list_pos(k,nlist)-0.001 
      enddo
   endif
   endif

   do i=1,3
      list_rpos(i,nlist)=0.0d0
      do j=1,3
         list_rpos(i,nlist)=list_rpos(i,nlist)+br(j,i)*list_pos(j,nlist)
      enddo
   enddo

   list_iat(nlist)=iat
   list_ind(nlist)=ind
   

   write(6,'(a,a,3i5)') 'atom added to list:',aname(iat),nlist,iat,ind
   write(6,*)  (list_pos(i,nlist),i=1,3)   
   write(6,*)  (list_rpos(i,nlist),i=1,3)   

 end subroutine addatom2list
