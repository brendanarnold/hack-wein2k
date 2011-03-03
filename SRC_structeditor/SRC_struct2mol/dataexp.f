subroutine dataexp(logncm,add,logfor)

  USE struct
  USE atom_list   
  
   USE povprop

  logical       :: logncm,add,logfor
  INTEGER       :: i,j,ind
  CHARACTER     :: atom*10
  integer       :: sc(3),nchar
  real*8        :: mm,fi,teta,pi,zero,rgba_tmp(3)
  real*8        :: mmx(ndif),mmy(ndif),mmz(ndif)
  real*8        :: xyz(3),xyz2(3),fortmp(3),force(3,ndif),symr(3,3),t(3)
  real*8        :: tetai(ndif),fii(ndif),mmi(ndif),qq
  integer       :: is1,is2,is3,nsc,is,isb
  real*8, allocatable :: vdwr_sc(:),rgb_sc(:,:),list_rpos_sc(:,:),mm_sc(:,:)
  character*10, allocatable :: aname_sc(:)
  integer, allocatable :: list_iat_sc(:),list_ind_sc(:)
  real*8  tr(3)
  integer movea,moveb,movec,nlist_o
  character*50 number
  logical vectors_equal


   call init_povprop

   pi=2.0d0*asin(1.0d0)
   zero=0.0d0

   read(3,*)
   read(3,*) 
   read(3,*) 
   read(3,*) 
   read(3,*) 
   read(3,*)
   read(3,*) 
   read(3,*) 
   read(3,*) 
   read(3,*) 
   read(3,*)
   read(3,*) (sc(i),i=1,3)
   read(3,*)
   read(3,*) 
   read(3,*) ncol
   do i=1,ncol
      read(3,*) rgb_name(i),(rgb_tmp(k,i),k=1,3)  
   enddo
   read(3,*)
   read(3,*) nvdwr
   do i=1,nvdwr
      read(3,*) vdwr_name(i),vdwr_tmp(i) 
   enddo
   read(3,*)
   read(3,*) 
   read(3,*) 
   read(3,*) 
   read(3,*) nbond
   do i=1,nbond
      read(3,*) bond_name(i,1), bond_name(i,2), bond_dis(i) 
   enddo
   read(3,*)
   read(3,*) 
   read(3,*) imagcol, fmagcol
   read(3,*) (rgba_tmp(k),k=1,3)
   read(3,*) 
   read(3,*) 

   mmx(1:ndif)=0.0d0
   mmy(1:ndif)=0.0d0
   mmz(1:ndif)=0.0d0

   if (logncm) then
      mmm=0.0d0 
      ind=0
      read(4,*)      
      read(4,*)  (qspir(i),i=1,3) 
      DO i=1,nat
         read(5,*) 
         read(5,*) rnic,rnic,mm
         mm=mm*2.0d0
         do im=1,mult(i)
            read(4,*) fi,teta
            write(6,*) i,mm
            fi=fi*pi/180.0d0
            teta=teta*pi/180.0d0
            ind=ind+1
            tetai(ind)=teta
            fii(ind)=fi
            mmi(ind)=mm
            mmz(ind)=mm*cos(teta)
            mmx(ind)=mm*sin(teta)*cos(fi)
            mmy(ind)=mm*sin(teta)*sin(fi)
         enddo
      enddo
      read(4,*)
   endif

   if (logfor) then
      ind1=1
      do i=1,nat
         read(8,*) fortmp(1:3)
         force(1:3,ind1)=fortmp(1:3)
!write(0,'(i5,3f10.5,a)') ind1,force(1:3,ind1),'*'
         xyz(1:3)=pos(1:3,ind1)          
!         write(0,*) mult(i)
         do j=2,mult(i)
            ind2=ind1+j-1
!            write(0,*) ind1,ind2
            do isym=1,nsym
               symr(1:3,1:3)=iz(1:3,1:3,isym)
               t(1:3)=tau(1:3,isym)
               symr=transpose(symr)
               xyz2(1:3)=pos(1:3,ind2) 
               xyz2=matmul(symr,xyz2)
               xyz2=xyz2+t+10.0d0
               xyz2(1:3)=mod(xyz2(1:3),1.0d0)
!               write(0,'(3f10.5)') xyz
!               write(0,'(3f10.5)') xyz2
!                write(0,'(3f10.5)') symr
!                write(0,'(3f10.5,/)') t

               if (vectors_equal(xyz,xyz2,5.0d-4)) then
                  force(1:3,ind2)=matmul(symr,fortmp)
!                  write(0,'(i5,3f10.5)') ind2,force(1:3,ind2)
                  goto 1
               endif
            enddo
            write(0,*) 'sym operation not found'
            write(0,*) i,j
            stop
1           continue  
         enddo
         ind1=ind1+mult(i)
      enddo
   endif

   nbond_in=nbond

   call generate_rgb
   call generate_vdwr
   call generate_bonds

   if (nbond.eq.0) then
      nbond=1
      bond(1,1)=1
      bond(1,2)=1
   endif

   nsc=sc(1)*sc(2)*sc(3)*nlist
   allocate (vdwr_sc(nsc),rgb_sc(3,nsc),list_rpos_sc(3,nsc),mm_sc(3,nsc),aname_sc(nsc))
   allocate (list_iat_sc(nsc),list_ind_sc(nsc))

   mm_sc=0.0d0
   nlist_o=nlist

   is=0
   isb=0
   nsc=0
   do is1=0,sc(1)-1      
      do is2=0,sc(2)-1      
         do is3=0,sc(3)-1
            nsc=nsc+1
            do i=1,nlist
               if (add) then
                  if ((is1.gt.1).and.(abs(list_pos(1,i)).lt.zero)) cycle
                  if ((is2.gt.1).and.(abs(list_pos(2,i)).lt.zero)) cycle
                  if ((is3.gt.1).and.(abs(list_pos(3,i)).lt.zero)) cycle
               endif
               is=is+1
               list_iat_sc(is)=list_iat(i)
               list_ind_sc(is)=list_ind(i)

               vdwr_sc(is)=vdwr(list_iat(i))
               rgb_sc(1:3,is)=rgb(1:3,list_iat(i))

               tr(1)=list_rpos(1,i)+is1*br(1,1)+is2*br(2,1)+is3*br(3,1)
               tr(2)=list_rpos(2,i)+is1*br(1,2)+is2*br(2,2)+is3*br(3,2)
               tr(3)=list_rpos(3,i)+is1*br(1,3)+is2*br(2,3)+is3*br(3,3)
               list_rpos_sc(1:3,is)=tr(1:3)


               if (nsc.eq.1) then
                  write(number,*) list_ind(i)
                  number=adjustl(number)
                  aname_sc(is)=trim(aname(list_iat(i)))//'_'//number
               else
                  aname_sc(is)=' '
               endif
               if (logncm) then
                  movea=nint(list_pos(1,i)-pos(1,list_ind(i)))
                  moveb=nint(list_pos(2,i)-pos(2,list_ind(i)))
                  movec=nint(list_pos(3,i)-pos(3,list_ind(i)))                 
                  qq=2.0d0*pi*(qspir(1)*(movea+is1)+qspir(2)*(moveb+is2)+qspir(3)*(movec+is3))
                  qq=qq*180.0d0/pi
                  mm_sc(3,is)=mmi(list_ind(i))*cos(tetai(list_ind(i)))
                  mm_sc(1,is)=mmi(list_ind(i))*sin(tetai(list_ind(i)))*cos(fii(list_ind(i))+qq)
                  mm_sc(2,is)=mmi(list_ind(i))*sin(tetai(list_ind(i)))*sin(fii(list_ind(i))+qq)
               endif
               if (logfor) then
                  mm_sc(1:3,is)=force(1:3,list_ind(i))
               endif
            enddo
         enddo
      enddo
   enddo
   br(1,1:3)= br(1,1:3)*sc(1)
   br(2,1:3)= br(2,1:3)*sc(2)
   br(3,1:3)= br(3,1:3)*sc(3)

   deallocate(list_iat,list_rpos,bond)
   maxnb=nsc*nbond*100
   allocate(list_iat(nsc*nlist),list_rpos(3,nsc*nlist),bond(maxnb,2))
   nlist=nsc*nlist
   list_iat=list_iat_sc
   list_rpos=list_rpos_sc
   nbond=nbond_in


   call generate_bonds

   if (nbond.eq.0) then
      nbond=1
      bond(1,1)=1
      bond(1,2)=1
   endif

   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  atomic radii'   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') 'object   1'   
   write(2,'(a,i5)') 'class array     type float     rank 0     items',is
   write(2,'(a)') 'data follows'
   write(2,'(10f15.8)') (vdwr_sc(i),i=1,is)
   write(2,'(a)') 'attribute "dep" string "positions"'
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  atomic colors'   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') 'object   2'   
   write(2,'(a,i5)') 'class array type float rank 1 shape 3 items',is
   write(2,'(a)') 'data follows'
   write(2,'(10f15.8)') ((rgb_sc(j,i),j=1,3),i=1,is)
   write(2,'(a)') 'attribute "dep" string "positions"'
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  atomic positions'   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') 'object   3'   
   write(2,'(a,i5)') 'class array type float rank 1 shape 3 items',is
   write(2,'(a)') 'data follows'
   write(2,'(10f15.8)') ((list_rpos_sc(j,i),j=1,3),i=1,is)
   write(2,'(a)') 'attribute "dep" string "positions"'
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  atomic momenta'   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') 'object   4'   
   write(2,'(a,i5)') 'class array type float rank 1 shape 3 items',is
   write(2,'(a)') 'data follows'
   write(2,'(10f15.8)') (mm_sc(1:3,i),i=1,is)
   write(2,'(a)') 'attribute "dep" string "positions"'
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  bonds'   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') 'object   5'   
   write(2,'(a,i5)') 'class array type int rank 1 shape 2 items',nbond
   write(2,'(a)') 'data follows'
   write(2,'(10i5)') (bond(i,1)-1,bond(i,2)-1,i=1,nbond)
   write(2,'(a)') 'attribute "ref" string "positions"'
   write(2,'(a)') 'attribute "element type" string "lines"'
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  unit cell'   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') 'object   6'   
   write(2,'(a)') 'class array type float rank 1 shape 3 items 8'
   write(2,'(a)') 'data follows'
   write(2,'(3f15.8)') zero,zero,zero
   write(2,'(3f15.8)') (br(1,i),i=1,3)  
   write(2,'(3f15.8)') (br(2,i),i=1,3)  
   write(2,'(3f15.8)') (br(1,i)+br(2,i),i=1,3)  
   write(2,'(3f15.8)') (br(3,i),i=1,3)  
   write(2,'(3f15.8)') (br(1,i)+br(3,i),i=1,3)  
   write(2,'(3f15.8)') (br(2,i)+br(3,i),i=1,3)  
   write(2,'(3f15.8)') (br(1,i)+br(2,i)+br(3,i),i=1,3)  
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  atom names'   
   write(2,'(a)') '#  '   
   write(2,'(a)') '#  '   

   call find_longest_name(aname_sc,is,nchar)

   write(2,'(a)') 'object   7'   
   write(2,'(a,i5,a,i5)') 'class array type string rank 1 shape',nchar+1,'   items',is
   write(2,'(a)') 'data follows'
   write(2,'(10(3a))') ('"',aname_sc(i)(1:nchar),'"  ',i=1,is)
   write(2,'(a)') 'attribute "dep" string "positions"'
   write(2,'(a)') '#  '   
   write(2,'(a)') 'object "ballstick"'
   write(2,'(a)') 'class field'
   write(2,'(a)') 'component "data" value          1'
   write(2,'(a)') 'component "colors" value         2'
   write(2,'(a)') 'component "positions" value          3'
   write(2,'(a)') 'component "connections" value          5'
   write(2,'(a)') 'component "box" value          6'
   write(2,'(a)') 'attribute "name" string "atoms"'
   write(2,'(a)') '#'
   write(2,'(a)') 'object "ballstick_ncm"'
   write(2,'(a)') 'class field'
   write(2,'(a)') 'component "positions" value          3'
   write(2,'(a)') 'component "data" value          4'
   write(2,'(a)') 'attribute "name" string "atoms"'
   write(2,'(a)') '#'
   write(2,'(a)') 'object "names"'
   write(2,'(a)') 'class field'
   write(2,'(a)') 'component "positions" value          3'
   write(2,'(a)') 'component "data" value          7'
   write(2,'(a)') 'attribute "name" string "atoms"'
   write(2,'(a)') '#'
   write(2,'(a)') ' end'

end subroutine dataexp

 subroutine find_longest_name(aname,nat,nchar)

 character*10 aname(nat)
 integer i,j,nat,nchar

 nchar=0
 do i=1,nat
    do j=10,1,-1
       if (aname(i)(j:j).ne.' ') then
          nchar=max(nchar,j)
          go to 1
       endif
    enddo
1   continue
 enddo
  
 end subroutine find_longest_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function vectors_equal(xyz,xyz2,tol)

    real*8 xyz(3), xyz2(3),tol
    real*8 dxyz(3)
    integer k

    dxyz(1:3)=abs(xyz(1:3)-xyz2(1:3))
    do k=1,3
       if (abs(mod(dxyz(k),1.0d0)).lt.tol) dxyz(k)=0.0d0
!       if (abs(dxyz(k)-1.0d0).lt.tol) dxyz(k)=0.0d0
    enddo
    vectors_equal=dxyz(1).lt.tol.and.dxyz(2).lt.tol.and.dxyz(3).lt.tol

  end function vectors_equal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
