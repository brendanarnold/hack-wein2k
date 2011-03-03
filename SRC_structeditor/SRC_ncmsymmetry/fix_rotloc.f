  subroutine fix_rotloc

    use struct
    use rotations
    use Ylm_rot

    integer i,ind
    real*8 rot(1:3,1:3) 
    logical s(48)                    
    real*8 rotold(3,3),rotnew(3,3)
    integer nps

    ind=1
    do i=1,nat
       rot(1,1)= cos(fi(ind))*cos(teta(ind))
       rot(2,1)= sin(fi(ind))*cos(teta(ind))
       rot(3,1)=             -sin(teta(ind))
       rot(1,2)=-sin(fi(ind))
       rot(2,2)= cos(fi(ind))
       rot(3,2)= 0.0d0
       rot(1,3)= cos(fi(ind))*sin(teta(ind))
       rot(2,3)= sin(fi(ind))*sin(teta(ind))
       rot(3,3)=              cos(teta(ind))          
       rotloc(1:3,1:3,i)=transpose(rot(1:3,1:3))
       ind=ind+mult(i)
    enddo

!!$    do i=1,nat       
!!$       call index_symop(i,s)
!!$       nps=npgopat(i)
!!$       rotold(1:3,1:3)=transpose(rotloc(1:3,1:3,i))
!!$       call class(nps,s,lattic,rotold,rotnew)
!!$       rotloc(1:3,1:3,i)=transpose(rotnew(1:3,1:3))
!!$    enddo

  end subroutine fix_rotloc

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine index_symop(i,s)

   use Ylm_rot, only : npgopat,pg_sym_oper
   use struct

   IMPLICIT REAL*8 (A-H,O-Z)

   logical s(48)
   integer i
   character*20 oname 
   real*8 rot(3,3), rott(3,3)
   logical  rotseq
   integer j,k
   real*8  x,y,z,x1,y1,z1

   s=.false.

   do j=1,npgopat(i)
      rot(1:3,1:3)=pg_sym_oper(i,j,1:3,1:3)
      do k=1,48

         x=1.0d0
         y=0.0d0
         z=0.0d0
         call symop(k,x,y,z,x1,y1,z1,oname,nz,lattic) 
         rott(1,1)=x1
         rott(2,1)=y1
         rott(3,1)=z1

         x=0.0d0
         y=1.0d0
         z=0.0d0
         call symop(k,x,y,z,x1,y1,z1,oname,nz,lattic) 
         rott(1,2)=x1
         rott(2,2)=y1
         rott(3,2)=z1

         x=0.0d0
         y=0.0d0
         z=1.0d0
         call symop(k,x,y,z,x1,y1,z1,oname,nz,lattic) 
         rott(1,3)=x1
         rott(2,3)=y1
         rott(3,3)=z1
         
         if (rotseq(rot,rott)) then
            s(k)=.true.
!!$            write(*,*) k
!!$            write(*,'(3f10.4)') rot
!!$            write(*,*)
            goto 1
         endif
      enddo    
1     continue
   enddo

 end subroutine index_symop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 logical function rotseq(rot,rott) 

   real*8 rot(3,3), rott(3,3)
   integer i,j
   real*8 zero
   logical ok

   ok=.true.
   zero=1.0d-8
   do i=1,3
      do j=1,3
         ok=ok.and.(abs(rot(i,j)-rott(i,j)).lt.zero)
      enddo
   enddo
   
   rotseq=ok

 end function rotseq
