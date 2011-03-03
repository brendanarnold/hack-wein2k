  subroutine make_point_groups

    use struct
    use Ylm_rot, only: pg_sym_oper,pgtinv,npgopat

    real*8, allocatable  :: ct(:,:)
    real*8   symr(3,3),t(3)
    integer i,j,ind0,ind1,ind2,nct,ict
    real*8  tol,xyz(3),xyz2(3)
    logical vectors_equal

    tol=1.0d-4

    if (lattic(1:1).eq.'P') then
       nct=1
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
    else if (lattic(1:1).eq.'F') then
       nct=4
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.5d0,0.0d0/)
       ct(1:3,3)=(/0.5d0,0.0d0,0.5d0/)
       ct(1:3,4)=(/0.0d0,0.5d0,0.5d0/)
    else if (lattic(1:1).eq.'B') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.5d0,0.5d0/)
    else if (lattic(1:1).eq.'H') then
       nct=1
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
    else if (lattic(1:1).eq.'R') then
       nct=1
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
    else if (lattic(1:3).eq.'CXY') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.5d0,0.0d0/)
    else if (lattic(1:3).eq.'CYZ') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.0d0,0.5d0,0.5d0/)
    else if (lattic(1:3).eq.'CXZ') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.0d0,0.5d0/)
    endif

    ind0=0

    write(6,'(//,a)')'-----------------------------'    
    write(6,'(a)')'   point groups for atoms'
    write(6,'(a,/)')'-----------------------------'    

    do i=1,nat

       write(6,*) 'atom=',i

       npgopat(i)=0

       if (i.gt.1) ind0=ind0+mult(i-1)

       ind1=ind0+1
       ind2=ind0+1

       do j=1,nsym

         symr(1:3,1:3)=iz(1:3,1:3,j)
         t(1:3)=tau(1:3,j)
         symr=transpose(symr)

         xyz(1:3)=pos(1:3,ind1)           
         xyz=matmul(symr,xyz)
         xyz=xyz+t+1.0d0
         xyz(1:3)=mod(xyz(1:3),1.0d0)
         
         do ict=1,nct
            xyz2(1:3)=(pos(1:3,ind2)+ct(1:3,ict))
            if (vectors_equal(xyz,xyz2,tol)) then
               ! point group operation
               npgopat(i)=npgopat(i)+1
               pg_sym_oper(i,npgopat(i),1:3,1:3)=symr(1:3,1:3)
               pgtinv(i,npgopat(i))=tinv(j)
               write(6,*) 'operation:',npgopat(i)
               write(6,'(10x,3i5)') int(transpose(pg_sym_oper(i,npgopat(i),1:3,1:3)))
               write(6,'(10x,3f10.5)') t
               write(6,'(10x,a,i3)') 'time inverion:',int(pgtinv(i,npgopat(i)))
               goto 1
            endif
         enddo    
1        continue

       enddo
    enddo

    do i=1,nat
       write(6,'(a,2i5)') 'number of operations for atom',i,npgopat(i)
       if (npgopat(i).le.0) then
          write(6,*)' error in make_point groups: npgopat(i).le.0'
          stop ' error in make_point groups: npgopat(i).le.0'
       endif
    enddo


  deallocate(ct)

  end subroutine make_point_groups


