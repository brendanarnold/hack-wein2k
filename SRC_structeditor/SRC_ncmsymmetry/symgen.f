 subroutine symgen
   
    use struct
    use rotations, only: sspiral

    real*8  rot(3,3),t(3),tt(3),a(3),brinv(3,3)
    real*8, allocatable::  tsymo(:,:), tsymos(:,:)
    integer, allocatable::  symo(:,:,:),symos(:,:,:)
    integer is,k,i,j,iat,nt,kk
    real*8  small
    real*8, allocatable ::  ts(:,:)
    real*8  x,y,z,x1,y1,z1
    character*20 oname 
    integer nz,mnsym,j1,j2,iat1,iat2,iat0
    real*8  xyz1(3),xyz2(3),sym1(3,3),sym2(3,3)
    real*8, allocatable  :: ct(:,:)
    integer nct   
    real*8 SYM(4,3,48)
    integer nsymp
    logical aone


    call inversa(br2_dir,brinv)
    CALL PGLSYM(br2_dir,sym,nsymp)

    mnsym=(nsymp+1)*nsymp*ndif*ndif
    allocate (symo(3,3,nsymp),tsymo(3,nsymp),ts(3,mnsym))
    allocate (symos(3,3,nsymp),tsymos(3,nsymp))

    small=1.0d-6
    symo=0.0d0
    tsymo=0.0d0

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


    is=1
    symo(1:3,1,1)=(/1,0,0/)
    symo(1:3,2,1)=(/0,1,0/)
    symo(1:3,3,1)=(/0,0,1/)
    tsymo(1:3,1)= (/0,0,0/)
    write(6,*) is
    write(6,'(3i5)') symo(1:3,1:3,is) 

    if (sspiral) goto 10

      do i=1,nsymp
         rot=sym(1:3,1:3,i)
         rot=transpose(rot)
         if(.not.ortho.and.lattic(1:3).ne.'CXZ') then
            rot=matmul(brinv,rot)
            rot=matmul(rot,br2_dir)
         endif
         if (.not.aone(rot)) then
            is=is+1
            symo(1:3,1:3,is)=nint(transpose(rot(1:3,1:3)))
            tsymo(1:3,i)=(/0.0d0,0.0d0,0.0d0/)
            write(6,*) is
            write(6,'(3i5)') symo(1:3,1:3,is) 
          endif
      end do

!!$    do k=2,48
!!$       
!!$       x=1.0d0
!!$       y=0.0d0
!!$       z=0.0d0
!!$       call symop(k,x,y,z,x1,y1,z1,oname,nz,lattic) 
!!$       rot(1,1)=x1
!!$       rot(2,1)=y1
!!$       rot(3,1)=z1
!!$       
!!$       x=0.0d0
!!$       y=1.0d0
!!$       z=0.0d0
!!$       call symop(k,x,y,z,x1,y1,z1,oname,nz,lattic) 
!!$       rot(1,2)=x1
!!$       rot(2,2)=y1
!!$       rot(3,2)=z1
!!$       
!!$       x=0.0d0
!!$       y=0.0d0
!!$       z=1.0d0
!!$       call symop(k,x,y,z,x1,y1,z1,oname,nz,lattic) 
!!$       rot(1,3)=x1
!!$       rot(2,3)=y1
!!$       rot(3,3)=z1
!!$       
!!$       if (oname.ne.'1') then
!!$          if (dbr.lt.0.00001) then 
!!$             is=is+1
!!$             symo(1:3,1:3,is)=transpose(rot(1:3,1:3))
!!$             tsymo(1:3,is)=(/0.0d0,0.0d0,0.0d0/)
!!$
!!$!write(0,*) is
!!$!write(0,*) rot(1:3,1:3)
!!$          endif
!!$       endif
!!$       
!!$    enddo
    
    iat0=0
    nt=0
    do i=1,nat
       if (i.gt.1) iat0=iat0+mult(i-1)
       do j1=1,mult(i)
          iat1=iat0+j1
          do j2=1,mult(i)
             iat2=iat0+j2
             do i1=1,is
                do i2=1,is          
                   xyz1(1:3)=pos(1:3,iat1)
                   xyz2(1:3)=pos(1:3,iat2)
                   sym1(1:3,1:3)=symo(1:3,1:3,i1)
                   sym1=transpose(sym1)
                   sym2(1:3,1:3)=symo(1:3,1:3,i2)
                   sym2=transpose(sym2)
                   t=matmul(sym1,xyz1)-matmul(sym2,xyz2)
                   t=t+10.0d0
                   t(1:3)=mod(t(1:3),1.0d0)
                   do k=1,3
                      if (abs(t(k)-1.0d0).lt.small) t(k)=0.0d0
                   enddo
                   if (dot_product(t,t).lt.small) goto 1
                   do k=1,nct
                      tt(1:3)=t(1:3)-ct(1:3,k)
                      if (dot_product(tt,tt).lt.small) goto 1
                   enddo
                   do k=1,nt
                      tt(1:3)=t(1:3)-ts(1:3,k)
                      if (dot_product(tt,tt).lt.small) goto 1
                   enddo
! 2004.02.04 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   do k=1,nt
                      do kk=1,nct
                         tt(1:3)=t(1:3)+ct(1:3,kk)-ts(1:3,k)
                         tt=tt+10
                         tt(1:3)=mod(tt(1:3),1.0d0)
                         if (dot_product(tt,tt).lt.small) goto 1
                      enddo
                   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   nt=nt+1
                   if (nt.gt.mnsym) then
                      write(6,*) 'symgen: nt.gt.mnsym'
                      write(6,'(a,2i8,3f10.5)') 'tran',nt,mnsym,t
                      stop
                   endif
                   ts(1:3,nt)=t(1:3)
1                  continue
                enddo
             enddo
          enddo
       enddo
    enddo

    mnsym=(nsymp)*(nt+1)
    symos=symo
    tsymos=tsymo 
    deallocate(symo,tsymo)       
    allocate (symo(3,3,mnsym),tsymo(3,mnsym))

    symo(1:3,1:3,1:nsymp)=symos(1:3,1:3,1:nsymp)
    tsymo(1:3,1:nsymp)=tsymos(1:3,1:nsymp)

    ns=is    
    do i=1,ns
       do j=1,nt
          is=is+1
          if (is.gt.mnsym) then
             write(0,*) 'symgen: is.gt.mnsym'
             stop
          endif 
          symo(1:3,1:3,is)=symo(1:3,1:3,i)
          tsymo(1:3,is)=ts(1:3,j)
       enddo
    enddo

10  continue

    deallocate (iz,tau,inum,tinv,ct)
    allocate(iz(3,3,is+1),tau(3,is+1),inum(is+1),tinv(is+1))

    iz(1:3,1:3,1:is)=symo(1:3,1:3,1:is)
    tau(1:3,1:is)=tsymo(1:3,1:is)
    do i=1,is
       inum(i)=i
    end do
    nsym=is
    iord=is
    
 end subroutine symgen
