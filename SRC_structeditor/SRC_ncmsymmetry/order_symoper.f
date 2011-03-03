subroutine order_symoper

 use struct
 use rotations, only: sspiral

 integer one(3,3),i,sym(3,3),j
 logical matrixeq,rvecteq
 real*8   notr(3),tr(3)

 one(1,1:3)=(/1,0,0/)
 one(2,1:3)=(/0,1,0/)
 one(3,1:3)=(/0,0,1/)
 notr(1:3)=(/0.0d0,0.0d0,0.0d0/)

 j=iord
 do i=iord,1,-1
    sym(1:3,1:3)=iz(1:3,1:3,i)
    tr=tau(1:3,i)
    if (.not.(matrixeq(one,sym).and.rvecteq(tr,notr))) then
       iz(1:3,1:3,j)=iz(1:3,1:3,i)
       tau(1:3,j)=tau(1:3,i)
       j=j-1
    endif
 enddo
 if (j.ne.1) then
    write(6,*) 'cos nie tak w order_symoper',j 
    stop 'cos nie tak w order_symoper'
 endif
 iz(1:3,1:3,1)=one(1:3,1:3)
 tau(1:3,1)=0.0d0
 
 if (sspiral) then
    nsym=1
    iord=1
 endif

end subroutine order_symoper

logical function matrixeq(sym1,sym2)

  integer sym1(3,3),sym2(3,3)
  integer i,j
  
  matrixeq=.true.
  do i=1,3
     do j=1,3
        if (sym1(i,j).ne.sym2(i,j)) then
           matrixeq=.false.
           return
        endif
     enddo
  enddo

end function matrixeq

  logical function rvecteq(v1,v2)

    real*8  v1(3),v2(3)
    real*8    small
    logical ok
    integer i

    small=1.0d-8

    ok=.true.
    do i=1,3
       if(abs(v1(i)-v2(i)).gt.small) ok=.false.
    enddo
    rvecteq=ok

  end function rvecteq
