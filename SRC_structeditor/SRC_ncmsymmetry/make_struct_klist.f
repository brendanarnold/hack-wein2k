 subroutine make_struct_klist

    use struct

    integer i

    do i=1,nsym
       if (tinv(i).eq.1) iz(1:3,1:3,i)=iz(1:3,1:3,i)*(-1) 
    enddo

 end subroutine make_struct_klist
