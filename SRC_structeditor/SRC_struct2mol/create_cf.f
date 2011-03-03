subroutine create_cf

  use struct
  
  integer i
  real*8 rnic,rcos

  rnic=0.0d0
  rcos=1.0d0

  do i=1,nat
     write(5,*) 
     write(5,'(3f10.5)') rnic,rnic,rcos
  enddo
  rewind(5)

end subroutine create_cf

