 subroutine determinant(sym,det)
  
  real*8   sym(3,3),det

  det=sym(1,1)*(sym(2,2)*sym(3,3)-sym(3,2)*sym(2,3))-&
      sym(1,2)*(sym(2,1)*sym(3,3)-sym(3,1)*sym(2,3))+&
      sym(1,3)*(sym(2,1)*sym(3,2)-sym(2,2)*sym(3,1))

 end subroutine determinant

 subroutine determinant_c(sym,det)
  
  complex*16   sym(3,3),det

  det=sym(1,1)*(sym(2,2)*sym(3,3)-sym(3,2)*sym(2,3))-&
      sym(1,2)*(sym(2,1)*sym(3,3)-sym(3,1)*sym(2,3))+&
      sym(1,3)*(sym(2,1)*sym(3,2)-sym(2,2)*sym(3,1))

end subroutine determinant_c



 subroutine vecprod(a,b,c)
   implicit real*8 (a-h)
   dimension a(3),b(3),c(3)
   
   c(1)=a(2)*b(3)-a(3)*b(2)
   c(2)=a(3)*b(1)-a(1)*b(3)
   c(3)=a(1)*b(2)-a(2)*b(1)

 end subroutine vecprod


 
 real*8 function dot(a,b)
   implicit real*8 (a-h)
   dimension a(3),b(3)
   dot=0
   do i=1,3
      dot=dot+a(i)*b(i)
   end do
 end function dot
 
