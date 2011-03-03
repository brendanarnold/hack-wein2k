   subroutine inversa(a,ainv)
     implicit real*8 (a-h,o-z)
     dimension a(3,3),ainv(3,3)
     det= a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
          +a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3) &
          -a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
     ainv(1,1) =(   a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
     ainv(2,1) =( - a(2,1) * a(3,3) + a(2,3) * a(3,1) ) / det
     ainv(3,1) =(   a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
     ainv(1,2) =( - a(1,2) * a(3,3) + a(1,3) * a(3,2) ) / det
     ainv(2,2) =(   a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
     ainv(3,2) =( - a(1,1) * a(3,2) + a(1,2) * a(3,1) ) / det
     ainv(1,3) =(   a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det
     ainv(2,3) =( - a(1,1) * a(2,3) + a(1,3) * a(2,1) ) / det
     ainv(3,3) =(   a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det
     return
   end subroutine inversa

