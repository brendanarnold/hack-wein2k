!.....calculates the inverse AINV of a matrix A...........................
!
      SUBROUTINE INVmat(A,AINV)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension A(3,3),AINV(3,3)
!
!     write(*,*) 'in invamt: a'
!     write(6,5060) ((A(I1,I2),I2=1,3),I1=1,3)

        det= a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
            +a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3) &
            -a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
!
      if (abs(det).lt.1.d-5) then
        write(6,*) 'determinant = 0'
        WRITE(6,5060) ((a(L,M), m=1,3),l=1,3)
        return            
      end if
!
      AINV(1,1)=(   A(2,2) * A(3,3) - A(2,3) * A(3,2) )/det
      AINV(2,1)=( - A(2,1) * A(3,3) + A(2,3) * A(3,1) )/det
      AINV(3,1)=(   A(2,1) * A(3,2) - A(2,2) * A(3,1) )/det
      AINV(1,2)=( - A(1,2) * A(3,3) + A(1,3) * A(3,2) )/det
      AINV(2,2)=(   A(1,1) * A(3,3) - A(1,3) * A(3,1) )/det
      AINV(3,2)=( - A(1,1) * A(3,2) + A(1,2) * A(3,1) )/det
      AINV(1,3)=(   A(1,2) * A(2,3) - A(1,3) * A(2,2) )/det
      AINV(2,3)=( - A(1,1) * A(2,3) + A(1,3) * A(2,1) )/det
      AINV(3,3)=(   A(1,1) * A(2,2) - A(1,2) * A(2,1) )/det

      WRITE(6,222) (a(1,M), m=1,3),(ainv(1,M), m=1,3)
      WRITE(6,222) (a(2,M), m=1,3),(ainv(2,M), m=1,3)
      WRITE(6,222) (a(3,M), m=1,3),(ainv(3,M), m=1,3)

 222  format(5x,3f10.7,4x,3f10.7)
 5060 FORMAT(20X,3F10.7)
      RETURN
      END
