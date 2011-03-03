!.........DEFINE INVERS SYMMETRY OPERATION............................
      SUBROUTINE INVERSSYMDEF(A,AINV)
      INTEGER A(3,3),AINV(3,3),DET
        det= a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
            +a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3) &
            -a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
      if (det.eq.0) then
                    write(6,*) 'NEGATIV DET'
        WRITE(6,222) ((a(L,M),ainv(L,M), m=1,3),l=1,3)
 222  FORMAT(3(5X,3I2,4X,3I2,/))
                    end if
      AINV(1,1)=(   A(2,2) * A(3,3) - A(2,3) * A(3,2) )/det
      AINV(2,1)=( - A(2,1) * A(3,3) + A(2,3) * A(3,1) )/det
      AINV(3,1)=(   A(2,1) * A(3,2) - A(2,2) * A(3,1) )/det
      AINV(1,2)=( - A(1,2) * A(3,3) + A(1,3) * A(3,2) )/det
      AINV(2,2)=(   A(1,1) * A(3,3) - A(1,3) * A(3,1) )/det
      AINV(3,2)=( - A(1,1) * A(3,2) + A(1,2) * A(3,1) )/det
      AINV(1,3)=(   A(1,2) * A(2,3) - A(1,3) * A(2,2) )/det
      AINV(2,3)=( - A(1,1) * A(2,3) + A(1,3) * A(2,1) )/det
      AINV(3,3)=(   A(1,1) * A(2,2) - A(1,2) * A(2,1) )/det
      RETURN
      END
