      
      
!********************************************************************
      SUBROUTINE MultScalMatrix(S, M)
!     multiply the matrix M with S -> M
      
      DOUBLE PRECISION M(3,3), S
      INTEGER I, J
      
      DO I=1, 3
         DO J=1, 3
            M(I, J) = M(I, J) * S
         ENDDO
      ENDDO
      
      RETURN
      END
      

