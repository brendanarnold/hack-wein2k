      SUBROUTINE ROTATE(VECTOR,ROTMAT,ROTVEC)
!
!        Arguments
!
      DOUBLE PRECISION   ROTVEC(3), VECTOR(3)
      DOUBLE PRECISION   ROTMAT(3,3)
!
!     ..................................................................
!
!        'ROTATE' performs a rotation of the vector from the general
!        cartesian coordination system into the local one of the
!        JATOM-th sphere.
!
!     ..................................................................
!
!        Local Scalars
!
      INTEGER            J, JCOORD
      DOUBLE PRECISION   DOTPRO
!
      DO 20 JCOORD = 1, 3
         DOTPRO = 0.0
         DO 10 J = 1, 3
            DOTPRO = DOTPRO + VECTOR(J)*ROTMAT(JCOORD,J)
   10    CONTINUE
         ROTVEC(JCOORD) = DOTPRO
   20 CONTINUE
!
      RETURN
!
!        End of 'ROTATE'
!
      END
