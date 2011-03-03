!BOP
! !ROUTINE: Rotation
! !INTERFACE:
    subroutine Rotation (V, Alpha, Beta, Gamma, Vout)
! !INPUT/OUTPUT PARAMETERS:
!   v    :   vector before rotation
!   vout :   vector after rotation
!   alpha,beta,gamma : Euler angles defining the rotation
! !DESCRIPTION:
!   Rotation of a vector, where the rotation is defined by Euler angles.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!   IN/OUTPUT
      real*8,intent(in)  :: Alpha, Beta, Gamma,V(3)
      real*8,intent(out) :: Vout(3)
!   LOCALS  : rotation matrix
      real*8 A(3,3), B(3,3), C(3,3), R(3,3), Rint(3,3)
      
!     reset matrices to A(i,j)=delta(i,j)
      call identity (A)
      call identity (B)
      call identity (C)
!     A : rotation matrice around z: (x,y,z)->(x',y',z'=z)
      A(1,1) = dcos(Alpha)
      A(1,2) = -dsin(Alpha)
      A(2,2) = dcos(Alpha)
      A(2,1) = dsin(Alpha)
!     B : rotation matrice around y': (x',y',z')->(x",y"=y',z")
      B(1,1) = dcos(Beta)
      B(1,3) = dsin(Beta)
      B(3,3) = dcos(Beta)
      B(3,1) = -dsin(Beta)
!     C: rotation matrice around z": (x",y",z")->(x"',y"',z"'=z") 
      C(1,1) = dcos(Gamma)
      C(1,2) = -dsin(Gamma)
      C(2,2) = dcos(Gamma)
      C(2,1) = dsin(Gamma)

!     We set in R the product A * B * C . Vout = R * V
      call ProductMatMat(A,B,Rint)
      call ProductMatMat(Rint,C,R)
      call ProductMatVect(R,V,Vout)

      return
      end

