!*********************************************************
!********* INPUT: Cubic lattice parameter a                                   
!*********        Strain factor eps (ct/at=(1+eps))
!*********       
!********* OUTPUT: lattice parameters for the volume conservative
!*********         tetragonal strained lattice at,ct
!*********
!*********
!*********

      SUBROUTINE c2te(a,eps,at,ct)
      REAL*8 a,eps,at,ct,u,v

      u=(1.d0+eps)**(-1.d0/3.d0)
      v=(1.d0+eps)**(2.d0/3.d0)

      at=a*u
      ct=a*v



      RETURN
      END
