!BOP
! !ROUTINE: CalculateInterMatrix
! !INTERFACE:
      SUBROUTINE CalculateInterMatrix (AA, BB, CC, Alpha, Beta, Gamma,  IntM)
! !USES:
	  use program_control,only : verbosity
! !INPUT/OUTPUT PARAMETERS:
!     aa,bb,cc         :  the lattice parameters of the crystal in a.u.
!     alpha,beta,gamma :  the lattice angles of the crystal in rad
!     IntM             :  transformation matrix
! !DESCRIPTION:
!     we have Q in the carthesian labor-basis. 
!     - Euler angles -> in the carthesian basis of the crystal,
!     - IntM -> in the basis of the reciprocal lattice
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  IN/OUTPUT	        
      real*8,intent(in)  :: AA, BB, CC, Alpha, Beta, Gamma
      real*8,intent(out) :: IntM(3, 3)
!  LOCALS      
!     Al, Bl, Cl are the lattice basis vectors 
      DOUBLE PRECISION Al(3), Bl(3), Cl(3)
!     Ar, Br, Cr are the reciprocal lattice vectors
      DOUBLE PRECISION Ar(3), Br(3), Cr(3)
!     a matrix and a vector for the calculs:
      DOUBLE PRECISION M(3, 3), VInt(3)
!     Sinus and cosinus:
      DOUBLE PRECISION CA, CB, CG, SA, SB, SG
      DOUBLE PRECISION Volume
      INTEGER I
      
      CA = DCOS (Alpha)
      SA = DSIN (Alpha)
      CB = DCOS (Beta)
      SB = DSIN (Beta)
      CG = DCOS (Gamma)
      SG = DSIN (Gamma)

!     we define the 3 lattice vectors in the crystal basis:
      Al(1) = AA * DSQRT(1-CB**2-((CG-CA*CB)/SA)**2)
      Al(2) = AA * (CG-CA*CB)/SA
      Al(3) = AA * CB
      
      Bl(1) = dble(0)
      Bl(2) = BB * SA
      Bl(3) = BB * CA
      
      Cl(1) = dble(0)
      Cl(2) = dble(0)
      Cl(3) = CC

      if(verbosity.ge.1) then
         write(6,'(/,a)') 'Output of subroutine CalculInterMatrix:'
         write(6,'(a)') 'lattice vectors:'
        !WRITE(6,'(3(f8.3,x))') 'Al = ', (Al(I), I=1, 3)
        !WRITE(6,'(3(f8.3,x))') 'Bl = ', (Bl(I), I=1, 3)
        !WRITE(6,'(3(f8.3,x))') 'Cl = ', (Cl(I), I=1, 3)
        ! jl ifort doesn't like this format:-/
         WRITE(6,*) 'Al = ', (Al(I), I=1, 3)
         WRITE(6,*) 'Bl = ', (Bl(I), I=1, 3)
         WRITE(6,*) 'Cl = ', (Cl(I), I=1, 3)
      endif
      
!     we calculate the reciprocal lattice vectors (in the crystal basis):
      CALL ProdVectVect (Bl, Cl, VInt)
      Volume = Vint(1) * Al(1) + Vint(2) * Al(2) + Vint(3) * Al(3)
      CALL ReciprocLatVect (Bl, Cl, Volume, Ar)
      CALL ReciprocLatVect (Cl, Al, Volume, Br)
      CALL ReciprocLatVect (Al, Bl, Volume, Cr)

      if(verbosity.ge.1) then
        write(6,'(a)') 'reciprocal lattice vectors:'
        !WRITE(6,'(3(f8.3,x))') 'Ar = ', (Ar(I), I=1, 3)
        !WRITE(6,'(3(f8.3,x))') 'Br = ', (Br(I), I=1, 3)
        !WRITE(6,'(3(f8.3,x))') 'Cr = ', (Cr(I), I=1, 3)
        !jl same as above :-//
        WRITE(6,*) 'Ar = ', (Ar(I), I=1, 3)
        WRITE(6,*) 'Br = ', (Br(I), I=1, 3)
        WRITE(6,*) 'Cr = ', (Cr(I), I=1, 3)
      endif
            
!     define the matrice M: M*V = V'; V in the reciprocal lattice basis,
!     and V' in the carthesian basis of the crystal.
      DO I=1, 3
         M(I,1) = Ar(I) !Ar
         M(I,2) = Br(I) !Br
         M(I,3) = Cr(I) !Cr
      ENDDO
      
!     calculate the inverse of the matrix M:
      CALL Inverse(M, IntM)

!     for this matrix IntM, IntM*V' = V

      RETURN
      END
      
