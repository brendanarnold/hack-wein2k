!BOP
! !ROUTINE: CalculateMatrices
! !INTERFACE:
      SUBROUTINE CalculateMatrices
! !USES:
      use input
	  use rotation_matrices
	  use struct
	  use program_control,only : verbosity
! !DESCRIPTION:
!     build the matrices GeneralM(i) for each equivalent atom i:
!     GeneralM(i) * Q1 = Q2
!     Q1 in the basis of the laboratory, and Q2 in the basis of the atom i.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
	  implicit none
      
      DOUBLE PRECISION K1(3), K2(3), K3(3)
      DOUBLE PRECISION Kout1(3), Kout2(3), Kout3(3)
      DOUBLE PRECISION Krot1(3), Krot2(3), Krot3(3)
      DOUBLE PRECISION IntM(3, 3), MM(3, 3), Det
      INTEGER IEqAt, I, J
      
	  write(6,'(/,a)') 'Output from subroutine CalculateMatrices: '
!     K are in the x, y, z direction, in the labor basis:
      K1(1) = dble(1)
      K1(2) = dble(0)
      K1(3) = dble(0)
      K2(1) = dble(0)
      K2(2) = dble(1)
      K2(3) = dble(0)
      K3(1) = dble(0)
      K3(2) = dble(0)
      K3(3) = dble(1)
      
!     calcule K in the carthesian basis of the crystal:
      CALL RotateCrist (K1, -AGamma, -ABeta, -AAlpha)
      CALL RotateCrist (K2, -AGamma, -ABeta, -AAlpha)
      CALL RotateCrist (K3, -AGamma, -ABeta, -AAlpha)

      if(verbosity.ge.1) then
        WRITE(6,'(a)') 'from the laboratory basis to the crystal basis :' 
        DO I=1, 3 
           WRITE(6,'(3(f8.3,x))') K1(I), K2(I), K3(I)
        ENDDO
	  endif
      
            
!     we put K in the reciprocal lattice basis -> Kout:
      CALL CalculateInterMatrix (ALAT(1), ALAT(2), ALAT(3),  &
           ALPHA(1), ALPHA(2), ALPHA(3), &
           IntM(1,1))

      if(verbosity.ge.1) then
        WRITE(6,'(a)') 'from the crystal basis to the reciprocal basis:' 
        DO I=1, 3 
           WRITE(6,'(3(f8.3,x))') (IntM(I, J), J=1, 3)
        ENDDO
	  endif
      
      CALL ProductMatVect (IntM, K1, Kout1)
      CALL ProductMatVect (IntM, K2, Kout2)
      CALL ProductMatVect (IntM, K3, Kout3)
      
      DO IEqAt = NEqAtDeb, NEqAt                 ! used to be : 1,NEqAt
!     for each atom, we rotate Kout in the basis of the atom:

         CALL ROTATE(Kout1, ROTIJ(1, 1, IEqAt), Krot1)
         CALL ROTATE(Kout2, ROTIJ(1, 1, IEqAt), Krot2)
         CALL ROTATE(Kout3, ROTIJ(1, 1, IEqAt), Krot3)

!     we put krot in a carthesian basis:
         CALL UseBr1(krot1, K1)
         CALL UseBr1(krot2, K2)
         CALL UseBr1(krot3, K3)
         
!     we rotate it (rotloc for all equivalents atoms:
         CALL ROTATE(K1, ROTLOC(1, 1, NATOM), Krot1)
         CALL ROTATE(K2, ROTLOC(1, 1, NATOM), Krot2)
         CALL ROTATE(K3, ROTLOC(1, 1, NATOM), Krot3)
         
         
!     we write the general transformation matrix GeneralM,
!     GeneralM_i * V_labor = V_atom_i
         WRITE(6,'(a,i4,a)') 'General Matrix for atom', IEqAt, ':'
         DO I=1, 3 
            GeneralM(I, 1, IEqAt) = Krot1(I)
            GeneralM(I, 2, IEqAt) = Krot2(I)
            GeneralM(I, 3, IEqAt) = Krot3(I)
            WRITE(6,'(3(f8.3,x))') (GeneralM(I, J, IEqAt), J=1, 3)
         ENDDO
         CALL Determinant(GeneralM(1, 1, IEqAt), Det)
		 if(verbosity.ge.1) WRITE(6,'(a,f12.4)') 'Determinant :', Det
      ENDDO
      
      RETURN
      END


