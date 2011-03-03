      SUBROUTINE Broad(prev,jemax,ene,XInput, XOutput)
!     a broadening with the Gauss function (S = Sigma = PreV): 
!                                  2
!                                 x
!                             -  ----
!                                   2    
!                  1             2 S   
!     Fct(x) = ----------   e
!                   +----
!              S  \/ 2 Pi
!
!             \____  ____/
!                  \/                                 2
!                = Factor1/dx   and Factor2 = 2 (S/dx)
!     Precision = PreV / dx. (dx = Delta Energy between 2 points)
!     Broad. Spectrum(x) = 
!     Sum (x'=-3*Precision -> 3*Precision) of Fct(x')*Spectrum(x+x')*dx,
!     where x' = DeltaJ in the procedure.
!
! -> WARNING: at the place where the spectrum is not calculated, we take 0
!
!      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'param.inc'
!      INCLUDE 'datael.inc'

      real*4 ene(jemax),XInput(jemax), XOutput(jemax)

!     Internal variables:
      real*8 Factor1, Factor2
      INTEGER Precision, DeltaJ
      pi=acos(-1.d0)
      MinJEnergy=1
      jemin=1
!-------reset !
	DO J=MinJEnergy, JEMAX
 	  XOutput(J)=0.d0
	ENDDo

!**      WRITE (*,*) 'MinJEnergy =', MinJEnergy
      
!      IF (PreV.GT.0.d0) THEN
         Precision = NINT(ABS(PreV/(ENE(JEMIN+1)-ENE(JEMIN))))
      IF (Precision.GT.0) THEN
         Factor1 = 1.d0 / DBLE(Precision) / SQRT(2.d0*PI)
         Factor2 = DBLE (2.d0 * Precision * Precision)
              write(6,*) 'smearing', Precision,Factor1 , Factor2
         DO J=MinJEnergy, JEMAX
            DO DeltaJ= MAX(-3*Precision, JEMIN-J),  &
                 MIN(3*Precision, JEMAX-J)
               XOutput(J) = XOutput(J)+ &
                    XInput(J+DeltaJ)* &
                    DEXP(-DBLE(DeltaJ*DeltaJ)/Factor2)
            ENDDO
            XOutput(J) = XOutput(J)*Factor1
         ENDDO
      else
	DO J=MinJEnergy, JEMAX
 	  XOutput(J)=XInput(j)
	ENDDo
      ENDIF

	RETURN 

      END
