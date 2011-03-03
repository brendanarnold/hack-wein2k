      SUBROUTINE ANGLE(XMS,THETA,PHI)
      USE param
      USE struct
      USE case
      IMPLICIT REAL*8 (A-H,O-Z)
!*******************************************************************
      LOGICAL          ORTHO
!
      COMMON/ORTH/     ORTHO
      DIMENSION XMS(3)

	PI=ACOS(-1.D0)
!---------------------------------------------------------------------
	write(6,*)'ortho=',ortho
	IF (ORTHO.OR.(LATTIC(1:1).EQ.'R')) THEN
        XA=AA*XMS(1)
	XB=BB*XMS(2)
	XC=CC*XMS(3)

	ELSE
      IF(LATTIC(1:1).EQ.'H') THEN
	XA=XMS(1)*AA*SQRT(3.D0)/2
	XB=AA*(XMS(2)-XMS(1)/2)
	XC=CC*XMS(3)

      ELSE

	IF (ABS(ALPHA(3)-PI/2).GT.1D-4) THEN
        XA=XMS(1)*AA*SIN(ALPHA(3))
        XB=XMS(1)*AA*COS(ALPHA(3))+BB*XMS(2)
        XC=CC*XMS(3)

        ELSE IF (ABS(ALPHA(2)-PI/2).GT.1D-4) THEN
	XA=XMS(1)*AA*SIN(ALPHA(2))
	XB=XMS(2)*BB
	XC=XMS(1)*AA*COS(ALPHA(2))+CC*XMS(3)

	ELSE
!	WRITE(6,*)'EXCHANGE THE LATTICE VECTORS, ALPHA(1) MUST BE PI/2'
!        END IF

	TT=ABS((ALPHA(1)-PI/2)*(ALPHA(2)-PI/2)*(ALPHA(3)-PI/2))
!	IF (TT.GT.1D-3) STOP 'TRICLINIC NOT IMPLEMENTED'
        write(6,*) ' Triclinic implemented, but never tested'
        cosg1=(cos(ALPHA(3))-cos(alpha(1))*cos(ALPHA(2)))/sin(alpha(1))/sin(alpha(2))
        gamma0=acos(cosg1)
!      from lapw5
!      BR2(1,1)=A(1)*1.0d0*sin(gamma0)*sin(beta)
!      BR2(1,2)=A(1)*1.0d0*cos(gamma0)*sin(beta)
!      BR2(1,3)=A(1)*1.0d0*cos(beta)    
!      BR2(2,1)=0.0d0              
!      BR2(2,2)=A(2)*1.0d0*sin(alpha)
!      BR2(2,3)=A(2)*1.0d0*cos(alpha)
!      BR2(3,1)=0.0d0              
!      BR2(3,2)=0.0d0              
!      BR2(3,3)=A(3)*1.0d0              
        XA=XMS(1)*sin(gamma0)*sin(alpha(2))*AA
        XB=XMS(1)*cos(gamma0)*sin(alpha(2))*AA + XMS(2)*sin(alpha(1))*BB
        XC=XMS(1)*cos(alpha(2))*AA + XMS(2)*cos(alpha(1))*BB + XMS(3)*CC
        END IF

      END IF
 	ENDIF
     
        XX=SQRT(XA**2+XB**2+XC**2)
        THETA=ACOS(XC/XX)
	IF (ABS(THETA).LT.1D-5) THEN 
	PHI=0.0
	ELSE
	XX=SQRT(XA**2+XB**2)
	PHI=ACOS(XA/XX)
	IF (ABS(XB).GT.1D-5) PHI=PHI*XB/ABS(XB)
	END IF

 	WRITE(6,*)'LATTICE=',LATTIC

!...... test of collinearity of the z-axis with spin quantization axis
        XX=SQRT(XA**2+XB**2+XC**2)
        do i=1,natom
	ia=iatom(i)
	sum=(rotloc(3,1,ia)*xa+rotloc(3,2,ia)*xb+rotloc(3,3,ia)*xc)/xx
	if (abs(sum-1.).gt.1d-5) then
	WRITE(6,*)'Z-AXIS NON-COLLINEAR WITH SPIN QUANTIZATION AXIS FOR ATOM',IA  
!       WRITE(6,*)'SPIN AXIS IS:',XA/XX,XB/XX,XC/XX
	WRITE(6,*)'Z_loc.S_z=',sum
!	WRITE(6,*) 'Z-AXIS NON-COLLINEAR WITH SPIN QUANTIZATION AXIS'
!	WRITE(6,*)'WARNING - SYMMETRIZATION MUST BE USED! CHECK IT!'
	end if
	end do
	END 
	



