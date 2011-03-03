      subroutine findequivs(pred)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax),pall(3,natmax*48)
      parameter (tol=0.25)
!     Find operations that regenerate reduced set sites
!
!     Expand the set
      call expandset(pred,pall)
!
      DO JATOM=1,NAT
!	 First equivalent atom
	 IF(JATOM .gt. 1)then
		IND=IND+MULT(JATOM-1)
	 else
		IND=1
	 endif
!	 Last
	 IUP=IND+MULT(JATOM)-1
!
         X00=pred(1,JATOM)
         Y00=pred(2,JATOM)
         Z00=pred(3,JATOM)
         X0=0
         Y0=0
         Z0=0
         NTHIS=0
!        Loop over equivalents
 	  DO M2=IND,IUP
            DO 10 I=1,IORD
            DO 10 I1=-3,3
            DO 10 I2=-3,3
            DO 10 I3=-3,3
!       Generate equivalent sites
            X=TAU(1,I)+I1
            Y=TAU(2,I)+I2
            Z=TAU(3,I)+I3
            DO J=1,3
                X=X+dble(IZ(J,1,I))*pall(J,M2)
                Y=Y+dble(IZ(J,2,I))*pall(J,M2)
                Z=Z+dble(IZ(J,3,I))*pall(J,M2)
            ENDDO
!       How far is it, in a.u.
            X1=ABS(X-X00)*alat(1)
            Y1=ABS(Y-Y00)*alat(2)
            Z1=ABS(Z-Z00)*alat(3)

            IF((X1.LT.tol).AND.(Y1.LT.tol).AND.(Z1.LT.tol)) THEN
!               It is close, so add to list
                NTHIS=NTHIS+1
                indeq(1,JATOM,NTHIS)=I
                indeq(2,JATOM,NTHIS)=M2
                indeq(3,JATOM,NTHIS)=I1
                indeq(4,JATOM,NTHIS)=I2
                indeq(5,JATOM,NTHIS)=I3
                nineq(JATOM)=NTHIS
            ENDIF

10        CONTINUE
          ENDDO
            write(6,1)jatom,nthis,iord
1           format('Atom ',i2,' has ',i3,' equivalents ', i2,' operations')
      ENDDO
      return
      end
      subroutine genequivs(pred,pose)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax),pall(3,natmax*48),pose(3,natmax*58)
      parameter (tol=0.25)
!     Generate symmetry equivalents of reduced position
!
!     Expand the set
      call expandset(pred,pall)
!
      II=0
      DO JATOM=1,NAT
         DO N=1,nineq(JATOM)
          II=II+1
          I=indeq(1,JATOM,N)
          M2=indeq(2,JATOM,N)
!       Generate equivalent sites
            X=TAU(1,I)+indeq(3,JATOM,N)
            Y=TAU(2,I)+indeq(4,JATOM,N)
            Z=TAU(3,I)+indeq(5,JATOM,N)
            DO J=1,3
                X=X+dble(IZ(J,1,I))*pall(J,M2)
                Y=Y+dble(IZ(J,2,I))*pall(J,M2)
                Z=Z+dble(IZ(J,3,I))*pall(J,M2)
            ENDDO
            POSE(1,II)=X
            POSE(2,II)=Y
            POSE(3,II)=Z
	ENDDO
      ENDDO
      return
      end
      subroutine daverage(pred)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax),pall(3,natmax*48)
      real*8 X,Y,Z,X0,Y0,Z0,t
      parameter (tol=0.25)
!     Generate symmetry equivalents of reduced position
!
!     Expand the set
      call expandset(pred,pall)
!
      DO JATOM=1,NAT
         X0=0
         Y0=0
         Z0=0
         DO N=1,nineq(JATOM)
          I=indeq(1,JATOM,N)
          M2=indeq(2,JATOM,N)
!       Generate equivalent sites
            X=TAU(1,I)+indeq(3,JATOM,N)
            Y=TAU(2,I)+indeq(4,JATOM,N)
            Z=TAU(3,I)+indeq(5,JATOM,N)
            DO J=1,3
                X=X+dble(IZ(J,1,I))*pall(J,M2)
                Y=Y+dble(IZ(J,2,I))*pall(J,M2)
                Z=Z+dble(IZ(J,3,I))*pall(J,M2)
            ENDDO
            X0=X0+X
            Y0=Y0+Y
            Z0=Z0+Z
	ENDDO
	t=1.D0/dble(nineq(JATOM))
	pred(1,JATOM)=X0*t
	pred(2,JATOM)=Y0*t
	pred(3,JATOM)=Z0*t
      ENDDO
      return
      end
