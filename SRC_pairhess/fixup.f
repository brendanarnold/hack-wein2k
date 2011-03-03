      subroutine fixup(pred)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax),pold(3,natmax)
      dimension pnew(3,48),tmp(3)
      parameter (tol=0.1D0, tol2=-1D-4, tol3=1.D0+tol2)
      save debug,icnt
      logical debug
      data debug/.true./,icnt/0/
! Patch up the positions to remove symmetry issues
! Note: two cycles are probably not needed, but why not?
      if(debug)write(79,*)'Icnt ',icnt
      if(debug)pold(1:3,1:nat)=pred(1:3,1:nat)
!      DO ITEST=1,2
      DO JATOM=1,nat
         NTHIS=0
         X00=pred(1,JATOM)
         Y00=pred(2,JATOM)
         Z00=pred(3,JATOM)
         X0=0
         Y0=0
         Z0=0
            DO I=1,IORD
!       Generate equivalent sites
            X=TAU(1,I)+1.D0
            Y=TAU(2,I)+1.D0
            Z=TAU(3,I)+1.D0
            DO J=1,3
                X=X+dble(IZ(J,1,I))*pred(J,JATOM)
                Y=Y+dble(IZ(J,2,I))*pred(J,JATOM)
                Z=Z+dble(IZ(J,3,I))*pred(J,JATOM)
            ENDDO
            pnew(1,i)=MOD(X,1.0D0)
            pnew(2,i)=MOD(Y,1.0D0)
            pnew(3,i)=MOD(Z,1.0D0)
            enddo

           DO L=1,IORD
            DO 10 I1=-1,1
            TMP(1)=pnew(1,L)+I1
            DO 10 I2=-1,1
            TMP(2)=pnew(2,L)+I2
            DO 10 I3=-1,1
            TMP(3)=pnew(3,L)+I3
            DO I=1,IORD
!       Generate equivalent sites
            X=TAU(1,I)
            Y=TAU(2,I)
            Z=TAU(3,I)
            DO J=1,3
                X=X+dble(IZ(J,1,I))*TMP(J)
                Y=Y+dble(IZ(J,2,I))*TMP(J)
                Z=Z+dble(IZ(J,3,I))*TMP(J)
            ENDDO
            X=MOD(X,1.0D0)
            Y=MOD(Y,1.0D0)
            Z=MOD(Z,1.0D0)
!       How far is it, in a.u.
            X1=ABS(X-X00)*alat(1)
            Y1=ABS(Y-Y00)*alat(2)
            Z1=ABS(Z-Z00)*alat(3)
            IF((X1.LT.TOL).AND.(Y1.LT.TOL).AND.(Z1.LT.TOL)) THEN
!               It is close, so add to average
                NTHIS=NTHIS+1
                X0=X0+X
                Y0=Y0+Y
                Z0=Z0+Z
            ENDIF
            ENDDO
10          CONTINUE
           ENDDO
       if(nthis.gt.0)then
!       Final averaged value
        TT=1.D0/dble(nthis)
        pred(1,JATOM)=X0*TT
        pred(2,JATOM)=Y0*TT
        pred(3,JATOM)=Z0*TT
        if(pred(1,JATOM) .lt. TOL2)pred(1,jatom)=pred(1,jatom)+1.D0
        if(pred(2,JATOM) .lt. TOL2)pred(2,jatom)=pred(2,jatom)+1.D0
        if(pred(3,JATOM) .lt. TOL2)pred(3,jatom)=pred(3,jatom)+1.D0
        if(pred(1,JATOM) .gt. TOL3)pred(1,jatom)=pred(1,jatom)-1.D0
        if(pred(2,JATOM) .gt. TOL3)pred(2,jatom)=pred(2,jatom)-1.D0
        if(pred(3,JATOM) .gt. TOL3)pred(3,jatom)=pred(3,jatom)-1.D0
        if(debug)write(79,79)jatom,nthis,(pred(k,jatom),pold(k,jatom),k=1,3)
79      format(2i4,6F12.8)
       endif
      ENDDO
!      ENDDO
      icnt=icnt+1
      if(icnt.gt.20)debug=.false.
      return
      end
!
      subroutine fixup2(pall)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pall(3,natmax),pold(3,natmax),nequiv(natmax),p(3,natmax)
      real*8 x,y,z,x0,y0,z0,tt
      parameter (tol=0.5)
! Patch up the positions to remove symmetry issues
! Note: two cycles are probably not needed, but why not?
      do i=1,index
        do j=1,3
        pold(j,i)=pall(j,i)
        enddo
      enddo
      write(6,1) 1
1     format('Phase ',i1,', applying symmetry operations to atom')
      DO ITEST=1,2
      DO JATOM=1,index
         NTHIS=0
         X00=pall(1,JATOM)
         Y00=pall(2,JATOM)
         Z00=pall(3,JATOM)
         X0=0
         Y0=0
         Z0=0
            DO 10 I=1,IORD
            do 10 I1=-2,2
            DO 10 I2=-2,2
            DO 10 I3=-2,2
!       Generate equivalent sites
            X=TAU(1,I)+I1
            Y=TAU(2,I)+I2
            Z=TAU(3,I)+I3
            DO J=1,3
                X=X+dble(IZ(J,1,I))*pall(J,JATOM)
                Y=Y+dble(IZ(J,2,I))*pall(J,JATOM)
                Z=Z+dble(IZ(J,3,I))*pall(J,JATOM)
            ENDDO
!       How far is it, in a.u.
            X1=ABS(X-X00)*alat(1)
            Y1=ABS(Y-Y00)*alat(2)
            Z1=ABS(Z-Z00)*alat(3)
!
            IF((X1.LT.tol).AND.(Y1.LT.tol).AND.(Z1.LT.tol)) THEN
!               It is close, so add to average
                NTHIS=NTHIS+1
                X0=X0+X
                Y0=Y0+Y
                Z0=Z0+Z
            ENDIF
10          CONTINUE
       if(nthis.gt.0)then
!       Final averaged value
        TT=1.D0/dble(nthis)
        pall(1,JATOM)=X0*TT
        pall(2,JATOM)=Y0*TT
        pall(3,JATOM)=Z0*TT
       endif
       nequiv(jatom)=nthis
      ENDDO
      ENDDO
      do j=1,index
                write(6,88)j,nequiv(j),(pall(k,j),k=1,3),(pall(k,j)-pold(k,j),k=1,3)
                do jj=1,3
                        p(jj,j)=pall(jj,j)
                enddo
      enddo
      write(6,1)2
      IND=0
      DO JATOM=1,NAT
!	 First member of this atom
	 IF(JATOM .gt. 1)then
		IND=IND+MULT(JATOM-1)
	 else
		IND=1
	 endif
!	 Last
	 IUP=IND+MULT(JATOM)-1
	 DO M1=IND,IUP

         NTHIS=0
         X00=pall(1,M1)
         Y00=pall(2,M1)
         Z00=pall(3,M1)
         X0=0
         Y0=0
         Z0=0
 	  DO M2=IND,IUP
            DO 20 I=1,IORD
            DO 20 I1=-2,2
            DO 20 I2=-2,2
            DO 20 I3=-2,2
!       Generate equivalent sites
            X=TAU(1,I)+I1
            Y=TAU(2,I)+I2
            Z=TAU(3,I)+I3
            DO J=1,3
                X=X+dble(IZ(J,1,I))*p(J,M2)
                Y=Y+dble(IZ(J,2,I))*p(J,M2)
                Z=Z+dble(IZ(J,3,I))*p(J,M2)
            ENDDO
!       How far is it, in a.u.
            X1=ABS(X-X00)*alat(1)
            Y1=ABS(Y-Y00)*alat(2)
            Z1=ABS(Z-Z00)*alat(3)

            IF((X1.LT.tol).AND.(Y1.LT.tol).AND.(Z1.LT.tol)) THEN
!               It is close, so add to average
                NTHIS=NTHIS+1
                X0=X0+X
                Y0=Y0+Y
                Z0=Z0+Z
            ENDIF
20          CONTINUE
            ENDDO

	  Nequiv(M1)=nthis
          if(nthis.gt.0)then
!           Final averaged value
            TT=1.D0/dble(nthis)
            pall(1,M1)=X0*TT
            pall(2,M1)=Y0*TT
            pall(3,M1)=Z0*TT
	  ENDIF
	ENDDO
      ENDDO
      do j=1,index
                write(6,88)j,nequiv(j),(pall(k,j),k=1,3),(pall(k,j)-pold(k,j),k=1,3)
                if(pall(1,J) .lt. 0.D0)pall(1,j)=pall(1,j)+1.D0
                if(pall(2,J) .lt. 0.D0)pall(2,j)=pall(2,j)+1.D0
                if(pall(3,J) .lt. 0.D0)pall(3,j)=pall(3,j)+1.D0
                if(pall(1,J) .gt. 1.D0)pall(1,j)=pall(1,j)-1.D0
                if(pall(2,J) .gt. 1.D0)pall(2,j)=pall(2,j)-1.D0
                if(pall(3,J) .gt. 1.D0)pall(3,j)=pall(3,j)-1.D0
      enddo
      return
 88    format('Site ',i3,' Ops ',i2,' New: ',3f13.8,' Diff: ',3D12.5)
      end
!
      subroutine fixup3(pred)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax),pold(3,natmax)
      real*8 X,Y,Z,X0,Y0,Z0,X00,Y00,Z00,tmp(3),TT
      parameter (tol=0.25)
! Patch up the positions to remove symmetry issues
!
      pold(1:3,1:nat)=pred(1:3,1:nat)
      IND=0
      DO M1=1,NAT
         NTHIS=0
         X00=pred(1,M1)
         Y00=pred(2,M1)
         Z00=pred(3,M1)
         X0=0
         Y0=0
         Z0=0
!         Apply symmetry operations to generate temporary sites
 	  DO I2=1,IORD
            X=TAU(1,I2)
            Y=TAU(2,I2)
            Z=TAU(3,I2)
            DO J=1,3
                X=X+dble(IZ(J,1,I2))*pred(J,M1)
                Y=Y+dble(IZ(J,2,I2))*pred(J,M1)
                Z=Z+dble(IZ(J,3,I2))*pred(J,M1)
            ENDDO
            tmp(1)=X
            tmp(2)=Y
            tmp(3)=Z
!           Reapply symmetry operations
            DO 10 I=1,IORD
!       Generate equivalent sites
            DO 10 IT1=-2,2
            DO 10 IT2=-2,2
            DO 10 IT3=-2,2
            X=TAU(1,I)+IT1
            Y=TAU(2,I)+IT2
            Z=TAU(3,I)+IT3
            DO J=1,3
                X=X+dble(IZ(J,1,I))*tmp(j)
                Y=Y+dble(IZ(J,2,I))*tmp(j)
                Z=Z+dble(IZ(J,3,I))*tmp(j)
            ENDDO
!       How far is it, in a.u.
            X1=ABS(X-X00)*alat(1)
            Y1=ABS(Y-Y00)*alat(2)
            Z1=ABS(Z-Z00)*alat(3)

            IF((X1.LT.tol).AND.(Y1.LT.tol).AND.(Z1.LT.tol)) THEN
!               It is close, so add to average
                NTHIS=NTHIS+1
                X0=X0+X
                Y0=Y0+Y
                Z0=Z0+Z
!                write(78,79)nthis,X,Y,Z
            ENDIF
10          CONTINUE
          ENDDO
          if(nthis.gt.0)then
!           Final averaged value
            TT=1.D0/dble(nthis)
            pred(1,M1)=X0*TT
            pred(2,M1)=Y0*TT
            pred(3,M1)=Z0*TT
!            write(78,78)M1,nthis,(pred(L,M1),pold(L,M1),L=1,3)
78          format(2i4,6D15.8)
79          format(i4,3D16.9)
	  ENDIF
      ENDDO
!
      return
!
      end

