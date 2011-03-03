	subroutine drmngwien(xin,n,iter,fret,nat,ende,fwert,mult,bond,step0)
!
!       Main calling routing for port3 code drmng
!       June 2004: Added "right" multiplicity and bond term
!       Seems to be stable
!
        implicit real*8 (a-h,o-z)
        dimension IV(60),xin(n,*),fwert(n,*),mult(*),bond(*)
	real*8,allocatable :: V(:),D(:),G(:),X(:)
        real*8,allocatable :: XOLD(:),GOLD(:),P(:)
        real*8,allocatable :: HESS(:), WTMP(:)
	logical ende,touch,hinit,init,first,getscale
        nok=0
        getscale=.true.
        first=.true.
        init=.true.
        f0=-0.5
	LIV=60
	LV=71+N*(N+13)/2+1
        LHESS=(N*(N+1))/2+4
	allocate (G(N),V(LV),D(N),X(N),HESS(LHESS))
        allocate (XOLD(N),GOLD(N),P(N),WTMP(N))
	ende=.false.
! Input values
	call DIVSET(2, IV, LIV, LV, V)
! D.... SCALE VECTOR, take as 1 (other values possible)
	do j=1,N
	   D(J)=1
	   X(J)=XIN(J,1)
           XOLD(J)=X(J)
           G(J)=0
	enddo
        do j=1,LHESS
                HESS(j)=0
        enddo
!       Initialize Hessian 
!       Is there a .minrestart file?
!       N.B., presumably inquire is standard
        inquire(file='.minrestart',exist=hinit)
        if(hinit)then
!       Yes, file is there so let's read it
        write(6,*)'Retrieving prior Hessian'
        open(unit=77,file='.minrestart')
        read(77,*)(HESS(J),J=1,(N*(N+1))/2)
        write(6,*) 'reading hessian from .minrestart file'
        close(unit=77)
        endif
!       These are some of the key arrays used by drmng
!       FX... FUNCTION VALUE (p).
!       G.... GRADIENT VECTOR.
!       IV... INTEGER VALUE ARRAY.
!       LIV.. LENGTH OF IV (AT LEAST 60).
!       LV... LENGTH OF V (AT LEAST 71 + N*(N+13)/2).
!       N.... NUMBER OF VARIABLES (COMPONENTS IN X AND G).
!       V.... FLOATING-POINT VALUE ARRAY.
!       X.... VECTOR OF PARAMETERS TO BE OPTIMIZED.
!
!       Set some default values (see dmng.f for details)
	iv(1)=12
	iv(25)=1
	iv(17)=200
	iv(18)=150
	iv(19)=1
	iv(20)=1
	iv(21)=6
	iv(22)=1
        itport=1
!	X convergence tolerance, needs to be > eps in position test
        V(33)=2d-6
!	2D norm for first step, settable by user, default 0.25
	V(35)=step0
!       Pass stored Hessian ?
        if(hinit)then
!               Have drmng initialize itself
                IV(1)=13
                call DRMNG(D, FX, G, IV, LIV, LV, N, V, X, HESS)
!               Store Cholesky factor
                L=iv(42)
		do j=1,(N*(N+1))/2
                        V(L+J-1)=HESS(J)
                enddo
                IV(25)=0
        else
!       No user Hessian provided
!       Code multiplicity into initial hessian
                I=0
                tmp=1
                do j=1,N,3
                        I=I+1
!       Include here the multiplicity and user-bond term
                        T=abs(mult(I))*bond(I)
                        P(J)=T
                        P(J+1)=T
                        P(J+2)=T
                        tmp=tmp*T*T*T
                enddo
!       Normalization for matrix, to unity (about right)
                tmp=(tmp**(1.D0/dble(N)))
                do j=1,N
!       sqrt -- Cholesky factor, 3 as bond strength estimate about OK
!       The "3.d0" could be changed to something smaller for soft systems
!       We cannot satisfy eveyone!
                        P(J)=sqrt(3.d0*P(J)/tmp)
                enddo
!       Copy diagonal elements to local copy of Hessian
                I=0
                DO J=1,N
                        I=I+J
                        HESS(I)=P(J)
                ENDDO
!       Let drmng initialize itself
                IV(1)=13
                call DRMNG(D, FX, G, IV, LIV, LV, N, V, X, HESS)
!       Code in Hessian
                L=iv(42)
                do j=1,(N*(N+1))/2
                        V(L+J-1)=HESS(J)
                enddo
!                WRITE(6,*)'Enter ',hess(1),V(L)
                IV(25)=0
        endif
!       These V values might be changeable
!       DECFAC -- how much to reduce radius     (default 0.5)
!       INCFAC -- how much to increase          (default 2.0)
!       V(INCFAC)=1.5 might be better
!       PARAMETER (DECFAC=22, INCFAC=23)
        write(6,*)':MIN Minimizing for ',n,' parameters'
10	continue
!       Save prior values (might be useful?)
!        write(6,*)'Saving old'
!        DO J=1,N
!           GOLD(J)=G(J)
!           XOLD(J)=X(J)
!        ENDDO

!	Call port routine
20      continue
!       write(6,*)':DRMNG Calling drmng',itport
!       Note: F0 is used to offset energies for "nice" E printout via :MIN
!       Changes relative to -0.5 are easier to see, but should not effect
!       the code performance
        FX=FX-F0
        L=iv(42)
!        write(6,*)'Deb into ',L,V(L)
!       Let drmng do everything for us
 	call DRMNG(D, FX, G, IV, LIV, LV, N, V, X, HESS)
        itport=itport+1	
!
!       We need to take action depending upon the return code of drmng
	if(iv(1).eq.1)then
!       IV(1) = 1 MEANS THE CALLER SHOULD SET FX TO F(X), THE FUNCTION VALUE
!             AT X, AND CALL DRMNG AGAIN, HAVING CHANGED NONE OF THE
!             OTHER PARAMETERS.  AN EXCEPTION OCCURS IF F(X) CANNOT BE
!             (E.G. IF OVERFLOW WOULD OCCUR), WHICH MAY HAPPEN BECAUSE
!             OF AN OVERSIZED STEP.  IN THIS CASE THE CALLER SHOULD SET
!             IV(TOOBIG) = IV(2) TO 1, WHICH WILL CAUSE DRMNG TO IG-
!             NORE FX AND TRY A SMALLER STEP.  THE PARAMETER NF THAT
!             DMNG PASSES TO CALCF (FOR POSSIBLE USE BY CALCG) IS A
!             COPY OF IV(NFCALL) = IV(6).
!
!       This should be trapped, and with the code going to 20
!
601     format(3f12.7,'  ',3f12.6)
!       Is there a stored energy/gradient for this position?
      	   fx=func2(x)
!       Setup F0 if first entry
           if(first)then
                f0=min(-0.5d0,fx+0.5d0)
                first=.false.
           endif
!	If it is not there, exit for a reverse call
      	   if (fx.eq.999.111d0) then
!       Check that the point is viable in terms of overlapping sphere
           do J=1,N
                P(J)=0
           enddo
           call traptouch(x,p,n,stpmax,nat,touch)
!!!              touch=.false.
           if(touch)then
!       Atoms touch
               write(6,*)':WARNING: Step size reduced due to ', &
                'overlapping spheres -- check RMT'
!       Code in flag that tell drmng step was unfeasible, and let it
!       find a new point
               IV(2)=1
               goto 20
           endif
           fret=fx
!       Store the Hessian in .min_hess for a possible restart
        open(unit=77,file='.min_hess')
        write(77,*)(HESS(J),J=1,(N*(N+1))/2)
        close(unit=77)
        I=0
!        DO J=1,N
!        I=I+J
!        GOLD(J)=HESS(I)
!        ENDDO
!       Evaluate diagonal terms of Hessian for user display
!       Most people won't know what to do with this...
        I=0
        do j1=1,n
                T=0
                do j2=1,j1
                        I=I+1
                        T=T+HESS(I)*HESS(I)
                enddo
                GOLD(j1)=T
        enddo
        write(6,*)'      X           Y           Z     ',  & 
        '        HX          HY          HZ'
           do j=1,n,3
               write(6,601)x(j),x(j+1),x(j+2),gold(j),gold(j+1),gold(j+2)
           enddo
           write(6,*)':ENE Returning for function',itport
           return
        endif
        write(6,*)':ENE Retrived energy',fx
        niter=niter+1
        write(6,*)'      X           Y           Z     ',  & 
        '        HX          HY          HZ'
!       Evaluate diagonal terms of Hessian
        I=0
        do j1=1,n
                T=0
                do j2=1,j1
                        I=I+1
                        T=T+HESS(I)*HESS(I)
                enddo
                GOLD(j1)=T
        enddo
!       Get the gradient as well, and check mod & RMS for monitoring
        GSUM=0
        GSUM2=0
        call dfunc2(x,g)
        do j=1,n,3
           write(6,601)x(j),x(j+1),x(j+2),gold(j),gold(j+1),gold(j+2)
             GSUM=GSUM+abs(G(j))+abs(G(J+1))+abs(G(j+2))
             GSUM2=GSUM2+G(j)**2 + G(j+1)**2 + G(j+2)**2
        enddo
        write(6,*)':GRAD Abs mean ',GSUM/n,' RMS ',sqrt(GSUM2/n)

        goto 10
	else if(iv(1).eq.2)then
!       IV(1) = 2 MEANS THE CALLER SHOULD SET G TO G(X), THE GRADIENT VECTOR
!             OF F AT X, AND CALL DRMNG AGAIN, HAVING CHANGED NONE OF
!             THE OTHER PARAMETERS EXCEPT POSSIBLY THE SCALE VECTOR D
!             WHEN IV(DTYPE) = 0.  THE PARAMETER NF THAT  DMNG PASSES
!             TO CALCG IS IV(NFGCAL) = IV(7).  IF G(X) CANNOT BE
!             EVALUATED, THEN THE CALLER MAY SET IV(TOOBIG) TO 0, IN
!             WHICH CASE DRMNG WILL RETURN WITH IV(1) = 65.
!
!       This should be trapped, and with the code going to 20
!
      	   fp=func2(x)
      	   if (fp.eq.999.111d0) then
!       Check that the point is viable
           do J=1,N
               P(J)=0
           enddo
           call traptouch(x,p,n,stpmax,nat,touch)
!!!        touch=.false.
           if(touch)then
!       Atoms touch
               IV(2)=0
               goto 20
           endif
           fret=fp
           write(6,*)':GRAD Returning for gradient',itport
!       We only get here if energy exists but not gradient (makes no sense)
!       Never seen the code get here, but anyway
           write(6,*)':WARNING: This should not happen!'
           return
        endif
!        write(6,*)':GRAD Retrived function, getting grad'
      	call dfunc2(x,g)
        write(6,*)'      X           Y           Z     ',  & 
        '       GX           GY          GZ'
!        GSUM=0
        do j=1,n,3
             write(6,601)x(j),x(j+1),x(j+2), &
             g(j),g(j+1),g(j+2)
!             GSUM=GSUM+abs(G(j))+abs(G(J+1))+abs(G(j+2))
        enddo
!        write(6,*)':GRAD Abs mean ',GSUM/n
!       Scaling ?
        nok=nok+1
        if(getscale.and.nok.gt.1)then
        GDIFF=0
        XDIFF=0
        GX=0
        SHS=0
!     Have a look at Hessian scaling, not very useful at the moment
!     Evaluate S*H (I think)
      I0 = 0
!      DO I = 1, N
!         YI = X(I)-XOLD(I)
!         WTMP(I) = ZERO
!         DO J = 1, I
!              IJ = I0 + J
!              WTMP(J) = WTMP(J) + YI*HESS(IJ)
!         ENDDO
!         I0 = I0 + I
!      ENDDO
!        DO j=1,n
!                t1=G(J)-GOLD(J)
!                t2=X(J)-XOLD(J)
!                GDIFF=GDIFF+t1*t1
!                XDIFF=XDIFF+t2*t2
!                GX=GX+t1*t2
!                SHS=SHS+WTMP(J)*WTMP(J)
!        enddo
!       Oren-Spedicato diagonal term (inverse form)
!        DOS=abs(GX/SHS)
!        write(6,*)':Plausible Hessian scaling ',DOS
!        getscale=.false.
        endif
!       Save prior values (might be useful?)
!        write(6,*)'Saving old'
        DO J=1,N
           GOLD(J)=G(J)
           XOLD(J)=X(J)
        ENDDO
	goto 10
!       These are termination conditions
!       In most cases the gradient term in case.inM will determine things
	else if(iv(1).eq.3)then
	   write(6,*)':MIN Position convergence achieved'
           CALL OUTERR('MINI','Position convergence achieved')
	else if(iv(1).eq.4)then
	   write(6,*)':MIN Relative energy convergence'
           CALL OUTERR('MINI','Relative energy convergence')
	else if(iv(1).eq.5)then
	   write(6,*)':MIN Position & energy convergence'
           CALL OUTERR('MINI','Position & energy convergence')
	else if(iv(1).eq.6)then
	   write(6,*)':MIN Absolute energy convergence'
!       We should not really get here
           write(6,*)':WARNING Please check that forces have converged'
           CALL OUTERR('MINI','Absolute energy convergence')
	else if(iv(1).eq.7.or.iv(1).eq.8)then
	   write(6,*)':MIN Probably converged to numerical accuracy'
           write(6,*)':WARNING you might have inconsistent Forces/Energies'
           CALL OUTERR('MINI','Probably convergenced -- check log')
	else if(iv(1).eq.9.or.iv(1).eq.10)then
	   write(6,*)':MIN Function evaluation limit reached'
           write(6,*)':WARNING Please restart anew from final position'
           CALL OUTERR('MINI','Function evaluation limit reached')
	else 
	   write(6,*)':MIN Something obscure happened'
           CALL OUTERR('MINI','Unknown error -- check case.outputM')
	   write(6,*)':MIN drmng return code',iv(1)
	endif
        fret=fx
	ende=.true.
!       Output final information
        write(6,*)':MIN BFGS minimizer drmng terminating'
        write(6,*)':MIN Final function value',V(10)
        IG=iv(28)-1
        write(6,*)'      X           Y           Z     ',  & 
        '       GX           GY          GZ'
!       The v values look wrong (not sure why)
        do j=1,n,3
           write(6,601)x(j),x(j+1),x(j+2),  &
             v(ig+j),v(ig+j+1),v(ig+j+2) 
        enddo
!       Maybe right ?
        IOUT=3*n*(iter)
!         write(*,*) 'iout1',iout,n,iter,ig
        do j=1,n
!           xin(j+iout)=x(j)
           xin(j,iter+1)=x(j)
           fwert(j,iter+1)=v(ig+j)
        enddo
	return
	end
