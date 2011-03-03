!BOP
! !ROUTINE: EnergyXSpectrum
! !INTERFACE:
      SUBROUTINE EnergyXSpectrum (PP, QQ, EnergyOut, wi,ind)
! !USES:
      use constants
	  use dimension_constants
      use spectra_hyperfine
	  use qvectors
	  use energygrid
	  use radial_functions
	  use input
	  use leftovers
	  use struct, only : nat,nats
	  use program_control, only : RelatQ
! !INPUT/OUTPUT PARAMETERS:
!     PP, QQ:      Radial functions of the core state (Large and small component).
!     EnergyOut:   Energy in eV of the scattered electron when the 
!                  core electron comes at the FERMI Level.
!     ind :        in which column of sdlm, ctr, x, etc. to put the results.; 1-3
!     wi  :        a prefactor multiplied to all arrays (contains branching ratio)
! !DESCRIPTION:
!     Calculates the double differential scattering cross section and integrates it
!     with respect to impuls transfer, yielding energy differential (partial) spectra.
!     Selection rules are taken into account.
!     Summation over equivalent positions is performed.
!     The input specifies whether j=l+1/2 or j=l-1/2 edge is calculated.

! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!   Added spin variables February 2005 (Kevin Jorissen)
!EOP
	  implicit none
 
 !    input variables :
      integer,intent(in) :: ind
      real*8,intent(in) ::  PP(NRAD+41), QQ(NRAD+41)
      real*8,intent(in) ::  EnergyOut, wi
      
!     external functions:
      real*8, external ::      ThreeJSymbol
      complex*16, external ::  DDlmLM

!     local variables:
      INTEGER l, ll, m, mm, s, ss    ! summation variables from the decomposition of the final states
      INTEGER IPos, IAtom,indx,J
      COMPLEX*16 Factor1, Factor2, Factor3, Factor4 
      COMPLEX*16 YQ((lc+lmax+1)*(lc+lmax+1))   ! spherical harmonics :   YQ(l²+l+1+m) = Ylm(Q)
      INTEGER Lambda, LambdaPr, mu, mupr   ! summation variables from the Rayleigh expansion of the plane wave operator
      real*8 JQ,JQPr,ei,prefac
      integer jshift             ! for the second edge

      
	  if (ind.eq.2) then
	    jshift=jdecal
	  else
	    jshift=0
	  endif

      DO J = JEMIN, JEMAX-JSHIFT
      EI = ENE(J)
!     only calculate if E>EFermi!
	  IF(EI.GE.0) THEN
         CALL QMesh (EnergyOut-EI)  ! energyout contains split, if necessary.
         
         DO l = 0, lmax
            DO ll = 0, lmax
			if(UseThisL(l).and.UseThisL(ll)) then
               DO IPos = 1, NPos
                  DO IAtom = NEqAtDeb, NEqAt
                     
                     CALL YLM (QV(1,IAtom,IPos), lc+max(l,ll), YQ)
				     if(relatQ) then
				       prefac=(QLenVClas(IAtom,IPos)**2-((energy-energyout+ei)/hbarc)**2)**2
					 else
					   prefac=QLenVClas(IAtom,IPos)**4
					 endif

                     DO Lambda = ABS(LC-l),LC+l
                        IF (UseThisLambda(Lambda).and.MOD(LC+l+Lambda,2).EQ.0) THEN
                        CALL RINT14(CFEIN1,CFEIN2, &
                             Uell(1,J,l),UPell(1,J,l),PP,QQ, &
                             JQ, QLenV(IAtom,IPos), Lambda, natom)
                        Factor1 = DCMPLX(dble(8)*wi*gammaa0**2*PI*dble(-1)**(l+ll) &       ! 4 to 16, KJ
                             * DSQRT ( DBLE( (2*Lambda+1) * (2*l+1) * (2*ll+1) ) ) &
                             * WeightV(IPos) * (2*LC+1) /prefac &
                             * JQ * ThreeJSymbol(LC,Lambda,l,0,0,0) )
                        
                        DO LambdaPr = ABS(LC-ll),LC+ll
                           IF (UseThisLambda(LambdaPr).and.MOD(LC+ll+LambdaPr,2).EQ.0) THEN
                           CALL RINT14(CFEIN1,CFEIN2, Uell(1,J,ll),UPell(1,J,ll), &
                                   PP,QQ,JQPr, QLenV(IAtom,IPos),LambdaPr, natom)
                           Factor2 = DCMPLX(Factor1 * DSQRT(DBLE(2*LambdaPr+1)) &
                                * DCMPLX(0,1)**(Lambda-LambdaPr) &
                                * JQPr*ThreeJSymbol(LC,LambdaPr,ll,0,0,0))
                           
                           DO Mu = -Lambda, Lambda
                              DO m = -l, l
                                 Factor3 = DCMPLX(ThreeJSymbol(LC,Lambda,l,-Mu-m,Mu,m)  &
                                   * Factor2 * DCONJG(YQ(Lambda*Lambda+Lambda+1+Mu)) )

                                 DO MuPr = -LambdaPr, LambdaPr
                                    mm = Mu+m-MuPr
                                    IF((-ll.le.mm).and.(ll.ge.mm)) THEN

                                      Factor4=dcmplx(0,0)
									  DO s = 1,nspin
									    DO ss = 1,nspin
                                           Factor4 = Factor4 + DCMPLX(ThreeJSymbol(LC,LambdaPr,ll,-MuPr-mm,MuPr,mm) &
                                                * YQ(LambdaPr*LambdaPr+LambdaPr+1+MuPr) &
                                                * DDlmLM(J,l,m,s,ll,mm,ss)/UellUell(J,l,ll) ) * Factor3
									    ENDDO   ! ss
									  ENDDO  ! s
                                      X(J+JSHIFT,ind)  = X(J+JSHIFT,ind)  + Factor4

!-----here we evaluate explicitely the cross terms and partial DOS contributions
! only s p d final states
	  if ((l.LE.2).AND.(ll.LE.2)) then
	    if ((l.EQ.ll).AND.(m.EQ.mm)) then
!         we start with the DOS contributions
	      indx=1+l*(l+1)+m
	      SDLM(indx,J+JSHIFT,ind)=SDLM(indx,J+JSHIFT,ind) + dble (Factor4)
! dipole-allowed L mixing
	    elseif (ABS(l-ll).EQ.2) then
	      if (ABS(m-mm).EQ.2) then		  ! 00 22
	        CTR(1,J+JSHIFT,ind)=CTR(1,J+JSHIFT,ind) + Factor4
	      else if (ABS(m-mm).EQ.1) then		  ! 00 21
	        CTR(2,J+JSHIFT,ind)=CTR(2,J+JSHIFT,ind) + Factor4
	      else if (ABS(m-mm).EQ.0) then		  ! 00 20
	        CTR(3,J+JSHIFT,ind)=CTR(3,J+JSHIFT,ind) + Factor4
	      end if	! all possible m' M' combinations for l'=0,L'=2
! l'=L', m'M' mixing (not for l'=L'=0)
	    elseif ((l.EQ.ll).AND.(l.GT.0)) then
! l'=L'=1
	      if (l.EQ.1) then
	        if (ABS(m-mm).EQ.2) then		! 11 1-1
	          CTR(4,J+JSHIFT,ind)=CTR(4,J+JSHIFT,ind) + Factor4
	        else if (ABS(m-mm).EQ.1) then		! 11 10
	          CTR(5,J+JSHIFT,ind)=CTR(5,J+JSHIFT,ind) + Factor4
	        end if  	! all possible m'M' combinations for l'=L'=1
! l'=L'=2
	      elseif (l.EQ.2) then
	        if (ABS(m-mm).EQ.4) then		! 22 2-2
	          CTR(6,J+JSHIFT,ind)=CTR(6,J+JSHIFT,ind) + Factor4
	        else if (ABS(m-mm).EQ.3) then		! 22 2-1
	          CTR(7,J+JSHIFT,ind)=CTR(7,J+JSHIFT,ind) + Factor4
	        else if (ABS(m-mm).EQ.2) then		! 21 2-1 or 22 20
	          if (ABS(m+mm).EQ.2) then		  ! 22 20
	            CTR(8,J+JSHIFT,ind)=CTR(8,J+JSHIFT,ind) + Factor4
	          else					  ! 21 2-1
	            CTR(9,J+JSHIFT,ind)=CTR(9,J+JSHIFT,ind) + Factor4
	          end if
	        else if (ABS(m-mm).EQ.1) then		! 22 21 or 21 20
	          if (ABS(m+mm).EQ.3) then		  ! 22 21
	            CTR(10,J+JSHIFT,ind)=CTR(10,J+JSHIFT,ind) + Factor4
	          else					  ! 21 20
	            CTR(11,J+JSHIFT,ind)=CTR(11,J+JSHIFT,ind) + Factor4
	          end if
	        end if  	! all possible m'M' combinations for l'=L'=2
	      end if      ! l'=L'=2	
	    end if  ! what to do if l,ll < 3

	  elseif(l.eq.3.and.ll.eq.3.and.m.eq.mm) then
	    sdlm(10,J+JSHIFT,ind)=sdlm(10,J+JSHIFT,ind) + Factor4
	  else   !  term is a CT with l=3 or ll=3
	    ctr(12,j+JSHIFT,ind)=ctr(12,j+JSHIFT,ind) + Factor4
      end if      ! making partial spectra and cross terms
!---------------------------------------------------------------------------	
                                    ENDIF !mm
                                 ENDDO !MuPr
                              ENDDO !m
                           ENDDO !Mu
                           ENDIF ! LC+ll+LambdaPr must be even ; selection rule for lambdapr
                        ENDDO !LambdaPr
                        ENDIF ! LC+l+Lambda must be even ; selection rule for lambda
                     ENDDO !Lambda
                     
                  ENDDO !IAtom
               ENDDO !IPos
            endif ! selection rule UseThesel's
            ENDDO !ll
         ENDDO !l
         
      ENDIF
      ENDDO !J
      
      RETURN
      END

