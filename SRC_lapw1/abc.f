      SUBROUTINE ABC (JATOM,J,jlo,lapw)
!
      use loabc, only : ALO, BLO, CLO, ELO, PLO, DPLO, PELO, &
                        DPELO, PEILO, PI12LO, PE12LO
      use atspdt, only  : INS, NL, LMMAX, LQIND, LM, LQNS, DP, DPE, E, P, PE, PEI, VNS1, VNS2, VNS3, GFAC
      use struk, only : POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      logical            lapw
      INTEGER            J,  JATOM, JLO
!
!     ..................................................................
!
!        ABC calculates the cofficients A,B,C of the local orbitals (LO)
!
!     ..................................................................
!
!        Local Parameters
!
      DOUBLE PRECISION   CUTOFF
      PARAMETER          (CUTOFF = 200.0D+0)
!
!        Local Scalars
!
      DOUBLE PRECISION   XAC, XBC, ALONORM
!
!        Intrinsic Functions
!
      INTRINSIC          SQRT
!

      if(lapw) then
         XAC = PLO(J,JATOM)*DPE(J+1,JATOM) -  &
              DPLO(J,JATOM)*PE(J+1,JATOM)
         XAC = XAC*RMT(JATOM)*RMT(JATOM)
         XBC = PLO(J,JATOM)*DP(J+1,JATOM) -  &
              DPLO(J,JATOM)*P(J+1,JATOM)
         XBC = -XBC*RMT(JATOM)*RMT(JATOM)
         CLO(J,jlo,JATOM) = XAC*(XAC + 2.0D+0*PI12LO(J,JATOM)) + &
              XBC*(XBC*PEI(J+1,JATOM) + 2.0D+0*PE12LO(J,JATOM)) + &
              1.0D+0
         CLO(J,jlo,JATOM) = 1.0D0/SQRT(CLO(J,jlo,JATOM))
         CLO(J,jlo,JATOM) = MIN(CLO(J,jlo,JATOM),CUTOFF)
         ALO(J,jlo,JATOM) = CLO(J,jlo,JATOM)*XAC
         BLO(J,jlo,JATOM) = CLO(J,jlo,JATOM)*XBC
      else
!.....APW definitions
         if(jlo.eq.1) then
           alonorm=sqrt(1.d0+(P(J+1,JATOM)/PE(J+1,JATOM))**2*PEI(J+1,JATOM))
            ALO(J,jlo,JATOM) = 1.d0 /alonorm 
            BLO(J,jlo,JATOM) = -P(J+1,JATOM)/PE(J+1,JATOM)/alonorm
            CLO(J,jlo,JATOM) = 0.d0
         else 
           xbc=-P(J+1,JATOM)/PLO(J,JATOM)
           xac=sqrt(1+xbc**2+2*xbc*PI12LO(J,JATOM))
            ALO(J,jlo,JATOM) = 1.d0/xac
            BLO(J,jlo,JATOM) = 0.d0
            CLO(J,jlo,JATOM) = xbc/xac
        endif
      endif
      WRITE(6,6000)J,ALO(J,jlo,JATOM),BLO(J,jlo,JATOM),CLO(J,jlo,JATOM)
!
      RETURN
!
 6000 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
!
!        End of 'ABC'
!
      END
