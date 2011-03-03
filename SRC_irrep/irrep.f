      PROGRAM IRREDUCIBLE_REPRESENTATIONS
      USE FELDER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!...............................................................................
!     Determines the irreducible representations (IR) of the eigenstates at a 
!     given k-state; both single-valued and the double-valued space groups are 
!     implemented. The IRs of the space groups are presented in terms of the
!     'relevant' IRs of the point groups. For some special k-points at the 
!     Brillouin zone surface in non-symmorphic crystals, the relevant IRs of 
!     the factor group G(k)/T(k) are presented, see Refs. 1-3. The program 
!     checks for additional denegeracy due to time-reversal symmetry. The code 
!     requires plane wave based basis functions, and it has been developed for 
!     the WIEN2k full-potential LAPW package.
!
!     Further information can be found in:
!     [1] J.F. Cornwell, Group Theory in Physics, 
!         Vol. 1, (Academic Press, London, 1984).
!     [2] J.F. Cornwell, Group Theory and Electronic Energy Bands in Solids, 
!         Selected Topics in Solid State Physics, editor: E.P. Wohlfarth, 
!         Vol. 10 (North-Holland, Amsterdam, 1969).
!     [3] H.-W. Streitwolf, Group Theory in Solid-State Physics,  
!         (Macdonald, London, 1971).
!     [4] W. Jones and N.H. March, Theoretical Solid State Physics, Vol. 2, 
!         Non-Equilibrium and Disorder (Wiley, New York, 1973).
!     [5] C.J. Bradley and A.P. Cracknell, The Mathematical Theory of Symmetry 
!         in Solids. Representation Theory for Point Groups and Space Groups,
!         (Clarendon Press, Oxford, 1972).
!
!
!     The computational method is according to:
!     [6] C. Persson, Electronic Structure of Intrinsic and Doped Silicon 
!         Carbide and Si, PhD thesis: ISBN 91-7219-442-1, Linkoping University, 
!         Sweden (1999).
!
!     The notation of the irreducible representations is from:
!     [7] G.F. Koster, J.O. Dimmock, R.G. Wheeler, and H. Statz, Properties of 
!         the Thirty-Two Point Groups (MIT, Cambridge Massachusetts, 1963).
!     [8] Schoenflies notation; See for instance S.L. Altmann and P. Herzig, 
!         Point-Group Theory Tables (Oxford University, Oxford, 1994).
!
!     Clas Persson,    Condensed Matter Theory Group, 
!                      Department of Physics, 
!                      Box 530, University of Uppsala, 
!                      SE-751 21 Uppsala, Sweden
!                      Email: Clas.Persson@fysik.uu.se
!
!     Uppsala, Sweden, September  2000
!     last updated:    October    2000   C. Persson
!                      November   2001   P. Blaha, F77 -> F90
!                      December   2001   C. Persson, tested for WIEN2k
!                      April      2002   C. Persson, output for band drawing
!...............................................................................
!
!.....input files:
!     *.vector      -  unit=10,9, eigenvalues and eigenfunctions 
!     *.struct      -  unit=20,   symmetry operations Pi={Ri|ti}
!
!.....output files:
!     *.irrrep      -  unit=5,    data file for band plotting (spaghetti)
!     *.outputir    -  unit=6,    irreducible representations (IRs).
!
!.....notation: 
!     Pi={Ri|ti}    -  crystallographic symmetry operations {rotation|transl};
!                      ti=taui+tm; 0<=taui<1 and tm is a lattice vectors.
!                      The product Pi*Pj = {Ri*Rj|Ri*tj+ti} 
!     inv(Pi)       -  inverse of Pi; inv(Pi) = {inv(Ri)|-inv(Ri)*ti}
!     Ri~           -  transpose of Ri; Ri~ not always inv(Ri) 
!     G             -  space group, consisting of the elements {Ri|ti} 
!     G(k)          -  space group of the allowed wave vector k 
!                      (little space group); those {Ri|ti} such that  
!                      k*inv(Ri)=k+K, where K is a reciprocal lattice vector.
!     Go            -  crystallographic point group of G, consisting of the 
!                      elements {Ri|0} 
!     Go(k)         -  point group of the allowed wave vector k
!     T(k)          -  translation group of allowed wave vector k; 
!                      those {1|tm} such  that exp(-i*k*tm)=1; 
!                      T(k) is an invariant (normal) subgroup to G(k)
!     F(k)          -  factor group G(k)/T(k), used for special k-points at 
!                      BZ surface in non-symmorphic crystals.
!     Fo(k)         -  point group of F(k)
!     IRp           -  p:th 'relevant' irreducible representation of 
!                      Go(k) or Fo(k). IRp is real if all representation 
!                      matrices are real. IRp can be equivalent to a real 
!                      representation (IRp ~ IRq_real), or equivalent to its 
!                      complex conjugate (IRp ~ (IRp)*), or essentially 
!                      complex (IRp !~ (IRp)* !~ IRq_real)
!                      
!
!.....basic theory, Ref. [5] p46-48:
!     <psi_l(k,r)|Pi*psi_l'(k,r)> = exp(-ik*ti)*IRp({Ri|0})_l,l' for k-points 
!                      obeying the Cornwell conditions (Ref [1] p239)
!     psi(k,r)      -  nf*SUM{ c(k+K)*fi(k+K,r) };  fi(k,r)=exp(ikr); 
!                      nf is the normalization factor
!     Pi*psi(k,r)   -  nf*SUM{ c(k+K)*exp( i(k+K)inv(Pi)r]) }
!                     =nf*SUM{ c(k+K)*exp( i(k+K)inv(Ri)(r-ti)) }
!                     =nf*SUM{ c(k+K)*exp(-i(k+K)inv(Ri)ti)*fi((k+K)inv(Ri),r)}
!                     =nf*SUM{ c((k+K')Ri) *exp(-i(k+K')ti) * fi(k+K',r) },
!                      where (k+K)inv(Ri)=k+K' => k+K=(k+G')Ri 
!                      In Ref. [6] p47 orthogonal matrices are presumed,
!                      i.e, Ri~ =inv(Ri)
!
!.....input:
!     FL(4)         -  flags:
!                      (1) true if complex eigenfunctions
!                      (2) true if spinors (i.e., with spin-orbit coupling)
!                      (3) true if crystallographic inversion symmetry  
!                      (4) true if spin-polarized
!     IZ,TAU,IIZ    -  symm. ops. {Ri|ti} read from case.struct. IIZ=inv(IZ)
!     A(*,*)        -  spin-up   part of eigenfunctions 
!     B(*,*)        -  spin-down part of eigenfunctions
!     SK=ISK/ISKDEN -  k-point; ISK and ISKDEN are integers
!     EN(NE)        -  eigenvalues in Rydberg  
!     KV(3,NV)      -  (k+K)*ISKDEN, where K is reciprocal lattice vector
! 
!...............................................................................
!.....transformations:
!     BR1(3,3)      -  lapw1 k-list coord. into Cartesian coord (recipr. space)
!     BR2(3,3)      -  reciprocal coord.   into Cartesian coord (recipr. space)
!     DR1,DR2       -  inv(BR1), inv(BR2)
!     DB1           -  DR2*BR1; from lapw1 k-list coord. into primitiv coord.
!
!.....calculated quantities:
!     RNAME,CNAME   -  name of rotaions and classes
!     RAN(4,*)      -  rotation direction of Ri
!     IAXC2(3)      -  main C2 axes
!     SU2(2,2,*)    -  spin rotation of double groups
!     JIJ(*,*)      -  equals n in Rn=Rj*Ri*inv(Rj)
!     PH(*,*)       -  phase factor  
!     LKG(IKG)      -  list of symm.op. in the little k-group G(k)
!     FGT           -  true if IR of G(k) cannot simple be associated to an IR
!                      of Go(k), which can happen for non-symmorphic crystals 
!                      if k is a special k-point at BZ surface. In this case 
!                      the local relevant IR is presented.
!     GAM(I,J)      -  IR matrix element (I,J) of Go(k) for a certain energy 
!                      and symm.op. I,J=1,...ND, where ND is the no of 
!                      degenerate states
!
!.....output:
!     XM(IKG,NE)    -  Characters for each eigenstate and each symm. ops.
!                      XM(*,*)=trace(GAM)
!
!.....common block:
!     COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
!     CTIR(i,j)     -  Character tables of the 32 point group
!                      String format for each irreducible representation
!     TTIR          -  the title of the classes
!     ZTIR          -  same as CTIR, but complex numbers instead of string
!     NTAB( 1)      -  Table number in G.F. Koster, et al., Ref. 7
!     NTAB( 2)      -  Page number for the table in Ref. 7
!     NTAB( 3)      -  Total no. of IRs 
!     NTAB( 4)      -  No. of IRs in the single group
!     JTAB( 1)      -  Table number in S.L. Altmann, et al., Ref. 8
!     JTAB( 2)      -  Page number for the table in Ref. 8
!     
!*******************************************************************
      COMPLEX*16       SU2(2,2,NSYM)
      CHARACTER*80     DEFFN,ERRFN
      CHARACTER*6      RNAM(NSYM),CNAM(NSYM)
      CHARACTER*3      GRPNAM
      LOGICAL          FGT
      LOGICAL          FL(FLMAX)
      DIMENSION        SK(3),ISK(3),IAXC2(3)
      DIMENSION        IZ(3,3,NSYM),IIZ(3,3,NSYM),TAU(3,NSYM)
      DIMENSION        LKG(NSYM),JIJ(NSYM,NSYM)
      DIMENSION        RAN(4,NSYM),ILC(NSYM),NROT(NSYM)
      DIMENSION        BR1(3,3),BR2(3,3)
      DIMENSION        DR1(3,3),DB1(3,3)
!
!*******************************************************************
!
      CALL OPNFS(DEFFN,ERRFN,FL)
      CALL ERRFLG(ERRFN,'ERROR IN IRREP')
      CALL CPUTIM(DTIME1)	
!
!.....initialization: read struct-file and define transformations
      CALL INIT(FL,NTYPE,NLO,IZ,IIZ,TAU,IORD,BR1,BR2,DR1,DB1)
!     
!.....properties of Ri: Refs [1] p25,55; [2] p10,197; [3] p16,169
      CALL RMPROP(IZ,IIZ,TAU,IORD,BR1,DR1,SU2,JIJ,ILC,RAN,RNAM)
!
!.....loop over k-points
      KKK=0
  100 KKK=KKK+1
!
!.....read the vector file.
      CALL KPTIN(FL,NTYPE,NLO,SK,NE,NV,ISK,ISKDEN,KKK)
      IF(KKK.LE.0) GOTO 101
!
!.....find elements of G(k): Refs [1] p235; [2] p89; [3] p79
      CALL KGROUP(FL,DB1,IZ,IIZ,IORD,ISK,ISKDEN,LKG,IKG,KKK,NV)
!     
!.....identify the point group Go(k)
      CALL PNTGRP(FL,LKG,IKG,ILC,JIJ,RNAM,GRPNAM,NCC,NROT)
!
!.....check the 'Cornwell condition': Ref [1] p239
      CALL CRWCND(ISK,ISKDEN,IZ,TAU,LKG,IKG,FGT)
!
!.....characterize according to Fo(k): Ref [1] p241; [2] p168
      IF(.NOT.FGT) THEN
!       THIS WILL BE IMPLEMENTED SOON     
!.......find the factor group F(k)=G(k)/T(k)
        CALL MDFPG(FL,FGT)
!.......character table of Fo(k)
!       CALL MPNTGRP()
!.......time-reversal symm: Ref [1] p158; [2] p149,210; [3] p177
        CALL TRSYMB(FL,IZ,ISK,ISKDEN,LKG,IKG,TAU,SK,IORD,GRPNAM,RAN,KKK)             
!
      ELSE
!.......determine classes
!        call flush(6)
        CALL CLASSE(FL,IZ,LKG,IKG,ILC,NROT,JIJ,RAN,RNAM,CNAM,GRPNAM,IAXC2)
!
!.......time-reversal symm: Ref [1] p158; [2] p149,210; [3] p177
        CALL TRSYMA(FL,ISK,ISKDEN,LKG,IKG,TAU,CNAM,SK,IORD,GRPNAM,RAN,IAXC2,KKK)
!
!.......rotate reciprocal lattice vector.
        CALL ROTKV(LKG,IKG,NV,ISK,ISKDEN,IZ,IIZ,TAU) 
!      
!.......calculate the characters: Ref [6] p46
        CALL CHRCT(FL,NE,NV,LKG,IKG,IORD,SU2,KKK)
      ENDIF
!
!.....identify the irreducible representations
      CALL WRTIR(FL,FGT,NE,LKG,IKG,CNAM,KKK)
!      
      DEALLOCATE(A,B,EE,KV,L,PH,XM)
      GOTO 100
  101 CALL CPUTIM(DTIME4)
      WRITE(6,500) 'TOTAL TIME:', DTIME4-DTIME1
!
      CALL ERRCLR(ERRFN)
      STOP 'IRREP END'
 500  FORMAT(/,80('*'),//,7X,A11,F7.2)
      END
      


