      SUBROUTINE PAROP(iat,nl,L,OP,ipr)
      use opr
!      use ldau
      use struct
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'param.inc'
      logical rel
! calculation of Slater parameters and orbital polarization parameter
      DIMENSION RHO(NRAD),VALUE(NRAD),vr(nrad)
      DIMENSION R(NRAD),F(nrad,0:6),d(3,0:6),ff(0:6)
      COMMON /UHELP/   A(NRAD),B(NRAD),AP(NRAD),BP(NRAD),AE(NRAD), &
                       BE(NRAD)
       jat=iatom(iat)
       jri=jrj(jat)
       dx=dxm(jat)
       rnot=ro(jat)
       z=zz(jat)
!.....set up radial mesh.
      DO 10 I=1,JRI
   10 R(I)=RMT(JAT)*DEXP(DFLOAT(I-JRI)*DX)
! calculate Slater and Racah parameters ab initio
! coefficients connecting F^k with F_k (Griffith:The Theory of Transition
! Metal Ions Table 4.4, p.77;  F_k=F^k/d_k**2
! d -electrons
      d(2,0)=1.d0
      d(2,2)=7.d0
      d(2,4)=21.d0
! f-electrons
      d(3,0)=1.d0
      d(3,2)=15.d0
      d(3,4)=33.d0
      d(3,6)=429.d0/5.d0
!
      REL=(.true.)
      CFEIN1 = 1.0D+0
      CFEIN2 = 1.d0/137.0359895d0**2
          fl=L
! Eop is in Ry, E in H
            e=eop(nl,iat)/2.d0
              do 2 i=1,jri
               vr(i)=vsp(i,jat)/2.d0
2             continue
! calculate the wave function at energy E
            call outwin(rel,vr,rnot,dx,jri,E,fl,uv,duv,nodel,z)
!
! renormalize the density within muffin tin sphere
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,IAT,RNOT,DX,JRI)
            if(ipr.gt.1)write(6,555)iat,l,E,ovlp
555   format(' iat, L, E, ovlp',2i3,2e14.6)
              do 3 i=1,jri
               rho(i)=(a(i)**2+cfein2*b(i)**2)/ovlp
3             continue
!
! calculate Slater integrals F^k using the integration from LAPW0
! including F^0
        do 20 l1=0,l
          k=2*l1
            do 21 i=1,jri
! r' < r
              do 22 j=1,i-1
                value(j)=rho(j)*r(j)**k/r(i)**(k+1)
22            continue
! r' > r
              do 23 j=i,jri
                value(j)=rho(j)*r(i)**k/r(j)**(k+1)
23            continue
             call charg2(r,dx,value,1,i,b1)
             call charg3(r,dx,value,i,jri,b2)
             Fk=value(1)*r(1)/2.d0+b1+b2
! calculate F_k=F^k/d_k**2
             F(i,k)=Fk/d(l,k)**2
21          continue
! perform the second integration to get integral values of F_k
          do j=1,jri
!           value(j)=f(j,k)     changed to Pavels version on 10.6.05
           value(j)=f(j,k)*rho(j)
          enddo
          call charg2(r,dx,value,1,jri,b1)
          b0=value(1)*r(1)/2.d0
          ff(k)=b0+b1
20      continue
!
! From Slater integrals F_k calculate Racah parameters, coef. 2 because
! Hartree -> Ry
! calculation of average Racah parameters
      if(L.eq.2)then
        ABRacah=2.d0*(FF(2)-5.d0*FF(4))
        write(6,*) &
      '     F_0         F_2         F_4          B   '
        write(6,501)ff(0),ff(2),ff(4),abracah
        OP=-ABRacah
      ELSE IF(L.EQ.3)THEN
! Racah parameter E^3 for f-electrons
        AE3=2.d0*(5.d0*FF(2)+6.d0*FF(4)-91.d0*FF(6))/3.d0
        write(6,502)ff(2),ff(4),ff(6),AE3
        OP=-AE3
      endif
! end of calculation of average Racah parameters
501   format(7f12.6)
502   format(' F_2, F_4, F_6, E3, OP para.',5f12.6)
        return
        end
