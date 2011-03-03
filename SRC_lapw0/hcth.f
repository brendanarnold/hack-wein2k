      SUBROUTINE hcth(dfdra, dfdza, dfdrb, dfdzb, dfdzab,  &
                        rhoa, rhob, za, zb, zab, energy, totalF_xc)

!    SUPPLIED TO THE ROUTINE:
!    
!    rhoa   -- value of rhoalpha at a given grid point 
!    rhob   -- value of rhobeta at a given grid point
!    za     -- zeta_alpha, as defined in the TH1 paper (JCP 108 2545), 
!              that is mod(grad(rhoalpha)), a scalar quantity.
!    zb     -- mod(grad(rhobeta)) 
!    zab    -- zeta_{alpha beta} as defined in the TH1 paper, that is
!              grad(rhoalpha).grad(rhobeta) 
!    energy -- a boolean variable deciding whether to compute the energy 
!              contribution at the point in space (true) or the
!              appropriate derivatives (false) needed for the KS matrix.  


!    RETURNED FROM THE ROUTINE:

!    totalF_xc --   the contribution to the energy at this point in space.
!    dfdra  -- partial functional derivative of F_xc with respect to 
!              rhoalpha
!    dfdra  -- partial functional derivative of F_xc with respect to   
!              rhobeta
!    dfdza  -- partial functional derivative of F_xc with respect to   
!              mod(grad(rhoalpha)), divided by za !!!!!!!!!
!    dfdzb  -- partial functional derivative of F_xc with respect to   
!              mod(grad(rhobeta)), divided by zb !!!!!!!!!
!    dfdzab -- partial functional derivative of F_xc with respect to   
!              grad(rhoalpha)*grad(rhobeta) -- is zero for this functional


      IMPLICIT real*8(a-h,o-z)

!fah      PARAMETER(max_pow_u = 2)

       PARAMETER(max_pow_u = 4)


      PARAMETER (PI=3.1415926535898D0)
!fah max_pow_u is equivalent to "m" in the Becke V paper, that is, the greatest 
!fah power of u appearing in the power expansion. 

      DIMENSION sol((max_pow_u+1)*3), F((max_pow_u+1)*3,5),  &
                FF((max_pow_u+1)*3,5,5), &
                F_xc((max_pow_u+1)*3)


!fah sol -- contains the coefficients of the terms in F_xc
!fah        convention: sol(1) = c_{x alpha, 0}, c_{x beta, 0}
!fah                    sol(2) = c_{c alpha alpha, 0}, c_{c beta beta, 0} 
!fah                    sol(3) = c_{c alpha beta, 0} 
!fah                    sol(4) = c_{x alpha, 1}, c_{x beta, 1}
!fah                    sol(5) = c_{c alpha alpha, 1}, c_{c beta beta, 1} 
!fah                    sol(6) = c_{c alpha beta, 1} 
!fah                           
!fah                           etc.
!fah 
!fah f(5) -- contains the partial first functional derivatives of F_xc with 
!fah respect to 
!fah the five quantities (IN THIS ORDER): ra, rb, za, zb, zab  
!fah 
!fah ff(5,5) contains the second derivatives with
!fah respect to the same five quantities

!fah F_xa -- contains the alpha exchange bit containing the various powers 
!fah         of u_{x alpha} (eq. (18) of Becke V paper) 
!fah F_xb --              beta       
!fah            u_{x beta} 
!fah F_caa -- contains the alpha parallel spin correlation bit with the powers
!fah          of u_{c alpha alpha} 
!fah F_cbb --              beta 
!fah             u_{c beta beta} 
!fah F_cab -- contains the anti-parallel spin correlation bit with the powers 
!fah          of u_{c alpha beta} 

!fah these transformed variables u will be defined and given short-cut names 
!fah below. 

      LOGICAL energy

!     Initialise

      dfdra  = 0.D0
      dfdrb  = 0.D0
      dfdza  = 0.D0
      dfdzb  = 0.D0
      dfdzab = 0.D0
      totalF_xc = 0.D0


      IF ((rhoa + rhob) .LT. 1.D-20) RETURN
!fah numerical cutoff: if the density is too low, its contribution is 
!fah neglectable. 

!fah here new solution needs to be inserted after the first run of the 
!fah fit program. 



! please refer to these coeff's as THCH1/iterate-e750-g500-v1-m4-n4
      sol( 1) =     0.109320D+01
      sol( 2) =     0.222601D+00
      sol( 3) =     0.729974D+00
      sol( 4) =    -0.744056D+00
      sol( 5) =    -0.338622D-01
      sol( 6) =     0.335287D+01
      sol( 7) =     0.559920D+01
      sol( 8) =    -0.125170D-01
      sol( 9) =    -0.115430D+02
      sol(10) =    -0.678549D+01
      sol(11) =    -0.802496D+00
      sol(12) =     0.808564D+01
      sol(13) =     0.449357D+01
      sol(14) =     0.155396D+01
      sol(15) =    -0.447857D+01

!fah here new solution needs to be inserted after the first run of the 
!fah fit program. 

      IF (energy) THEN
        CALL deriv(F,FF, &
                   F_xc, &
                   .FALSE., &
                   rhoa,rhob,za,zb, &
                       max_pow_u, &
                   sol)
!fah this logical indicates whether or not the
!fah derivatives are required 
 
        DO n = 0, max_pow_u 
          totalF_xc = totalF_xc +  &
                      F_xc (n*3 + 1) + &
                      F_xc (n*3 + 2) + &
                      F_xc (n*3 + 3)
        ENDDO

      ELSE

        CALL deriv(F,FF, &
                   F_xc, &
                   .TRUE., &
                   rhoa,rhob,za,zb, &
                       max_pow_u, &
                   sol)
!fah this logical indicates whether or not the
!fah derivatives are required 


        DO n = 1, (max_pow_u+1)*3 
          dfdra = dfdra + F(n,1) 
          dfdrb = dfdrb + F(n,2) 
          dfdza = dfdza + F(n,3) / za
          dfdzb = dfdzb + F(n,4) / zb  
        ENDDO
!fah big thanks to NCH: cadpac requires df/(za * dza), NOT 
!fah                                    df/dza 
        dfdzab = 0.D0 
        
      ENDIF

      RETURN
      END
!fah------------------------------------------------------------ 

      FUNCTION F_xs (n, c_xs, rhos, zs)
      IMPLICIT real*8(a-h,o-z)      
!fah cheating about the array size here to avoid introduction of 
!fah another variable: 
      DIMENSION c_xs(1)
      PARAMETER (PI=3.1415926535898D0)
!fah this is one term of the exchange bit of F_xc; n is the power of u in that
!fah term of the power expansion.
      F_xs = (-3.D0*c_xs(1)* &
          (3.D0/Pi)**(1.D0/3.D0)* &
          rhos**(4.D0/3.D0)* &
          ((0.004D0*zs**2.D0)/ &
             (0.004D0*zs**2.D0 + &
               rhos**(8.D0/3.D0)))**n)/ &
        (2.D0*2.D0**(2.D0/3.D0))

      IF (rhos.lt.1.D-10) F_xs = 0.D0
!fah to make safe for the case of beta part of Hydrogen.

      END 

      FUNCTION dF_xs_by_drhos (n, c_xs, rhos, zs)
      IMPLICIT real*8(a-h,o-z)      
      DIMENSION c_xs(1)
      PARAMETER (PI=3.1415926535898D0)
!fah  computes the derivative of the term with u^n of the exchange part of
!fah  F_xc with respect to rho of the same spin.
!fah  n     -- the power of u involved in this term
!fah  c_xs  -- the coefficient c_xs(n) of the term of spin s with the
!fah           power n of u; is NOT passed over as an array.
!fah  rhos -- rhosigma, that is, either rhoalpha or rhobeta
!fah  zs    -- mod(grad(rhosigma)), again for alpha or beta
       dF_xs_by_drhos =  -(c_xs(1)*(6.D0/Pi)** &
            (1.D0/3.D0)* &
           rhos**(1.D0/3.D0)* &
           ((0.004D0*zs**2.D0)/ &
              (rhos** &
                 (8.D0/3.D0) &
                 + 0.004D0*zs**2.D0))**n) &
          + (2.D0*c_xs(1)*0.004D0*n* &
           (6.D0/Pi)** &
            (1.D0/3.D0)* &
           rhos**3.D0*zs**2.D0* &
           ((0.004D0*zs**2.D0)/ &
              (rhos** &
                 (8.D0/3.D0) &
                 + 0.004D0*zs**2.D0))** &
            (-1.D0 + n))/ &
         (rhos**(8.D0/3.D0)+ &
            0.004D0*zs**2.D0)**2.D0
      END

      FUNCTION dF_xs_by_dzs (n, c_xs, rhos, zs)
      IMPLICIT real*8(a-h,o-z)      
      DIMENSION c_xs(1)    
      PARAMETER (PI=3.1415926535898D0)
!fah  idem, but derivation with respect to zs
!fah  see above (function dF_xs_by_drhos) for definition of the
!fah  other variables
      dF_xs_by_dzs = (-3*c_xs(1)*n* &
          (3.D0/Pi)** &
           (1.D0/3.D0)* &
          rhos**(4.D0/3.D0)* &
          ((0.004D0*zs**2.D0)/ &
             (rhos** &
                (8.D0/3.D0)+ &
               0.004D0*zs**2.D0))** &
           (-1.D0 + n)* &
          ((-2.D0*0.004D0**2.D0*zs**3.D0)/ &
             (rhos** &
                 (8.D0/3.D0) &
                 + 0.004D0*zs**2.D0)**2.D0 &
             + (2.D0*0.004D0*zs)/ &
             (rhos** &
                (8.D0/3.D0)+ &
               0.004D0*zs**2.D0)))/ &
        (2.D0*2.D0**(2.D0/3.D0))
      END

      FUNCTION d2F_xs_by_drhos_drhos (n, c_xs, rhos, zs)
      IMPLICIT real*8(a-h,o-z)      
      DIMENSION c_xs(1)
      PARAMETER (PI=3.1415926535898D0)
!fah  this does the second derivative as specified by the name of the
!fah  function.
      d2F_xs_by_drhos_drhos = -((c_xs(1)*(2.D0/Pi)** &
             (1.D0/3.D0)* &
            ((0.004D0*zs**2.D0)/ &
               (rhos** &
                  (8.D0/3.D0) &
                  + 0.004D0*zs**2.D0))**n &
             *((1.D0 - 10.D0*n + 16.D0*n**2.D0)* &
               rhos** &
                (16.D0/3.D0) + &
              2*0.004D0*(1.D0 - 13.D0*n)* &
               rhos** &
                (8.D0/3.D0)* &
               zs**2 + &
              0.004D0**2*zs**4))/ &
          (3**(2.D0/3.D0)* &
            rhos** &
             (2.D0/3.D0)* &
            (rhos** &
                (8.D0/3.D0) + &
               0.004D0*zs**2)**2))
      END

      FUNCTION d2F_xs_by_drhos_dzs (n, c_xs, rhos, zs)
      IMPLICIT real*8(a-h,o-z)      
      DIMENSION c_xs(1)
      PARAMETER (PI=3.1415926535898D0)
!fah  this does the second derivative as specified by the name of the
!fah  function.
      d2F_xs_by_drhos_dzs = (2.D0*c_xs(1)*n* &
          (6.D0/Pi)** &
           (1.D0/3.D0)* &
          rhos**3.D0* &
          ((-1.D0 + 2.D0*n)* &
             rhos**(8.D0/3.D0) &
              - 3.D0*0.004D0*zs**2.D0)* &
          ((0.004D0*zs**2.D0)/ &
             (rhos** &
                (8.D0/3.D0) + &
               0.004D0*zs**2.D0))**n)/ &
        (zs*(rhos** &
              (8.D0/3.D0) + &
             0.004D0*zs**2.D0)**2.D0)
      END

      FUNCTION d2F_xs_by_dzs_dzs (n, c_xs, rhos, zs)
      IMPLICIT real*8(a-h,o-z)      
      DIMENSION c_xs(1)
      PARAMETER (PI=3.1415926535898D0)
!fah  this does the second derivative as specified by the name of the
!fah  function.
      d2F_xs_by_dzs_dzs = (-3.D0*c_xs(1)*n* &
          (3.D0/Pi)** &
           (1.D0/3.D0)* &
          rhos**4.D0* &
          ((-1.D0 + 2.D0*n)* &
             rhos**(8.D0/3.D0) &
              - 3.D0*0.004D0*zs**2.D0)* &
          ((0.004D0*zs**2.D0)/ &
             (rhos** &
                (8.D0/3.D0) + &
               0.004D0*zs**2.D0))**n)/ &
        (2.D0**(2.D0/3.D0)* &
          zs**2.D0* &
          (rhos** &
              (8.D0/3.D0) + &
             0.004D0*zs**2.D0)**2.D0)
      END

!fah------------------------------------------------------------ 

      SUBROUTINE deriv(F,FF, &
                   F_xc, &
                   derivestuff, &
                   rhoa,rhob,za,zb, &
                       max_pow_u, &
                   sol)
!fah this logical indicates whether or not the 
!fah derivatives are required (important for cadpac use only) &

      IMPLICIT NONE
!fah this makes sure I have no typo's. 

      INTEGER max_pow_u
      REAL*8 F_xc((max_pow_u + 1)*3)
      REAL*8 F((max_pow_u+1)*3,5), FF((max_pow_u+1)*3,5,5)
      REAL*8 rhoa, rhob, za, zb

      COMMON/special/h_atom
      LOGICAL h_atom
      LOGICAL derivestuff 
      REAL*8 sol((max_pow_u+1)*3)
    
!fah these are the first derivatives of the different transformed variables 
!fah u with respect to rhoa, rhob, za and zb. These different derivatives 
!fah with respect to these 4 quantities named above are stored in these 
!fah arrays.
      
      REAL*8 dF_xa(4)
      REAL*8 dF_xb(4)
      REAL*8 dF_caa(4)
      REAL*8 dF_cbb(4)
      REAL*8 dF_cab(4)
!fah these are the first derivatives of the terms of F_xc with respect to 
!fah the 4 quantities. the index
!fah runs over the particular partial derivatives of each term.  
!fah More explicitly: these are the partial functional derivatives of 
!fah F_??? with respect to rhoa, rhob, za and zb. 

      REAL*8 d2F_xa(4,4)
      REAL*8 d2F_xb(4,4)
      REAL*8 d2F_caa(4,4)
      REAL*8 d2F_cbb(4,4)
      REAL*8 d2F_cab(4,4)

      REAL*8 Pi 
      PARAMETER (Pi = 3.1415926535898D0)
      REAL*8 rho
      REAL*8 s_a2, s_b2, s_avg2, u_caa, u_cbb, u_cab
      REAL*8 du_caa_by_drhoa, du_caa_by_dza, du_cbb_by_drhob 
      REAL*8 du_cbb_by_dzb, du_cab_by_drhoa, du_cab_by_drhob 
      REAL*8 du_cab_by_dza, du_cab_by_dzb, du_caa_by_drhoa_dza 
      REAL*8 du_caa_by_dza_dza, du_cbb_by_dzb_dzb
      REAL*8 du_cbb_by_drhob_dzb, du_cab_by_drhoa_dza 
      REAL*8 du_cab_by_drhoa_dzb, du_cab_by_drhob_dza 
      REAL*8 du_cab_by_drhob_dzb, du_cab_by_dza_dza, du_cab_by_dza_dzb 
      REAL*8 du_cab_by_dzb_dzb 
      REAL*8 rsa, rsa12, rsa32, rsa21, rsb, rsb12, rsb32, rsb21 
      REAL*8 rsab, rsab12, rsab32, rsab21 
      REAL*8 drsa_by_drhoa, drsb_by_drhob, drsab_by_drhoa
      REAL*8 drsab_by_drhob 
      REAL*8 zeta, dzeta_by_drhoa, dzeta_by_drhob 
      REAL*8 fzeta, dfzeta_by_dzeta, e_crsa1, e_crsb1, e_crsab1 
      REAL*8 e_crsab0, a_crsab
      REAL*8 e_crsabzeta, de_crsa1_by_drsa, de_crsb1_by_drsb 
      REAL*8 da_crsab_by_drsab, de_crsab0_by_drsab 
      REAL*8 de_crsab1_by_drsab, de_crsabzeta_by_drsab 
      REAL*8 de_crsabzeta_by_dzeta, e_caa, e_cbb, e_cab,  &
       de_caa_by_drhoa, de_cbb_by_drhob, de_cab_by_drhoa,  &
       de_cab_by_drhob, &
       c_naa, c_nbb, c_nab
      REAL*8 F_xs ! this is a function which is called. 
      REAL*8 dF_xs_by_drhos, dF_xs_by_dzs,  &
       d2F_xs_by_drhos_dzs, d2F_xs_by_dzs_dzs


      INTEGER i, j, k, n


      DO j = 1, 4
        DO n = 1, (max_pow_u+1)*3
          F(n,j) = 0.D0
!fah  later on, n has a different meaning: n as power of u, not 
!fah  as number of the coefficient. 
        ENDDO
        dF_xa(j) = 0.D0
        dF_xb(j) = 0.D0
        dF_caa(j) = 0.D0
        dF_cbb(j) = 0.D0
        dF_cab(j) = 0.D0
        DO k = 1, 4
          DO n = 1, (max_pow_u+1)*3
            FF(n,j,k) = 0.D0
          ENDDO
          d2F_xa(j,k) = 0.D0
          d2F_xb(j,k) = 0.D0
          d2F_caa(j,k) = 0.D0
          d2F_cbb(j,k) = 0.D0
          d2F_cab(j,k) = 0.D0
        ENDDO
      ENDDO
      DO j = 1, (max_pow_u+1)*3
        F_xc(j) = 0.D0
      ENDDO 

!fah --------------------------------------------------------------

!fah call the expensive correlation parts here just once, and store their
!fah values in a temporary variable. Then compute the actual F_c derivatives
!fah with the various powers of u.  

      rho = rhoa + rhob

      s_a2 = za**2.D0 / rhoa**(8.D0/3.D0)
      s_b2 = zb**2.D0 / rhob**(8.D0/3.D0)
      s_avg2 = 0.5D0*(s_a2 + s_b2)

      u_caa = 0.2D0*s_a2/(1.D0+0.2D0*s_a2) 
      u_cbb = 0.2D0*s_b2/(1.D0+0.2D0*s_b2) 
      u_cab = 0.006D0*s_avg2/(1+0.006D0*s_avg2)

      rsa = ((3/Pi)**(1.D0/3.D0)* &
          (1/rhoa)**(1.D0/3.D0))/ &
        2**(2.D0/3.D0)
      rsa12 = rsa**(1.D0/2.D0)
      rsa32 = rsa**(3.D0/2.D0)
      rsa21 = rsa**2.D0

      rsb = ((3/Pi)**(1.D0/3.D0)* &
          (1/rhob)**(1.D0/3.D0))/ &
        2**(2.D0/3.D0)
      rsb12 = rsb**(1.D0/2.D0)
      rsb32 = rsb**(3.D0/2.D0)
      rsb21 = rsb**2.D0

      rsab = ((3/Pi)**(1.D0/3.D0)* &
          (1/rho)**(1.D0/3.D0))/ &
        2**(2.D0/3.D0)
      rsab12 = rsab**(1.D0/2.D0)
      rsab32 = rsab**(3.D0/2.D0)
      rsab21 = rsab**2.D0

      zeta = (rhoa-rhob)/rho

      fzeta = (-2 + (1 - zeta)** &
           (4.D0/3.D0) + &
          (1 + zeta)**(4.D0/3.D0))/ &
        (-2 + 2*2**(1.D0/3.D0))

      e_crsa1 = -0.03108999999999999* &
        dlog(1 + 32.16468317787069/ &
           (14.1189*rsa12 + &
             6.1977*rsa + 3.3662*rsa32 + &
             0.6251699999999999*rsa21))* &
        (1 + 0.20548*rsa)

      e_crsb1 = -0.03108999999999999* &
        dlog(1 + 32.16468317787069/ &
           (14.1189*rsb12 + &
             6.1977*rsb + 3.3662*rsb32 + &
             0.6251699999999999*rsb21))* &
        (1 + 0.20548*rsb)

      e_crsab1 = -0.03108999999999999* &
        dlog(1 + 32.16468317787069/ &
           (14.1189*rsab12 + &
             6.1977*rsab + 3.3662*rsab32 + &
             0.6251699999999999*rsab21))* &
        (1 + 0.20548*rsab)

      e_crsab0 = -0.062182*dlog(1 + &
          16.0818243221511/ &
           (7.595699999999999*rsab12 + &
             3.5876*rsab + &
             1.6382*rsab32 + &
             0.49294*rsab21))* &
        (1 + 0.2137*rsab)

      a_crsab = 0.03377399999999999* &
        dlog(1 + 29.60857464321667/ &
           (10.35699999999999*rsab12 + &
             3.623099999999999*rsab + &
             0.88026*rsab32 + &
             0.49671*rsab21))* &
        (1 + 0.11125*rsab)

      e_crsabzeta = e_crsab0+a_crsab*fzeta*(1-zeta**4)/1.709921 + &
        (e_crsab1-e_crsab0)*fzeta*zeta**4
  
!fah       print*, 'deriv2', e_crsab0, a_crsab, fzeta, zeta, e_crsab1

      e_caa = rhoa*e_crsa1
      e_cbb = rhob*e_crsb1
      e_cab = rho*e_crsabzeta - rhoa*e_crsa1 - rhob*e_crsb1

!fah       print*, 'deriv1:',rho, e_crsabzeta, e_crsa1, e_crsb1, rhoa, rhob


!fah derive all this stuff only if we need it, not if only F_c itself 
!fah is required. 
      IF (derivestuff.EQV..TRUE.) THEN  

      du_caa_by_drhoa = (-8*0.2D0*za**2*rhoa**(5.D0/3.D0))/ &
        (3.*(0.2D0*za**2 +  &
             rhoa**(8.D0/3.D0))**2)
      du_caa_by_dza = (2*0.2D0*za*rhoa**(8.D0/3.D0))/ &
        (0.2D0*za**2 + rhoa**(8.D0/3.D0))**2

      du_cbb_by_drhob = (-8*0.2D0*zb**2*rhob**(5.D0/3.D0))/ &
        (3.*(0.2D0*zb**2 +  &
             rhob**(8.D0/3.D0))**2)
      du_cbb_by_dzb = (2*0.2D0*zb*rhob**(8.D0/3.D0))/ &
        (0.2D0*zb**2 + rhob**(8.D0/3.D0))**2

      du_cab_by_drhoa = (-16*0.006D0*za**2*rhoa**(5.D0/3.D0)* &
          rhob**(16.D0/3.D0))/ &
        (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) +  &
             0.006D0*za**2*rhob**(8.D0/3.D0) +  &
             2*rhoa**(8.D0/3.D0)* &
              rhob**(8.D0/3.D0))**2) 
      du_cab_by_drhob = (-16*0.006D0*zb**2*rhob**(5.D0/3.D0)* &
          rhoa**(16.D0/3.D0))/ &
        (3.*(0.006D0*za**2*rhob**(8.D0/3.D0) +  &
             0.006D0*zb**2*rhoa**(8.D0/3.D0) +  &
             2*rhob**(8.D0/3.D0)* &
              rhoa**(8.D0/3.D0))**2) 
      du_cab_by_dza = (4*0.006D0*za*rhoa**(8.D0/3.D0)* &
          rhob**(16.D0/3.D0))/ &
        (0.006D0*zb**2*rhoa**(8.D0/3.D0) +  &
           0.006D0*za**2*rhob**(8.D0/3.D0) +  &
           2*rhoa**(8.D0/3.D0)* &
            rhob**(8.D0/3.D0))**2
      du_cab_by_dzb = (4*0.006D0*zb*rhoa**(16.D0/3.D0)* &
          rhob**(8.D0/3.D0))/ &
        (0.006D0*zb**2*rhoa**(8.D0/3.D0) +  &
           0.006D0*za**2*rhob**(8.D0/3.D0) +  &
           2*rhoa**(8.D0/3.D0)* &
            rhob**(8.D0/3.D0))**2

!fah Second derivatives are not required by cadpac. 
!fah   du_caa_by_drhoa_dza = (16*0.2D0*za*rhoa**(5.D0/3.D0)*
!fah -    (0.2D0*za**2 - rhoa**(8.D0/3.D0)))/
!fah -  (3.*(0.2D0*za**2 + 
!fah -       rhoa**(8.D0/3.D0))**3)
!fah   du_cbb_by_drhob_dzb = (16*0.2D0*zb*rhob**(5.D0/3.D0)*
!fah -    (0.2D0*zb**2 - rhob**(8.D0/3.D0)))/
!fah -  (3.*(0.2D0*zb**2 + 
!fah -       rhob**(8.D0/3.D0))**3)
!fah
!fah   du_caa_by_dza_dza = (2*0.2D0*rhoa**(8.D0/3.D0)*
!fah -    (-3*0.2D0*za**2 + rhoa**(8.D0/3.D0))
!fah -    )/
!fah -  (0.2D0*za**2 + rhoa**(8.D0/3.D0))**3
!fah   du_cbb_by_dzb_dzb = (2*0.2D0*rhob**(8.D0/3.D0)*
!fah -    (-3*0.2D0*zb**2 + rhob**(8.D0/3.D0))
!fah -    )/
!fah -  (0.2D0*zb**2 + rhob**(8.D0/3.D0))**3
!fah
!fah   du_cab_by_drhoa_dza = (-32*0.006D0*rhoa**(5.D0/3.D0)*
!fah -    (0.006D0*za*zb**2*
!fah -       rhoa**(8.D0/3.D0)*
!fah -       rhob**(16.D0/3.D0) - 
!fah -      0.006D0*za**3*rhob**8 + 
!fah -      2*za*rhoa**(8.D0/3.D0)*rhob**8))/
!fah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -       2*rhoa**(8.D0/3.D0)*
!fah -        rhob**(8.D0/3.D0))**3) 
!fah   du_cab_by_drhoa_dzb = (64*0.006D0**2*za**2*zb*
!fah -    rhoa**(13.D0/3.D0)*
!fah -    rhob**(16.D0/3.D0))/
!fah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -       2*rhoa**(8.D0/3.D0)*
!fah -        rhob**(8.D0/3.D0))**3) 
!fah   du_cab_by_drhob_dza = (64*0.006D0**2*za*zb**2*
!fah -    rhoa**(16.D0/3.D0)*
!fah -    rhob**(13.D0/3.D0))/
!fah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -       2*rhoa**(8.D0/3.D0)*
!fah -        rhob**(8.D0/3.D0))**3) 
!fah   du_cab_by_drhob_dzb = (-32*0.006D0*rhob**(5.D0/3.D0)*
!fah -    (-(0.006D0*zb**3*rhoa**8) + 
!fah -      0.006D0*za**2*zb*
!fah -       rhoa**(16.D0/3.D0)*
!fah -       rhob**(8.D0/3.D0) + 
!fah -      2*zb*rhoa**8*rhob**(8.D0/3.D0)))/
!fah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -       2*rhoa**(8.D0/3.D0)*
!fah -        rhob**(8.D0/3.D0))**3) 
!fah   du_cab_by_dza_dza = (4*0.006D0*(0.006D0*zb**2*
!fah -       rhoa**(16.D0/3.D0)*
!fah -       rhob**(16.D0/3.D0) - 
!fah -      3*0.006D0*za**2*rhoa**(8.D0/3.D0)*
!fah -       rhob**8 + 2*rhoa**(16.D0/3.D0)*rhob**8
!fah -      ))/
!fah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -     2*rhoa**(8.D0/3.D0)*
!fah -      rhob**(8.D0/3.D0))**3 
!fah   du_cab_by_dza_dzb = (-16*0.006D0**2*za*zb*
!fah -    rhoa**(16.D0/3.D0)*
!fah -    rhob**(16.D0/3.D0))/
!fah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -     2*rhoa**(8.D0/3.D0)*
!fah -      rhob**(8.D0/3.D0))**3 
!fah   du_cab_by_dzb_dzb = (4*0.006D0*rhoa**(16.D0/3.D0)*
!fah -    rhob**(8.D0/3.D0)*
!fah -    (-3*0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -      0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -      2*rhoa**(8.D0/3.D0)*
!fah -       rhob**(8.D0/3.D0)))/
!fah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
!fah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
!fah -     2*rhoa**(8.D0/3.D0)*
!fah -      rhob**(8.D0/3.D0))**3 

      drsa_by_drhoa = -((1/rhoa)**(4.D0/3.D0)/ &
          (6**(2.D0/3.D0)*Pi**(1.D0/3.D0)))
      drsb_by_drhob = -((1/rhob)**(4.D0/3.D0)/ &
          (6**(2.D0/3.D0)*Pi**(1.D0/3.D0)))
      drsab_by_drhoa = -((1/rho)**(4.D0/3.D0)/ &
          (6**(2.D0/3.D0)*Pi**(1.D0/3.D0))) 
      drsab_by_drhob = drsab_by_drhoa 

      dzeta_by_drhoa = 2*rhob/rho**2
      dzeta_by_drhob = -2*rhoa/rho**2

      dfzeta_by_dzeta = ((-4*(1 - zeta)**(1.D0/3.D0))/ &
           3. + (4*(1 + zeta)** &
              (1.D0/3.D0))/3.)/ &
        (-2 + 2*2**(1.D0/3.D0))

      de_crsa1_by_drsa = (1.*(1 + 0.20548*rsa)* &
           (6.1977 + 7.05945/rsa12 + 1.25034*rsa +  &
             5.0493*rsa12))/ &
         ((6.1977*rsa + 14.1189*rsa12 +  &
              0.6251699999999999*rsa21 + 3.3662*rsa32)** &
            2*(1 + 32.16468317787069/ &
              (6.1977*rsa + 14.1189*rsa12 +  &
                0.6251699999999999*rsa21 + 3.3662*rsa32) &
             )) - 0.006388373199999999* &
         dlog(1 + 32.16468317787069/ &
            (6.1977*rsa + 14.1189*rsa12 +  &
              0.6251699999999999*rsa21 + 3.3662*rsa32)) 

      de_crsb1_by_drsb = (1.*(1 + 0.20548*rsb)* &
           (6.1977 + 7.05945/rsb12 + 1.25034*rsb +  &
             5.0493*rsb12))/ &
         ((6.1977*rsb + 14.1189*rsb12 +  &
              0.6251699999999999*rsb21 + 3.3662*rsb32)** &
            2*(1 + 32.16468317787069/ &
              (6.1977*rsb + 14.1189*rsb12 +  &
                0.6251699999999999*rsb21 + 3.3662*rsb32) &
             )) - 0.006388373199999999* &
         dlog(1 + 32.16468317787069/ &
            (6.1977*rsb + 14.1189*rsb12 +  &
              0.6251699999999999*rsb21 + 3.3662*rsb32)) 

      da_crsab_by_drsab = (-1.*(1 + 0.11125*rsab)* &
           (3.623099999999999 +  &
             5.178499999999999/rsab12 +  &
             0.99342*rsab + 1.32039*rsab12))/ &
         ((1 + 29.60857464321667/ &
              (3.623099999999999*rsab +  &
                10.35699999999999*rsab12 +  &
                0.49671*rsab21 + 0.88026*rsab32))* &
           (3.623099999999999*rsab +  &
              10.35699999999999*rsab12 +  &
              0.49671*rsab21 + 0.88026*rsab32)**2) +  &
        0.003757357499999999* &
         dlog(1 + 29.60857464321667/ &
            (3.623099999999999*rsab +  &
              10.35699999999999*rsab12 +  &
              0.49671*rsab21 + 0.88026*rsab32)) 

      de_crsab0_by_drsab = (1.*(1 + 0.2137*rsab)* &
           (3.5876 + 3.797849999999999/rsab12 +  &
             0.98588*rsab + 2.4573*rsab12))/ &
         ((3.5876*rsab + 7.595699999999999*rsab12 +  &
              0.49294*rsab21 + 1.6382*rsab32)**2* &
           (1 + 16.0818243221511/ &
              (3.5876*rsab + 7.595699999999999*rsab12 +  &
                0.49294*rsab21 + 1.6382*rsab32))) -  &
        0.01328829339999999* &
         dlog(1 + 16.0818243221511/ &
            (3.5876*rsab + 7.595699999999999*rsab12 +  &
              0.49294*rsab21 + 1.6382*rsab32))

      de_crsab1_by_drsab = (1.*(1 + 0.20548*rsab)* &
           (6.1977 + 7.05945/Sqrt(rsab) +  &
             1.25034*rsab + 5.0493*rsab12))/ &
         ((6.1977*rsab + 14.1189*rsab12 +  &
              0.6251699999999999*rsab21 + 3.3662*rsab32) &
             **2*(1 + 32.16468317787069/ &
              (6.1977*rsab + 14.1189*rsab12 +  &
                0.6251699999999999*rsab21 +  &
                3.3662*rsab32))) -  &
        0.006388373199999999* &
         dlog(1 + 32.16468317787069/ &
            (6.1977*rsab + 14.1189*rsab12 +  &
              0.6251699999999999*rsab21 + 3.3662*rsab32) &
           )

      de_crsabzeta_by_drsab = 1.124999956683108*(1 - zeta**4)* &
         (-2 + (1 - zeta)**(4.D0/3.D0) +  &
           (1 + zeta)**(4.D0/3.D0))* &
         da_crsab_by_drsab +  &
        de_crsab0_by_drsab +  &
        fzeta*zeta**4* &
         (- de_crsab0_by_drsab +  &
            de_crsab1_by_drsab  )

      de_crsabzeta_by_dzeta = 1.499999942244144*(-1 + zeta**4)* &
         ((1 - zeta)**(1.D0/3.D0) -  &
           (1 + zeta)**(1.D0/3.D0))*a_crsab -  &
        4.499999826732434*zeta**3* &
         (-2 + (1 - zeta)**(4.D0/3.D0) +  &
           (1 + zeta)**(4.D0/3.D0))*a_crsab +  &
        (2*zeta**4*((1 - zeta)**(1.D0/3.D0) -  &
             (1 + zeta)**(1.D0/3.D0))* &
           (e_crsab0 - e_crsab1))/ &
         (3.*(-1 + 2**(1.D0/3.D0))) +  &
        4*fzeta*(-e_crsab0 +  &
           e_crsab1)*zeta**3

!fah this is with application of the chain rule; I keep it that general
!fah because this way, I only have to define one "G". 
      de_caa_by_drhoa = e_crsa1 + rhoa*de_crsa1_by_drsa*drsa_by_drhoa 
      de_cbb_by_drhob = e_crsb1 + rhob*de_crsb1_by_drsb*drsb_by_drhob 

      de_cab_by_drhoa = -e_crsa1 +  &
        e_crsabzeta -  &
        rhoa*de_crsa1_by_drsa*  &
        drsa_by_drhoa +  &
        rho*(de_crsabzeta_by_drsab* &
        drsab_by_drhoa +  &
        de_crsabzeta_by_dzeta* &
        dzeta_by_drhoa)

      de_cab_by_drhob = -e_crsb1 +  &
        e_crsabzeta -  &
        rhob*de_crsb1_by_drsb*  &
        drsb_by_drhob +  &
        rho*(de_crsabzeta_by_drsab* &
        drsab_by_drhob +  &
        de_crsabzeta_by_dzeta* &
        dzeta_by_drhob)

      ENDIF ! if (derivestuff) 


!fah Here starts the big outer loop over the powers u 

      DO n = 0, max_pow_u 
        c_naa = sol((n*3) + 2)
        c_nbb = c_naa
        c_nab = sol((n*3) + 3) 

!fah construction of the F_xc itself
!fah -------------------------------
        IF (rhoa.GT.1.D-10) THEN
          F_xc(n*3+1) = F_xs (n, sol((n*3) + 1), rhoa, za)
          F_xc(n*3+2) = e_caa*u_caa**n*c_naa
        ENDIF

        IF (rhob.GT.1.D-10) THEN
          F_xc(n*3+1) = F_xc(n*3+1)+F_xs(n, sol((n*3) + 1), rhob, zb)
          F_xc(n*3+2) = F_xc(n*3+2)+e_cbb*u_cbb**n*c_nbb
        ENDIF

        IF (rhoa.GT.1.D-10 .AND. rhob.GT.1.D-10) THEN
          F_xc(n*3+3) = e_cab*u_cab**n*c_nab
        ENDIF

!fah       print*, 'in deriv:', e_cab, u_cab, c_nab

!fah    First Derivatives
!fah ---------------------

        dF_xa(1) = dF_xs_by_drhos (n, sol((n*3) + 1), rhoa, za) 
        dF_xa(2) = 0 
        dF_xa(3) = dF_xs_by_dzs (n, sol((n*3) + 1), rhoa, za)
        dF_xa(4) = 0 

        dF_xb(1) = 0 
        dF_xb(2) = dF_xs_by_drhos (n, sol((n*3) + 1), rhob, zb)
        dF_xb(3) = 0 
        dF_xb(4) = dF_xs_by_dzs (n, sol((n*3) + 1), rhob, zb)

        dF_caa(1) = c_naa*u_caa**n*de_caa_by_drhoa +   &
         c_naa*n*e_caa*u_caa**(-1 + n)* &
         du_caa_by_drhoa
        dF_caa(2) = 0.D0
        dF_caa(3) = c_naa*n*e_caa*u_caa**(-1 + n)* &
        du_caa_by_dza
        dF_caa(4) = 0.D0

        dF_cbb(1) = 0.D0
        dF_cbb(2) = c_nbb*u_cbb**n*de_cbb_by_drhob +  &
         c_nbb*n*e_cbb*u_cbb**(-1 + n)* &
         du_cbb_by_drhob
        dF_cbb(3) = 0.D0
        dF_cbb(4) = c_nbb*n*e_cbb*u_cbb**(-1 + n)* &
        du_cbb_by_dzb


        dF_cab(1) = c_nab*u_cab**n* &
         de_cab_by_drhoa +  &
         c_nab*n*e_cab* &
         u_cab**(-1 + n)* &
         du_cab_by_drhoa 
        dF_cab(2) = c_nab*u_cab**n* &
         de_cab_by_drhob + &
         c_nab*n*e_cab* &
         u_cab**(-1 + n)* &
         du_cab_by_drhob
        dF_cab(3) = c_nab*n*e_cab* &
        u_cab**(-1 + n)* &
        du_cab_by_dza
        dF_cab(4) = c_nab*n*e_cab* &
        u_cab**(-1 + n)* &
        du_cab_by_dzb

!fah Second Derivatives
!fah ------------------

!fah         d2F_xa(1,1) = d2F_xs_by_drhos_drhos (n, sol((n*3) + 1), 
!fah      &                                         rhoa, za)
!fah see comment below, for the (2,2) term. 
!fah    d2F_xa(1,2) = 0 
!fah    d2F_xa(1,3) = d2F_xs_by_drhos_dzs (n, sol((n*3) + 1), rhoa, za)
!fah    d2F_xa(1,4) = 0 
!fah    d2F_xa(2,2) = 0  
!fah    d2F_xa(2,3) = 0  
!fah    d2F_xa(2,4) = 0 
!fah    d2F_xa(3,3) = d2F_xs_by_dzs_dzs (n, sol((n*3) + 1), rhoa, za)
!fah    d2F_xa(3,4) = 0 
!fah    d2F_xa(4,4) = 0 

!fah for alpha spin, elements are non-zero when both indices are odd; 
!fah for beta spin, elements are non-zero when both indices are even. 
!fah the matrix is symmetric, and the upper triangle contains the 
!fah 10 elements given above and below. 

!fah    d2F_xb(1,1) = 0
!fah    d2F_xb(1,2) = 0  
!fah    d2F_xb(1,3) = 0           
!fah    d2F_xb(1,4) = 0 
!fah        d2F_xb(2,2) = d2F_xs_by_drhos_drhos (n, sol((n*3) + 1), 
!fah     &                                         rhob, zb)
!fah this term is NOT zero, but needs not be evaluated since we don't 
!fah need it for the construction of v (cf. routine "va" in the fit 
!fah program) 
!fah    d2F_xb(2,3) = 0  
!fah    d2F_xb(2,4) = d2F_xs_by_drhos_dzs (n, sol((n*3) + 1), rhob, zb)
!fah    d2F_xb(3,3) = 0
!fah    d2F_xb(3,4) = 0 
!fah    d2F_xb(4,4) = d2F_xs_by_dzs_dzs (n, sol((n*3) + 1), rhob, zb)


!fah    d2F_caa(1,1) = !=0, but not needed 
!fah    d2F_caa(1,2) = 0.D0 (not needed)  
!fah    d2F_caa(1,3) = c_naa*n*u_caa**(-1 + n)*
!fah -   de_caa_by_drhoa*
!fah -   du_caa_by_dza + 
!fah -  c_naa*(-1 + n)*n*e_caa*
!fah -   u_caa**(-2 + n)*
!fah -   du_caa_by_dza*
!fah -   du_caa_by_drhoa + 
!fah -  c_naa*n*e_caa*u_caa**(-1 + n)*
!fah -   du_caa_by_drhoa_dza
!fah    d2F_caa(1,4) = 0.D0
!fah    d2F_caa(2,2) = 0.D0 (not needed)
!fah    d2F_caa(2,3) = 0.D0 
!fah    d2F_caa(2,4) = 0.D0 
!fah    d2F_caa(3,3) = c_naa*n*e_caa*u_caa**(-2 + n)*
!fah -  ((-1 + n)*du_caa_by_dza**
!fah -      2 + u_caa*
!fah -     du_caa_by_dza_dza)
!fah    d2F_caa(3,4) = 0.D0 
!fah    d2F_caa(4,4) = 0.D0 


!fah    d2F_cbb(1,1) = 0.D0 (not needed)
!fah    d2F_cbb(1,2) = 0.D0 (not needed)
!fah    d2F_cbb(1,3) = 0.D0
!fah    d2F_cbb(1,4) = 0.D0
!fah    d2F_cbb(2,2) = !=0, but not needed 
!fah    d2F_cbb(2,3) = 0.D0
!fah    d2F_cbb(2,4) = c_nbb*n*u_cbb**(-1 + n)*
!fah -   de_cbb_by_drhob*
!fah -   du_cbb_by_dzb +
!fah -  c_nbb*(-1 + n)*n*e_cbb*
!fah -   u_cbb**(-2 + n)*
!fah -   du_cbb_by_dzb*
!fah -   du_cbb_by_drhob +
!fah -  c_nbb*n*e_cbb*u_cbb**(-1 + n)*
!fah -   du_cbb_by_drhob_dzb
!fah    d2F_cbb(3,3) = 0.D0
!fah    d2F_cbb(3,4) = 0.D0
!fah    d2F_cbb(4,4) =  c_nbb*n*e_cbb*u_cbb**(-2 + n)*
!fah -  ((-1 + n)*du_cbb_by_dzb**
!fah -      2 + u_cbb*
!fah -     du_cbb_by_dzb_dzb)

!fah    d2F_cab(1,1) = not needed
!fah    d2F_cab(1,2) = not needed
!fah    d2F_cab(1,3) = c_nab*n*u_cab**(-2 + n)*
!fah -  ((-1 + n)*e_cab*
!fah -     du_cab_by_dza
!fah -      *du_cab_by_drhoa + 
!fah -    u_cab*(de_cab_by_drhoa*
!fah -        du_cab_by_dza + e_cab*
!fah -        du_cab_by_drhoa_dza
!fah -  ))
!fah    d2F_cab(1,4) = c_nab*n*u_cab**(-2 + n)*
!fah -  (  (-1 + n)*e_cab*
!fah -     du_cab_by_dzb
!fah -      *du_cab_by_drhoa + 
!fah -    u_cab*
!fah -     (de_cab_by_drhoa*
!fah -        du_cab_by_dzb + 
!fah -       e_cab*
!fah -        du_cab_by_drhoa_dzb))
!fah    d2F_cab(2,2) = not needed
!fah    d2F_cab(2,3) = c_nab*n*u_cab**(-2 + n)*
!fah -  ((-1 + n)*e_cab*
!fah -     du_cab_by_dza
!fah -      *du_cab_by_drhob +
!fah -    u_cab*
!fah -     (de_cab_by_drhob*
!fah -        du_cab_by_dza + 
!fah -       e_cab*
!fah -        du_cab_by_drhob_dza))
!fah    d2F_cab(2,4) = c_nab*n*u_cab**(-2 + n)*
!fah -  ((-1 + n)*e_cab*
!fah -     du_cab_by_dzb
!fah -      *du_cab_by_drhob + 
!fah -    u_cab*(de_cab_by_drhob*
!fah -        du_cab_by_dzb + e_cab*
!fah -        du_cab_by_drhob_dzb ))
!fah    d2F_cab(3,3) = c_nab*n*e_cab*
!fah -  u_cab**(-2 + n)*
!fah -  ((-1 + n)*du_cab_by_dza**2 + 
!fah -    u_cab*
!fah -     du_cab_by_dza_dza)
!fah    d2F_cab(3,4) = c_nab*n*e_cab*
!fah -  u_cab**(-2 + n)*
!fah -  ((-1 + n)*du_cab_by_dzb*
!fah -     du_cab_by_dza
!fah -       + u_cab*
!fah -     du_cab_by_dza_dzb)
!fah    d2F_cab(4,4) = c_nab*n*e_cab*
!fah -  u_cab**(-2 + n)*
!fah -  ((-1 + n)*du_cab_by_dzb**2 + 
!fah -    u_cab*
!fah -     du_cab_by_dzb_dzb)
!fah
!fah here, the second derivatives are completed (Schwartz's rule: 
!fah df/(dadb) = df/(dbda) 
!fah    DO i = 1, 4
!fah      DO j = i, 4
!fah        d2F_xa(j,i) = d2F_xa(i,j)  
!fah        d2F_xb(j,i) = d2F_xb(i,j)
!fah        d2F_caa(j,i) = d2F_caa(i,j) 
!fah        d2F_cbb(j,i) = d2F_cbb(i,j) 
!fah        d2F_cab(j,i) = d2F_cab(i,j) 
!fah      ENDDO
!fah    ENDDO

!fah test for zero densities (as in beta part of H atom):
        IF (rhob.LT.1.D-10) THEN 
          DO i = 1, 4
            dF_xb(i) = 0.D0
            dF_cbb(i) = 0.D0
            dF_cab(i) = 0.D0
!fah        DO j = 1, 4
!fah          d2F_xb(i,j) = 0.D0
!fah          d2F_cbb(i,j) = 0.D0
!fah          d2F_cab(i,j) = 0.D0
!fah        ENDDO 
          ENDDO 
        ENDIF

        IF (rhoa.LT.1.D-10) THEN
          DO i = 1, 4
            dF_xa(i) = 0.D0
            dF_caa(i) = 0.D0
            dF_cab(i) = 0.D0
!fah        DO j = 1, 4
!fah          d2F_xa(i,j) = 0.D0
!fah          d2F_caa(i,j) = 0.D0
!fah          d2F_cab(i,j) = 0.D0
!fah        ENDDO
          ENDDO
        ENDIF


!fah Sum up all the partial derivatives with respect to the same function
!fah of terms containing different powers of u with the help of the big outer 
!fah loop: 

!fah have the partial derivative 

        DO i = 1, 4
          F(n*3+1,i) = dF_xa(i) + dF_xb(i) 
          F(n*3+2,i) = dF_caa(i) +dF_cbb(i) 
          F(n*3+3,i) = dF_cab(i) 
!fah      DO j = 1, 4
!fah        FF(n*3+1,i,j) = d2F_xa(i,j) + d2F_xb(i,j) 
!fah        FF(n*3+2,i,j) = d2F_caa(i,j) + d2F_cbb(i,j) 
!fah        FF(n*3+3,i,j) = d2F_cab(i,j)  
!fah      ENDDO 
        ENDDO 

!fah these partial derivatives have not been computed because they are
!fah zero since we don't have a gradrhoagradrhob term in the Becke V functional
        F(n*3+1,5) = 0
        F(n*3+2,5) = 0
        F(n*3+3,5) = 0
!fah    DO i = 1, 5
!fah      FF(n*3+1,i,5) = 0
!fah      FF(n*3+2,i,5) = 0
!fah      FF(n*3+3,i,5) = 0
!fah      FF(n*3+1,5,i) = 0
!fah      FF(n*3+2,5,i) = 0
!fah      FF(n*3+3,5,i) = 0
!fah    ENDDO

      ENDDO

      RETURN

      END

!fah-----------------------------------------------------------

