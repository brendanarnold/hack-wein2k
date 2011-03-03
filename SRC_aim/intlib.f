subroutine avint ( ftab, xtab, ntab, a, b, result )
!
!***********************************************************************
!
!! AVINT estimates the integral of unevenly spaced data.
!
!
!  Discussion:
!
!    The method uses overlapping parabolas and smoothing.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    P E Hennion,
!    Algorithm 77,
!    Interpolation, Differentiation and Integration,
!    Communications of the Association for Computing Machinery,
!    Volume 5, page 96, 1962.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), the function values,
!    FTAB(I) = F(XTAB(I)).
!
!    Input, real*8 XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, integer NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real*8 A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real*8 B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer ntab
!
  real*8 a
  real*8 atemp
  real*8 b
  real*8 btemp
  real*8 ca
  real*8 cb
  real*8 cc
  real*8 ctemp
  real*8 ftab(ntab)
  integer i
  integer ihi
  integer ilo
  integer ind
  real*8 result
  real*8 sum1
  real*8 syl
  real*8 term1
  real*8 term2
  real*8 term3
  real*8 x1
  real*8 x2
  real*8 x3
  real*8 xtab(ntab)
!
  if ( ntab < 3 ) then
    write ( 6, '(a)' ) ' '
    write ( 6, '(a)' ) 'AVINT - Fatal error!'
    write ( 6, '(a,i6)' ) '  NTAB is less than 3.  NTAB = ', ntab
    results=0
    return
!    stop
  end if
 
  do i = 2, ntab
 
    if ( xtab(i) <= xtab(i-1) ) then
      write ( 6, '(a)' ) ' '
      write ( 6, '(a)' ) 'AVINT - Fatal error!'
      write ( 6, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
      write ( 6, '(a,i6)' ) '  Here, I = ', I
      write ( 6, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
      write ( 6, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
      result=0
      return
!      stop
    end if
 
  end do
 
  result = 0.0D+00
 
  if ( a == b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Warning!'
    write ( *, '(a)' ) '  A = B, integral=0.'
    return
  end if
!
!  If A > B, temporarily switch A and B, and store sign.
!
  if ( a > b ) then
    syl = b
    b = a
    a = syl
    ind = -1
  else
    syl = a
    ind = 1
  end if
!
!  Bracket A and B between XTAB(ILO) and XTAB(IHI).
!
  ilo = 1
  ihi = ntab

  do i = 1, ntab
    if ( xtab(i) >= a ) then
      exit
    end if
    ilo = ilo + 1
  end do

  ilo = max ( 2, ilo )
  ilo = min ( ilo, ntab-1 )

  do i = 1, ntab
    if ( b >= xtab(i) ) then
      exit
    end if
    ihi = ihi - 1
  end do
  
  ihi = min ( ihi, ntab-1 )
  ihi = max ( ilo, ihi-1 )
!
!  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
!
  sum1 = 0.0D+00
 
  do i = ilo, ihi
 
    x1 = xtab(i-1)
    x2 = xtab(i)
    x3 = xtab(i+1)
 
    term1 = ftab(i-1) / ((x1-x2)*(x1-x3))
    term2 = ftab(i) / ((x2-x1)*(x2-x3))
    term3 = ftab(i+1) / ((x3-x1)*(x3-x2))
 
    atemp = term1 + term2 + term3
    btemp = -(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3
    ctemp = x2*x3*term1+x1*x3*term2+x1*x2*term3
 
    if ( i <= ilo ) then
      ca = atemp
      cb = btemp
      cc = ctemp
    else
      ca = 0.5D+00 * ( atemp + ca )
      cb = 0.5D+00 * ( btemp + cb )
      cc = 0.5D+00 * ( ctemp + cc )
    end if
 
    sum1 = sum1 &
          + ca * ( x2**3 - syl**3 ) / 3.0D+00 &
          + cb * 0.5D+00 * ( x2**2 - syl**2 ) &
          + cc * ( x2 - syl )
 
    ca = atemp
    cb = btemp
    cc = ctemp
 
    syl = x2
 
  end do
 
  result = sum1 &
        + ca * ( b**3 - syl**3 ) / 3.0D+00 &
        + cb * 0.5D+00 * ( b**2 - syl**2 ) &
        + cc * ( b - syl )
!
!  Restore original values of A and B, reverse sign of integral
!  because of earlier switch.
!
  if ( ind /= 1 ) then
    ind = 1
    syl = b
    b = a
    a = syl
    result = -result
  end if
 
  return
end
subroutine cadre ( func, a, b, abserr, relerr, error, result, ind )
!
!***********************************************************************
!
!! CADRE estimates the integral of F(X) from A to B.
!
!
!  Discussion:
!
!    CADRE is the Cautious Adaptive Romberg Extrapolator.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    Carl DeBoor and J R Rice,
!    CADRE: An algorithm for numerical quadrature,
!    Mathematic Software, pages 417-449,
!    Academic Press, New York, 1971.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC, the name of the function to be integrated.
!    The user must declare the name an external parameter in the calling
!    program, write a function routine of the form FUNCTION FUNC(X) which
!    evaluates the function at X, and pass the name of the function
!    in FUNC.
!
!    Input, real*8 A, the lower limit of integration.
!
!    Input, real*8 B, the upper limit of integration.
!
!    Input, real*8 ABSERR, the absolute error tolerance.
!
!    Input, real*8 RELERR, the relative error tolerance.
!
!    Output, real*8 ERROR, an estimate of the absolute error.
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
!    Output, integer IND, reliability indicator.
!    If IND <= 2, RESULT is very reliable.  Higher values of
!    IND indicate less reliable values of RESULT.
!
  implicit real*8 (a-h,o-z)
!
  integer, parameter :: mxstge = 30
  integer, parameter :: maxtbl = 10
  integer, parameter :: maxts = 2049
!
  real*8 a
  real*8 abserr
  real*8 ait(maxtbl)
  logical aitken
  real*8, parameter :: aitlow = 1.1D+00
  real*8, parameter :: aittol = 0.1D+00
  real*8 astep
  real*8 b
  real*8 beg
  real*8 begin(mxstge)
  real*8 bma
  real*8 curest
  real*8 dif(maxtbl)
  real*8 diff
  real*8 end
  real*8 ergoal
  real*8 erra
  real*8 errer
  real*8 error
  real*8 errr
  real*8 est(mxstge)
  real*8 fbeg
  real*8 fbeg2
  real*8 fend
  real*8 fextm1
  real*8 fextrp
  real*8 finis(mxstge)
  real*8 fn
  real*8 fnsize
  real*8, external :: func
  logical h2conv
  real*8 h2next
  real*8 h2tfex
  real*8, parameter :: h2tol = 0.15D+00
  real*8 hovn
  integer i
  integer ibeg
  integer ibegs(mxstge)
  integer iend
  integer ii
  integer iii
  integer ind
  integer istage
  integer istep
  integer istep2
  integer it
  integer l
  integer lm1
  integer n
  integer n2
  integer nnleft
  real*8 prever
  real*8 r(maxtbl)
  logical reglar
  logical reglsv(mxstge)
  real*8 relerr
  real*8 result
  logical right
  real*8 rn(4)
  real*8 rnderr
  real*8 sing
  real*8 singnx
  real*8 slope
  real*8 stage
  real*8 step
  real*8 stepmn
  real*8 sum1
  real*8 sumabs
  real*8 t(maxtbl,maxtbl)
  real*8 tabs
  real*8 tabtlm
  real*8, parameter :: tljump = 0.01D+00
  real*8 ts(2049)
  real*8 vint
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  begin(1:mxstge) = 0.0D+00
  est(1:mxstge) = 0.0D+00
  finis(1:mxstge) = 0.0D+00
  ibegs(1:mxstge) = 0
  reglsv(1:mxstge) = .false.
 
  vint = 0.0D+00
 
  rn(1:4) = (/ 0.7142005D+00, 0.3466282D+00, 0.8437510D+00, 0.1263305D+00 /)
 
  rnderr = epsilon ( rnderr )
  result = 0.0D+00
  error = 0.0D+00
  ind = 1
  bma = abs ( b - a )
  errr = min ( 0.1D+00, max ( abs ( relerr ), 10.0D+00*rnderr) )
  erra = abs ( abserr )
  stepmn = max ( bma / 2**mxstge, max ( bma, abs ( a ), abs ( b ) ) * rnderr )
  stage = 0.5
  istage = 1
  curest = 0.0D+00
  fnsize = 0.0D+00
  prever = 0.0D+00
  reglar = .false.
  beg = a
  fbeg = func(beg) / 2.0D+00
  ts(1) = fbeg
  ibeg = 1
  end = b
  fend = func(end) / 2.0D+00
  ts(2) = fend
  iend = 2
 
10 continue
 
  right = .false.
 
20 continue

  step = end - beg
  astep = abs ( step )
 
  if ( astep < stepmn ) then
    ind = 5
    result = curest + vint
    return
  end if
 
  t(1,1) = fbeg+fend
  tabs = abs ( fbeg ) + abs ( fend )
  l = 1
  n = 1
  h2conv = .false.
  aitken = .false.
  go to 40
 
30 continue
 
40 continue
 
  lm1 = l
  l = l+1
  n2 = n*2
  fn = n2
  istep = (iend-ibeg)/n

  if ( istep > 1 ) then
    go to 60
  end if

  ii = iend
  iend = iend + n

  if ( iend > maxts ) then
    go to 440
  end if

  hovn = step / fn
 
  iii = iend
  do i = 1, n2, 2
    ts(iii) = ts(ii)
    ts(iii-1) = func(end-i*hovn)
    iii = iii-2
    ii = ii-1
  end do
 
  istep = 2
 
60 continue
 
  istep2 = ibeg+istep/2
 
  sum1 = 0.0D+00
  sumabs = 0.0D+00
  do i = istep2, iend, istep
    sum1 = sum1 + ts(i)
    sumabs = sumabs + abs ( ts(i) )
  end do
 
  t(l,1) = t(l-1,1) / 2.0D+00 + sum1 / fn
  tabs = tabs / 2.0D+00 + sumabs / fn
 
  n = n2
  it = 1
  vint = step * t(l,1)
  tabtlm = tabs * rnderr
  fnsize = max ( fnsize, abs ( t(l,1) ) )
  ergoal = max ( astep * rnderr * fnsize, &
    stage * max ( erra , errr * abs ( curest+vint ) ) )
  fextrp = 1.0D+00
  do i = 1, lm1
    fextrp = fextrp * 4.0D+00
    t(i,l) = t(l,i) - t(l-1,i)
    t(l,i+1) = t(l,i) + t(i,l) / ( fextrp - 1.0 )
  end do
 
  errer = astep * abs ( t(1,l) )
  if ( l > 2 ) go to 90
  if ( abs ( t(1,2) ) <= tabtlm ) go to 290
  go to 40
 
90 continue
 
  do i = 2, lm1

    if ( abs ( t(i-1,l) ) > tabtlm ) then
      diff = t(i-1,lm1) / t(i-1,l)
    else
      diff = 0.0D+00
    end if

    t(i-1,lm1) = diff

  end do
 
  if ( abs ( 4.0 - t(1,lm1) ) <= h2tol ) go to 130
  if ( t(1,lm1) == 0.0 ) go to 120
  if ( abs ( 2.0 - abs ( t(1,lm1) ) ) < tljump ) go to 280
  if (l==3) go to 30
  h2conv = .false.
  if ( abs ( ( t(1,lm1) - t(1,l-2) ) / t(1,lm1) ) <= aittol ) go to 160
 
  if ( .not. reglar .and. l == 4 ) go to 30
 
120 continue
 
  if ( errer <= ergoal ) go to 310
  go to 380

130 continue

  if ( .not. h2conv ) then
    aitken = .false.
    h2conv = .true.
  end if

140 continue

  fextrp = 4.0D+00

150 continue

  it = it+1
  vint = step * t(l,it)
  errer = abs ( step / (fextrp-1.0) * t(it-1,l))
  if ( errer <= ergoal ) go to 340
  if ( it == lm1 ) go to 270
  if ( t(it,lm1) == 0.0 ) go to 150
  if ( t(it,lm1) <= fextrp ) go to 270

  if ( abs ( t(it,lm1) / 4.0 - fextrp ) / fextrp < aittol ) then
    fextrp = fextrp*4.0D+00
  end if

  go to 150
 
160 continue

  if ( t(1,lm1) < aitlow ) then
    go to 380
  end if
 
  if ( .not. aitken ) then
    h2conv = .false.
    aitken = .true.
  end if
 
170 continue

  fextrp = t(l-2,lm1)
  if ( fextrp > 4.5 ) go to 140
  if ( fextrp < aitlow ) go to 380

  if ( abs ( fextrp - t(l-3,lm1) ) / t(1,lm1) > h2tol ) then
    go to 380
  end if

  sing = fextrp
  fextm1 = fextrp - 1.0D+00

  ait(1) = 0.0D+00
  do i = 2, l
    ait(i) = t(i,1) + (t(i,1)-t(i-1,1)) / fextm1
    r(i) = t(1,i-1)
    dif(i) = ait(i) - ait(i-1)
  end do

  it = 2

190 continue

  vint = step*ait(l)

200 continue

  errer = errer / fextm1
 
  if ( errer <= ergoal ) then
    ind = max ( ind, 2 )
    go to 340
  end if
 
210 continue

  it = it+1
  if ( it == lm1 ) go to 270

  if ( it <= 3 ) then
    h2next = 4.0D+00
    singnx = 2.0D+00 * sing
  end if

  if ( h2next < singnx ) go to 230
  fextrp = singnx
  singnx = 2.0D+00 * singnx
  go to 240

230 continue

  fextrp = h2next
  h2next = 4.0D+00 * h2next

240 continue
 
  do i = it, lm1
    if ( abs ( dif(i+1) ) > tabtlm ) then
      r(i+1) = dif(i) / dif(i+1)
    else
      r(i+1) = 0.0D+00
    end if
  end do
 
  h2tfex = -h2tol*fextrp
  if ( r(l) - fextrp < h2tfex ) go to 270
  if ( r(l-1) - fextrp < h2tfex ) go to 270
  errer = astep * abs ( dif(l) )
  fextm1 = fextrp - 1.0D+00
  do i = it, l
    ait(i) = ait(i)+dif(i) / fextm1
    dif(i) = ait(i)-ait(i-1)
  end do
 
  go to 190
 
270 continue

  fextrp = max(prever/errer,aitlow)
  prever = errer
  if (l<5) go to 40
  if (l-it>2.and.istage<mxstge) go to 370
  if (errer/fextrp**(maxtbl-l)<ergoal) go to 40
  go to 370
 
280 continue

  if ( errer > ergoal ) go to 370
  diff = abs ( t(1,l) ) * 2.0D+00 * fn
  go to 340
 
290 continue

  slope = (fend-fbeg) * 2.0D+00
  fbeg2 = fbeg * 2.0D+00
 
  do i = 1, 4
    diff = abs ( func ( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
    if ( diff > tabtlm ) go to 330
  end do
 
  go to 340
 
310 continue

  slope = (fend-fbeg)*2.0D+00
  fbeg2 = fbeg*2.0D+00
  i = 1
 
320 continue

  diff = abs ( func(beg+rn(i)*step) - fbeg2 - rn(i) * slope )
 
330 continue

  errer = max ( errer, astep * diff )
  if (errer > ergoal) go to 380
  i = i+1
  if ( i <= 4 ) go to 320
  ind = 3
 
340 continue

  result = result+vint
  error = error+errer
 
350 continue

  if (right) go to 360
  istage = istage-1
  if (istage==0) return
  reglar = reglsv(istage)
  beg = begin(istage)
  end = finis(istage)
  curest = curest-est(istage+1)+vint
  iend = ibeg-1
  fend = ts(iend)
  ibeg = ibegs(istage)
  go to 400
 
360 continue

  curest = curest+vint
  stage = stage*2.0D+00
  iend = ibeg
  ibeg = ibegs(istage)
  end = beg
  beg = begin(istage)
  fend = fbeg
  fbeg = ts(ibeg)
  go to 10
 
370 continue

  reglar = .true.
 
380 continue
 
  if ( istage == mxstge ) then
    ind = 5
    result = curest+vint
    return
  end if
 
390 continue

  if (right) go to 410
  reglsv(istage+1) = reglar
  begin(istage) = beg
  ibegs(istage) = ibeg
  stage = stage/2.0D+00

400 continue

  right = .true.
  beg = (beg+end)/2.0D+00
  ibeg = (ibeg+iend)/2
  ts(ibeg) = ts(ibeg) / 2.0D+00
  fbeg = ts(ibeg)
  go to 20

410 continue

  nnleft = ibeg-ibegs(istage)
  if (end+nnleft>=maxts) go to 440
  iii = ibegs(istage)
  ii = iend
  do i = iii, ibeg
    ii = ii+1
    ts(ii) = ts(i)
  end do
 
  do i = ibeg, ii
    ts(iii) = ts(i)
    iii = iii+1
  end do
 
  iend = iend+1
  ibeg = iend-nnleft
  fend = fbeg
  fbeg = ts(ibeg)
  finis(istage) = end
  end = beg
  beg = begin(istage)
  begin(istage) = end
  reglsv(istage) = reglar
  istage = istage+1
  reglar = reglsv(istage)
  est(istage) = vint
  curest = curest+est(istage)
  go to 10

440 continue

  ind = 4

460 continue

  result = curest+vint

  return
end
subroutine chinsp ( func, a, b, epsin, epsout, result )
!
!***********************************************************************
!
!! CHINSP estimates an integral using a modified Clenshaw-Curtis scheme.
!
!
!  Discussion:
!
!    The integral is approximated by Chebyshev polyonomials over each
!    subinterval.  These are integrated to give the approximate integral.
!    If the error estimate is unsatisfactory, the integration is repeated
!    with smaller intervals.
!
!    The internal parameter NUPPER is currently set to 9,
!    corresponding to 1024 subintervals for the unfolded integral,
!    and 1025 function evaluations.  This parameter may be changed
!    if necessary.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    T Havie
!    BIT 9 (1969), pages 338-350.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!      FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real*8 A, the lower limit of integration.
!
!    Input, real*8 B, the upper limit of integration.
!
!    Input, real*8 EPSIN, the relative error tolerance.
!
!    Output, real*8 EPSOUT, estimated integration error.
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer, parameter :: nupper = 9
!
  real*8 a
  real*8 a0
  real*8 a1
  real*8 a2
  real*8 acof(257)
  real*8 alf
  real*8 alfnj
  real*8 alfno
  real*8 b
  real*8 bcof(257)
  real*8 bet
  real*8 betnj
  real*8 betno
  real*8 bounds
  real*8 ccof(513)
  real*8 cof
  real*8 cofmax
  real*8 const1
  real*8 const2
  real*8 deln
  real*8 deltan
  real*8 epsin
  real*8 epsout
  real*8 error
  real*8 etank
  real*8, external :: func
  real*8 gamman
  real*8 hnstep
  integer i
  integer index
  integer j
  integer k
  integer ksign
  integer n
  integer ncof
  integer nhalf
  integer nn
  real*8, parameter :: one = 1.0D+00
  real*8 r1
  real*8 r2
  real*8 result
  real*8 rk
  real*8 rn
  real*8 rnderr
  real*8 rounde
  real*8 tend
  real*8 tnew
  real*8 triarg
  real*8 umid
  real*8 wmean
  real*8 xmin
  real*8 xplus
  real*8 xsink
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
!
!  ROUNDE = RNDERR*(R1+R2*N), where R1, R2 are two empirical constants.
!
!  Set coefficients in formula for accumulated roundoff error.
!  N is the current number of function values used.
!
  rnderr = epsilon ( 1.0D+00 )
 
  r1 = 1.0D+00
  r2 = 2.0D+00
  error = epsin
!
!  Integration interval parameters.
!
  alf = 0.5D+00 * ( b - a )
  bet = 0.5D+00 * ( b + a )
!
!  Parameters for trigonometric recurrence relations.
!
  triarg = atan ( 1.0D+00 )
  alfno = -1.0D+00
!
!  Parameters for integration stepsize and loops.
!
  rn = 2.0D+00
  n = 2
  nhalf = 1
  hnstep = 1.0D+00
!
!  Initial calculation for the end-point approximation.
!
  const1 = 0.5D+00 * ( func(a) + func(b) )
  const2 = func(bet)
  acof(1) = 0.5D+00 * (const1+const2)
  acof(2) = 0.5D+00 * (const1-const2)
  bcof(2) = acof(2)
  tend = 2.0D+00 * ( acof(1) - acof(2) / 3.0D+00 )
!
!  Start actual calculations.
!
  do i = 1, nupper
!
!  Compute function values.
!
    const1 = -sin(triarg)
    const2 = 0.5D+00 * alfno / const1
    alfno = const1
    betno = const2
    gamman = 1.0D+00 - 2.0D+00 * alfno**2
    deltan = -2.0D+00 * alfno * betno
    bcof(1) = 0.0D+00
 
    do j = 1, nhalf
      alfnj = gamman * const1 + deltan*const2
      betnj = gamman * const2 - deltan*const1
      xplus = alf * alfnj+bet
      xmin = -alf * alfnj+bet
      ccof(j) = func(xplus) + func(xmin)
      bcof(1) = bcof(1) + ccof(j)
      const1 = alfnj
      const2 = betnj
    end do
 
    bcof(1) = 0.5D+00 * hnstep * bcof(1)
!
!  Calculation of first B-coefficient finished compute the higher
!  coefficients if NHALF greater than one.
!
    if ( nhalf <= 1 ) go to 60
    const1 = one
    const2 = 0.0D+00
    ncof = nhalf-1
    ksign = -1
 
    do k = 1, ncof
!
!  Compute trigonometric sum for B-coefficient.
!
      etank = gamman * const1 - deltan*const2
      xsink = gamman * const2 + deltan*const1
      cof = 2.0D+00 * ( 2.0D+00 * etank**2 - 1.0D+00 )
      a2 = 0.0D+00
      a1 = 0.0D+00
      a0 = ccof(nhalf)
 
      do j = 1, ncof
        a2 = a1
        a1 = a0
        index = nhalf-j
        a0 = ccof(index) + cof * a1 - a2
      end do
 
      bcof(k+1) = hnstep * (a0-a1) * etank
      bcof(k+1) = ksign * bcof(k+1)
      ksign = -ksign
      const1 = etank
      const2 = xsink
 
    end do
!
!  Calculation of B-coefficients finished.
!
!  Compute new modified mid-point approximation when the interval
!  of integration is divided in N equal sub intervals.
!
60  continue
 
    umid = 0.0D+00
    rk = rn
    nn = nhalf+1
    do k = 1, nn
      index = nn+1-k
      umid = umid+bcof(index)/(rk**2-one)
      rk = rk-2.0D+00
    end do
 
    umid = -2.0D+00 * umid
!
!  Compute new C-coefficients for end-point approximation and largest
!  absolute value of coefficients.
!
    nn = n+2
    cofmax = 0.0D+00
 
    do j = 1, nhalf
      index = nn-j
      ccof(j) = 0.5D+00 * (acof(j)+bcof(j))
      ccof(index) = 0.5D+00 * (acof(j)-bcof(j))
      const1 = abs ( ccof(j) )
      cofmax = max ( cofmax, const1 )
      const1 = abs ( ccof(index) )
      cofmax = max ( cofmax, const1 )
    end do
 
    ccof(nhalf+1) = acof(nhalf+1)
!
!  Calculation of new coefficients finished.
!
!  Compute new end-point approximation when the interval of
!  integration is divided in 2N equal sub intervals.
!
    wmean = 0.5D+00 * (tend+umid)
    bounds = 0.5D+00 * (tend-umid)
    deln = 0.0D+00
    rk = 2.0D+00 * rn
    do j = 1, nhalf
      index = n+2-j
      deln = deln+ccof(index) / (rk**2-one)
      rk = rk-2.0D+00
    end do
 
    deln = -2.0D+00 * deln
    tnew = wmean+deln
    epsout = abs ( bounds / tnew )

    if ( cofmax < rnderr ) then
      go to 160
    end if

    rounde = rnderr*(r1+r2*rn)
    if ( epsout < rounde ) epsout = rounde
    if ( error < rounde ) error = rounde
    if ( epsout > error ) go to 160
!
!  Required accuracy obtained or the maximum number of function
!  values used without obtaining the required accuracy.
!
120 continue
 
    n = 2*n+1
    tend = alf*(tend+deln)
    umid = alf*(umid+deln)
    deln = alf*deln
    result = alf*tnew
    return
!
!  If I = NUPPER then the required accuracy is not obtained.
!
160 continue
 
    if ( i == nupper ) go to 120
 
    acof(1:n) = ccof(1:n)
    acof(n+1) = ccof(n+1)
    bcof(n+1) = ccof(n+1)
    tend = tnew
    nhalf = n
    n = 2 * n
    rn = 2.0D+00 * rn
    hnstep = 0.5D+00 * hnstep
    triarg = 0.5D+00 * triarg
 
  end do
 
  return
end
subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )
!
!***********************************************************************
!
!! CSPINT estimates the integral of a tabulated function.
!
!
!  Discussion:
!
!    The routine is given the value of a function F(X) at a set of 
!    nodes XTAB, and estimates
!
!      INTEGRAL (A to B) F(X) DX
!
!    by computing the cubic natural spline S(X) that interpolates
!    F(X) at the nodes, and then computing
!
!      INTEGRAL (A to B) S(X) DX
!
!    exactly.
!
!    Other output from the program includes the definite integral
!    from X(1) to X(I) of S(X), and the coefficients necessary for
!    the user to evaluate the spline S(X) at any point.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), contains the tabulated values of
!    the function, FTAB(I) = F(XTAB(I)).
!
!    Input, real*8 XTAB(NTAB), contains the points at which the
!    function was evaluated.  The XTAB's must be distinct and
!    in ascending order.
!
!    Input, integer NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real*8 A, lower limit of integration.
!
!    Input, real*8 B, upper limit of integration.
!
!    Output, real*8 Y(3,NTAB), will contain the coefficients
!    of the interpolating natural spline over each subinterval.
!
!    For XTAB(I) <= X <= XTAB(I+1),
!
!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!                   + Y(2,I)*(X-XTAB(I))**2
!                   + Y(3,I)*(X-XTAB(I))**3
!
!    Output, real*8 E(NTAB), E(I) = the definite integral from
!    XTAB(1) to XTAB(I) of S(X).
!
!    Workspace, real*8 WORK(NTAB).
!
!    Output, real*8 RESULT, the estimated value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer ntab
!
  real*8 a
  real*8 b
  real*8 e(ntab)
  real*8 ftab(ntab)
  integer i
  integer j
  real*8 r
  real*8 result
  real*8 s
  real*8 term
  real*8 u
  real*8 work(ntab)
  real*8 xtab(ntab)
  real*8 y(3,ntab)
!
  if ( ntab < 3 ) then
    write ( 6, '(a)' ) ' '
    write ( 6, '(a)' ) 'CSPINT - Fatal error!'
    write ( 6, '(a,i6)' ) '  NTAB must be at least 3, but input NTAB = ',ntab
!    stop  'Error in Cubic Spline Integration'
    result=0
    return
  end if
 
  do i = 1, ntab-1
 
    if ( xtab(i+1) <= xtab(i) ) then
      write ( 6, '(a)' ) ' '
      write ( 6, '(a)' ) 'CSPINT - Fatal error!'
      write ( 6, '(a)' ) '  Nodes not in strict increasing order.'
      write ( 6, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
!      write ( 6, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
!      write ( 6, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
      write ( 6, *) 'Npoints ',ntab
        do j=1,ntab
                write(6,*)'Input #',j,xtab(j),ftab(j)
        enddo
      result=0
      return
!      stop 'Error in Cubic Spline Integration'
    end if
 
  end do
 
  s = 0.0D+00
  do i = 1, ntab-1
    r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
    y(2,i) = r - s
    s = r
  end do
 
  result = 0.0D+00
  s = 0.0D+00
  r = 0.0D+00
  y(2,1) = 0.0D+00
  y(2,ntab) = 0.0D+00
 
  do i = 2, ntab-1
    y(2,i) = y(2,i)+r*y(2,i-1)
    work(i) = 2.0D+00 * ( xtab(i-1) - xtab(i+1) ) - r * s
    s = xtab(i+1) - xtab(i)
    r = s / work(i)
  end do
 
  do j = 2, ntab-1
    i = ntab+1-j
    y(2,i) = ((xtab(i+1)-xtab(i))*y(2,i+1)-y(2,i)) / work(i)
  end do
 
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    r = y(2,i+1)-y(2,i)
    y(3,i) = r / s
    y(2,i) = 3.0D+00 * y(2,i)
    y(1,i) = (ftab(i+1)-ftab(i)) / s-(y(2,i)+r)*s
  end do
 
  e(1) = 0.0D+00
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    term = (((y(3,i)* 0.25D+00 *s+y(2,i) / 3.0 ) *s+y(1,i)* 0.5 )*s+ftab(i))*s
    e(i+1) = e(i) + term
  end do
!
!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!
  r = a
  u = 1.0D+00
 
  do j = 1, 2
!
!  The endpoint is less than or equal to XTAB(1).
!
    if ( r <= xtab(1) ) then
      result = result-u*((r-xtab(1))*y(1,1)*0.5D+00 +ftab(1))*(r-xtab(1))
!
!  The endpoint is greater than or equal to XTAB(NTAB).
!
    else if ( r >= xtab(ntab) ) then

      result = result-u*(e(ntab)+(r-xtab(ntab))*(ftab(ntab)+ &
        0.5D+00 *(ftab(ntab-1)+(xtab(ntab)-xtab(ntab-1))*y(1,ntab-1)) &
        *(r-xtab(ntab))))
!
!  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!
    else
      do i = 1,ntab-1
 
        if ( r <= xtab(i+1) ) then
          r = r-xtab(i)
          result = result-u*(e(i)+(((y(3,i)*0.25*r+y(2,i)/3.0)*r &
            +y(1,i)*0.5D+00 )*r+ftab(i))*r)
          go to 120
        end if
 
      end do
 
    end if
 
  120   continue
 
    u = -1.0D+00
    r = b
 
  end do
 
  return
end
subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )
!
!***********************************************************************
!
!! CUBINT approximates an integral using cubic interpolation of data.
!
!
!  Discussion:
!
!    The integral to be approximated is
! 
!      INTEGRAL (XTAB(IB) to XTAB(IA)) F(X) DX
!
!    The routine estimates the error in integration.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    P E Gill and G F Miller
!    An algorithm for the integration of unequally spaced data,
!    Comput J, Number 15, 1972, pages 80-83.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, real*8 XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, integer NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, integer IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real*8 RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real*8 ERROR, an estimate of the error in
!    integration.
!
  implicit real*8 (a-h,o-z)
!
  integer ntab
!
  real*8 c
  real*8 d1
  real*8 d2
  real*8 d3
  real*8 error
  real*8 ftab(ntab)
  real*8 h1
  real*8 h2
  real*8 h3
  real*8 h4
  integer i
  integer ia
  integer ib
  integer ind
  integer it
  integer j
  integer k
  real*8 r1
  real*8 r2
  real*8 r3
  real*8 r4
  real*8 result
  real*8 s
  real*8 term
  real*8 xtab(ntab)
!
  result = 0.0D+00
  error = 0.0D+00
 
  if ( ia == ib ) then
    return
  end if
 
  if ( ntab < 4 ) then
    write ( 6, '(a)' ) ' '
    write ( 6, '(a)' ) 'CUBINT - Fatal error!'
    write ( 6, '(a,i6)' ) '  NTAB must be at least 4, but input NTAB = ',ntab
    result=0
    return
!    stop
  end if
 
  if ( ia < 1 ) then
    write ( 6, '(a)' ) ' '
    write ( 6, '(a)' ) 'CUBINT - Fatal error!'
    write ( 6, '(a,i6)' ) '  IA must be at least 1, but input IA = ',ia
    result=0
    return
!    stop
  end if
 
  if ( ia > ntab ) then
    write ( 6, '(a)' ) ' '
    write ( 6, '(a)' ) 'CUBINT - Fatal error!'
    write ( 6, '(a,i6)' ) '  IA must be <= NTAB, but input IA=',ia
    result=0
    return
!    stop
  end if
 
  if ( ib < 1 ) then
    write ( 6, '(a)' ) ' '
    write ( 6, '(a)' ) 'CUBINT - Fatal error!'
    write ( 6, '(a,i6)' ) '  IB must be at least 1, but input IB = ',ib
    result=0
    return
!    stop
  end if
 
  if ( ib > ntab ) then
    write ( 6, '(a)' ) ' '
    write ( 6, '(a)' ) 'CUBINT - Fatal error!'
    write ( 6, '(a,i6)' ) '  IB must be <= NTAB, but input IB=',ib
    result=0
    return
!    stop
  end if
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  if ( ia > ib ) then
    ind = -1
    it = ib
    ib = ia
    ia = it
  else
    ind = 1
  end if
 
  s = 0.0D+00
  c = 0.0D+00
  r4 = 0.0D+00
  j = ntab-2
  if ( ia < ntab-1 .or. ntab == 4 ) then
    j=max(3,ia)
  end if

  k = 4
  if ( ib > 2 .or. ntab == 4 ) then
    k=min(ntab,ib+2)-1
  end if
 
  do i = j, k
 
    if ( i <= j ) then
 
      h2 = xtab(j-1)-xtab(j-2)
      d3 = (ftab(j-1)-ftab(j-2)) / h2
      h3 = xtab(j)-xtab(j-1)
      d1 = (ftab(j)-ftab(j-1)) / h3
      h1 = h2+h3
      d2 = (d1-d3)/h1
      h4 = xtab(j+1)-xtab(j)
      r1 = (ftab(j+1)-ftab(j)) / h4
      r2 = (r1-d1) / (h4+h3)
      h1 = h1+h4
      r3 = (r2-d2) / h1
 
      if ( ia <= 1 ) then
        result = h2 * (ftab(1)+h2*(0.5*d3-h2*(d2/6.0-(h2+h3+h3)*r3/12.)))
        s = -h2**3 * (h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0
      end if
 
    else
 
      h4 = xtab(i+1)-xtab(i)
      r1 = (ftab(i+1)-ftab(i))/h4
      r4 = h4+h3
      r2 = (r1-d1)/r4
      r4 = r4+h2
      r3 = (r2-d2)/r4
      r4 = (r3-d3)/(r4+h1)
 
    end if
 
    if ( i > ia .and. i <= ib ) then
 
      term = h3*((ftab(i)+ftab(i-1))*0.5-h3*h3*(d2+r2+(h2-h4)*r3) / 12.0 )
      result = result+term
      c = h3**3*(2.0D+00 *h3*h3+5.*(h3*(h4+h2) + 2.0 * h2 * h4 ) ) / 120.0D+00
      error = error+(c+s)*r4
 
      if ( i /= j ) then
        s = c
      else
        s = s+c+c
      end if
 
    else
 
      error = error+r4*s
 
    end if
 
    if ( i >= k ) then
 
      if ( ib >= ntab ) then
        term = h4*(ftab(ntab) - h4*(0.5*r1+h4*(r2/6.0 +(h3+h3+h4)*r3/12.)))
        result = result + term
        error = error - h4**3 * r4 * &
          ( h4 * ( 3.0 * h4 + 5.0 * h2 ) &
          + 10.0 * h3 * ( h2 + h3 + h4 ) ) / 60.0D+00
      end if
 
      if ( ib >= ntab-1 ) error=error+s*r4
    else
      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    end if
 
  end do
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  if ( ind /= 1 ) then
    it = ib
    ib = ia
    ia = it
    result = -result
    error = -error
  end if
 
  return
end
subroutine filon_cos ( ftab, a, b, ntab, t, result )
!
!***********************************************************************
!
!! FILON_COS uses Filon's method on integrals with a cosine factor.
!
!
!  Discussion:
!
!    The integral to be approximated has the form:
!
!      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
!
!    where T is user specified.
!
!    The function is interpolated over each subinterval by
!    a parabolic arc.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    pages 890-891.
!
!    S M Chase and L D Fosdick,
!    Algorithm 353, Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, 1969, pages 457-458.
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    19 February 2002
!
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), contains the value of the function
!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!
!    Input, real*8 A, B, the limits of integration.
!
!    Input, integer NTAB, the number of data points.
!    NTAB must be odd, and greater than 1.
!
!    Input, real*8 T, the multiplier of the X argument of the cosine.
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer ntab
!
  real*8 a
  real*8 alpha
  real*8 b
  real*8 beta
  real*8 c2n
  real*8 c2nm1
  real*8 cost
  real*8 ftab(ntab)
  real*8 gamma
  real*8 h
  real*8 result
  real*8 sint
  real*8 t
  real*8 theta
  real*8 xtab(ntab)
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( ntab <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
    write ( *, '(a)' ) '  NTAB < 2'
    write ( *, '(a,i6)' ) '  NTAB = ', ntab
    stop
  end if
 
  if ( mod ( ntab, 2 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
    write ( *, '(a)' ) '  NTAB must be odd.'
    write ( *, '(a,i6)' ) '  NTAB = ', ntab
    stop
  end if
!
!  Set up a vector of the NTAB X values.
! 
  call rvec_even ( a, b, ntab, xtab )

  h = ( b - a ) / dble ( ntab - 1 )
  theta = t * h
  sint = sin ( theta )
  cost = cos ( theta )

  alpha = ( theta**2 + theta * sint * cost &
    - 2.0D+00 * sint**2 ) / theta**3

  beta = ( 2.0D+00 * theta + 2.0D+00 * theta * cost**2 &
    - 4.0D+00 * sint * cost ) / theta**3

  gamma = 4.0D+00 * ( sint - theta * cost ) / theta**3
  
  c2n = sum ( ftab(1:ntab:2) * cos ( t * xtab(1:ntab:2) ) ) &
    - 0.5D+00 * ( ftab(ntab) * cos ( t * xtab(ntab) ) &
                + ftab(1) * cos ( t * xtab(1) ) )

  c2nm1 = sum ( ftab(2:ntab-1:2) * cos ( t * xtab(2:ntab-1:2) ) )
 
  result = h * ( &
      alpha * ( ftab(ntab) * sin ( t * xtab(ntab) ) & 
              - ftab(1) * sin ( t * xtab(1) ) ) &
    + beta * c2n &
    + gamma * c2nm1 )

  return
end
subroutine filon_sin ( ftab, a, b, ntab, t, result )
!
!***********************************************************************
!
!! FILON_SIN uses Filon's method on integrals with a sine factor.
!
!
!  Discussion:
!
!    The integral to be approximated has the form
!
!      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
!
!    where T is user specified.
!
!    The function is interpolated over each subinterval by
!    a parabolic arc.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    pages 890-891.
!
!    S M Chase and L D Fosdick,
!    Algorithm 353, Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, 1969, pages 457-458.
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    19 February 2002
!
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), contains the value of the function
!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!
!    Input, real*8 A, B, the limits of integration.
!
!    Input, integer NTAB, the number of data points, including the
!    endpoints.  NTAB must be odd, and greater than 1.
!
!    Input, real*8 T, multiplier of the X argument of the sine.
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer ntab
!
  real*8 a
  real*8 alpha
  real*8 b
  real*8 beta
  real*8 cost
  real*8 ftab(ntab)
  real*8 gamma
  real*8 h
  real*8 result
  real*8 s2n
  real*8 s2nm1
  real*8 sint
  real*8 t
  real*8 theta
  real*8 xtab(ntab)
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( ntab <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
    write ( *, '(a)' ) '  NTAB < 2'
    write ( *, '(a,i6)' ) '  NTAB = ',ntab
    stop
  end if
 
  if ( mod ( ntab, 2 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
    write ( *, '(a)' ) '  NTAB must be odd.'
    write ( *, '(a,i6)' ) '  NTAB = ',ntab
    stop
  end if
!
!  Set up a vector of the NTAB X values.
! 
  call rvec_even ( a, b, ntab, xtab )

  h = ( b - a ) / dble ( ntab - 1 )
  theta = t * h
  sint = sin ( theta )
  cost = cos ( theta )
 
  alpha = ( theta**2 + theta * sint * cost &
    - 2.0D+00 * sint**2 ) / theta**3

  beta = ( 2.0D+00 * theta + 2.0D+00 * theta * cost**2 &
    - 4.0D+00 * sint * cost ) / theta**3

  gamma = 4.0D+00 * ( sint - theta * cost ) / theta**3
   
  s2n = sum ( ftab(1:ntab:2) * sin ( t * xtab(1:ntab:2) ) ) &
    - 0.5D+00 * ( ftab(ntab) * sin ( t * xtab(ntab) ) &
                + ftab(1) * sin ( t * xtab(1) ) )

  s2nm1 = sum ( ftab(2:ntab-1:2) * sin ( t * xtab(2:ntab-1:2) ) )

  result = h * ( &
      alpha * ( ftab(1) * cos ( t * xtab(1) ) &
              - ftab(ntab) * cos ( t * xtab(ntab) ) ) &
    + beta * s2n &
    + gamma * s2nm1 )
 
  return
end
real*8 function gamma ( x )
!
!*******************************************************************************
!
!! GAMMA calculates the Gamma function for a real*8 argument X.
!
!
!  Definition:
!
!    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
!
!  Recursion:
!
!    GAMMA(X+1) = X * GAMMA(X)
!
!  Special values:
!
!    GAMMA(0.5) = SQRT(PI)
!    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for X .GE. 12 are from reference 2.
!    The accuracy achieved depends on the arithmetic system, the
!    compiler, the intrinsic functions, and proper selection of the
!    machine-dependent constants.
!
!  Machine-dependent constants:
!
!    BETA: radix for the floating-point representation.
!    MAXEXP: the smallest positive power of BETA that overflows.
!    XBIG: the largest argument for which GAMMA(X) is representable
!      in the machine, i.e., the solution to the equation
!      GAMMA(XBIG) = BETA**MAXEXP.
!    XINF: the largest machine representable floating-point number;
!      approximately BETA**MAXEXP.
!    EPS: the smallest positive floating-point number such that
!      1.0+EPS .GT. 1.0.
!    XMININ: the smallest positive floating-point number such that
!      1/XMININ is machine representable.
!
!    Approximate values for some important machines are:
!
!                               BETA       MAXEXP        XBIG
!
!    CRAY-1         (S.P.)        2         8191        966.961
!    Cyber 180/855
!      under NOS    (S.P.)        2         1070        177.803
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)        2          128        35.040
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)        2         1024        171.624
!    IBM 3033       (D.P.)       16           63        57.574
!    VAX D-Format   (D.P.)        2          127        34.844
!    VAX G-Format   (D.P.)        2         1023        171.489
!
!                               XINF         EPS        XMININ
!
!    CRAY-1         (S.P.)   5.45d+2465   7.11d-15    1.84d-2466
!    Cyber 180/855
!      under NOS    (S.P.)   1.26d+322    3.55d-15    3.14d-294
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)   3.40d+38     1.19d-7     1.18d-38
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
!    IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
!    VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
!    VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!  Reference:
!
!    W J Cody,
!    "An Overview of Software Development for Special Functions",
!    Lecture Notes in Mathematics, 506,
!    Numerical Analysis Dundee, 1975,
!    G. A. Watson (ed.),
!    Springer Verlag, Berlin, 1976.
!
!    Hart et al,
!    Computer Approximations,
!    Wiley and sons, New York, 1968.
!
!  Author:
!
!    W. J. Cody and L. Stoltz,
!    Applied Mathematics Division,
!    Argonne National Laboratory,
!    Argonne, Illinois, 60439.
!
!  Parameters:
!
!    Input, real*8 X, the argument of the function.
!
!    Output, real*8 GAMMA, the value of the function.  The program
!    returns the value XINF for singularities or when overflow would occur.
!    The computation is believed to be free of underflow and overflow.
!
  implicit real*8 (a-h,o-z)
!
  real*8, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728d-03, &
     8.4171387781295d-04, &
    -5.952379913043012d-04, &
     7.93650793500350248d-04, &
    -2.777777777777681622553d-03, &
     8.333333333333333331554247d-02, &
     5.7083835261d-03 /)
  real*8, parameter :: EPS = 1.19d-07
  real*8 fact
!  real*8 gamma
  integer i
  integer n
  real*8, parameter, dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811d+00, &
     2.47656508055759199108314d+01, &
    -3.79804256470945635097577d+02, &
     6.29331155312818442661052d+02, &
     8.66966202790413211295064d+02, &
    -3.14512729688483675254357d+04, &
    -3.61444134186911729807069d+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real*8, parameter :: PI = &
    3.14159265358979323846264338327950288419716939937510D+00
  real*8, parameter, dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353d+01, &
     3.15350626979604161529144d+02, &
    -1.01515636749021914166146d+03, &
    -3.10777167157231109440444d+03, &
     2.25381184209801510330112d+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456d+05, &
    -1.15132259675553483497211d+05 /)
  real*8, parameter :: SQRTPI = 0.9189385332046727417803297d+00
  real*8 sum1
  real*8 x
  real*8, parameter :: XBIG = 35.040D+00
  real*8 xden
  real*8, parameter :: XINF = 3.4d+38
  real*8, parameter :: XMININ = 1.18d-38
  real*8 xnum
  real*8 y
  real*8 y1
  real*8 ysq
  real*8 z
!
  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    gamma = y - y1

    if ( gamma /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - PI / sin ( PI * gamma )
      y = y + 1.0D+00

    else

      gamma = XINF
      return

    end if

  end if
!
!  Argument < EPS
!
  if ( y < EPS ) then

    if (y >= XMININ ) then
      gamma = 1.0D+00 / y
    else
      gamma = XINF
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0D+00 < argument < 1.0D+00
!
    if ( y < 1.0D+00 ) then
      z = y
      y = y + 1.0D+00
!
!  1.0D+00 < argument < 12.0, reduce argument if necessary.
!
    else
      n = int ( y ) - 1
      y = y - dble ( n )
      z = y - 1.0D+00
    end if
!
!  Evaluate approximation for 1.0D+00 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    gamma = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0D+00 < argument < 1.0.
!
    if ( y1 < y ) then
      gamma = gamma / y1
!
!  Adjust result for case  2.0D+00 < argument < 12.0.
!
    else if ( y1 > y ) then

      do i = 1, n
        gamma = gamma * y
        y = y + 1.0D+00
      end do

    end if
!
!  Evaluate for 12 <= argument.
!
  else

    if ( y <= XBIG ) then

      ysq = y**2
      sum1 = c(7)
      do i = 1, 6
        sum1 = sum1 / ysq + c(i)
      end do
      sum1 = sum1 / y - y + SQRTPI
      sum1 = sum1 + ( y - 0.5D+00 ) * log ( y )
      gamma = exp ( sum1 )

    else

      gamma = XINF
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    gamma = - gamma
  end if

  if ( fact /= 1.0D+00 ) then
    gamma = fact / gamma
  end if

  return
end
subroutine gaus8 ( func, a, b, err, result, ierr )
!
!***********************************************************************
!
!! GAUS8 estimates the integral of a function.
!
!
!  Discussion:
!
!    GAUS8 integrates real*8 functions of one variable over finite
!    intervals using an adaptive 8-point Legendre-Gauss
!    algorithm.
!
!    GAUS8 is intended primarily for high accuracy integration or
!    integration of smooth functions.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Author:
!
!    R E Jones,
!    Sandia National Laboratory,
!    Los Alamos, New Mexico
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC, name of external function to be integrated.
!    This name must be in an external statement in the calling program.
!    FUNC must be a function of one real*8 argument.  The value
!    of the argument to FUNC is the variable of integration
!    which ranges from A to B.
!
!    Input, real*8 A, the lower limit of integration.
!
!    Input, real*8 B, the upper limit of integration.
!
!    Input/output, real*8 ERR.
!    On input, ERR is a requested pseudorelative error
!    tolerance.  Normally pick a value of ABS ( ERR ) so that
!    STOL < ABS ( ERR ) <= 1.0E-3 where STOL is the single precision
!    unit roundoff  = R1MACH(4).
!
!    RESULT will normally have no more error than
!    ABS ( ERR ) times the integral of the absolute value of
!    FUN(X).  Usually, smaller values for ERR yield more
!    accuracy and require more function evaluations.
!
!    A negative value for ERR causes an estimate of the
!    absolute error in RESULT to be returned in ERR.  Note that
!    ERR must be a variable (not a constant) in this case.
!    Note also that the user must reset the value of ERR
!    before making any more calls that use the variable ERR.
!
!    On output, ERR will be an estimate of the absolute error
!    in RESULT if the input value of ERR was negative.  ERR is
!    unchanged if the input value of ERR was non-negative.
!
!    The estimated error is solely for information to the user
!    and should not be used as a correction to the computed integral.
!
!    Output, real*8 RESULT, the computed value of the integral.
!
!    Output, integer IERR, a status code.
!
!    Normal Codes:
!
!     1 RESULT most likely meets requested error tolerance, or A = B.
!    -1 A and B are too nearly equal to allow normal
!        integration.  RESULT is set to zero.
!
!     Abnormal Code:
!
!     2 RESULT probably does not meet requested error tolerance.
!
  implicit real*8 (a-h,o-z)
!
  real*8 a
  real*8 aa(30)
  real*8 ae
  real*8 anib
  real*8 area
  real*8 b
  real*8 c
  real*8 ce
  real*8 ee
  real*8 ef
  real*8 eps
  real*8 err
  real*8 est
  real*8, external :: func
  real*8 g8
  real*8 gl
  real*8 glr
  real*8 gr(30)
  real*8 h
  real*8 hh(30)
  integer i1mach
  integer, save :: icall = 0
  integer ierr
  integer k
  integer, save :: kml = 6
  integer, save :: kmx = 5000
  integer l
  integer lmn
  integer lmx
  integer lr(30)
  integer mxl
  integer nbits
  integer nib
  integer, save :: nlmn = 1
  integer nlmx
  real*8 d1mach
  real*8 result
  real*8 tol
  real*8 vl(30)
  real*8 vr
  real*8, save :: w1
  real*8, save :: w2
  real*8, save :: w3
  real*8, save :: w4
  real*8 x
  real*8, save :: x1
  real*8, save :: x2
  real*8, save :: x3
  real*8, save :: x4
!
  data x1, x2, x3, x4/ &
          1.83434642495649805d-01,     5.25532409916328986d-01, &
          7.96666477413626740d-01,     9.60289856497536232d-01/
!
  data w1, w2, w3, w4/ &
          3.62683783378361983d-01,     3.13706645877887287d-01, &
          2.22381034453374471d-01,     1.01228536290376259d-01/
!
!  Warning!  Statement function!
!
  g8(x,h) = h*((w1*(func(x-x1*h) + func(x+x1*h)) &
             +w2*(func(x-x2*h) + func(x+x2*h))) &
            +(w3*(func(x-x3*h) + func(x+x3*h)) &
             +w4*(func(x-x4*h) + func(x+x4*h))))
!
  if ( a == b ) then
    err = 0.0D+00
    result = 0.0D+00
    return
  end if
 
  if ( icall /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAUS8 - Fatal error!'
    write ( *, '(a)' ) '  GAUS8 was called recursively.'
    stop
  end if

  icall = 1
  k = i1mach(11)
  anib = d1mach(5) * real(i1mach(11)) / 0.30102000d0
  nbits = int(anib)
  nlmx = min ( 30, (nbits*5)/8 )
  result = 0.0D+00
  ierr = 1
  ce = 0.0D+00
  result = 0.0D+00
  lmx = nlmx
  lmn = nlmn
 
  if ( b /= 0.0d0 ) then

    if ( sign ( 1.d0, b ) * a <= 0.0D+00 ) then
      go to 10
    end if

    c = abs ( 1.0d0 - a / b )
    if ( c > 0.1d0 ) go to 10
    if ( c <= 0.0D+00 ) go to 140
    anib = 0.5d0 - log(c) / 0.69314718d+00
    nib = int(anib)
    lmx = min(nlmx,nbits-nib-7)
    if ( lmx < 1 ) go to 130
    lmn = min ( lmn, lmx )

  end if
 
10    continue
 
  tol = max ( abs ( err ), 2.0d0**(5-nbits)) / 2.0D+00
  if ( err == 0.0d0 ) then
    tol=sqrt( epsilon ( 1.0D+00 ) )
  end if

  eps = tol
  hh(1) = (b-a) / 4.0D+00
  aa(1) = a
  lr(1) = 1
  l = 1
  est = g8 ( aa(l) + 2.0d0*hh(l),2.0d0*hh(l) )
  k = 8
  area = abs ( est )
  ef = 0.5d0
  mxl = 0
!
!  Compute refined estimates, estimate the error, etc.
!
20 continue
 
  gl = g8 ( aa(l)+hh(l),hh(l) )
  gr(l) = g8(aa(l)+3.0d0*hh(l),hh(l))
  k = k + 16
  area = area + ( abs ( gl ) + abs ( gr(l) ) - abs ( est ) )
 
  glr = gl + gr(l)
  ee = abs ( est - glr ) * ef
  ae = max ( eps * area, tol * abs ( glr ) )

  if ( ee - ae <= 0.0D+00 ) then
    go to 40
  else
    go to 50
  end if
 
30 continue
 
  mxl = 1
 
40 continue
 
  ce = ce + (est-glr)
 
  if ( lr(l) <= 0 ) then
    go to 60
  else
    go to 80
  end if
!
!  Consider the left half of this level
!
50 continue

  if ( k > kmx ) lmx = kml
  if ( l >= lmx ) go to 30
  l = l + 1
  eps = eps * 0.5D+00
  ef = ef / sqrt ( 2.0D+00 )
  hh(l) = hh(l-1) * 0.5D+00
  lr(l) = -1
  aa(l) = aa(l-1)
  est = gl
  go to 20
!
!  Proceed to right half at this level
!
60 continue

  vl(l) = glr

70 continue

  est = gr(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 4.0D+00 * hh(l)
  go to 20
!
!  Return one level
!
80 continue

  vr = glr

90 continue

  if ( l <= 1 ) go to 120
  l = l - 1
  eps = eps * 2.0D+00
  ef = ef * sqrt ( 2.0D+00 )
 
  if ( lr(l) <= 0 ) then
    vl(l) = vl(l+1) + vr
    go to 70
  else
    vr = vl(l+1) + vr
    go to 90
  end if
!
!  Exit
!
120   continue
 
  result = vr
  if ( mxl == 0 .or. abs ( ce ) <= 2.0D+00 * tol * area ) go to 140
  ierr = 2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAUS8 - Warning!'
  write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
  icall = 0

  if ( err < 0.0D+00 ) then
    err = ce
  end if

  return
 
130   continue
 
  ierr = -1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAUS8 - Warning!'
  write ( *, '(a)' ) '  A and B are too close to carry out integration.'
 
140   continue
 
  icall = 0
 
  if ( err < 0.0D+00 ) then
    err = ce
  end if
 
  return
end
subroutine hiordq ( n, y, delt, work, result )
!
!***********************************************************************
!
!! HIORDQ approximates the integral of a function using equally spaced data.
!
!
!  Discussion:
!
!    The method applies the trapezoidal rule to various subsets of the
!    data, and then applies Richardson extrapolation.
!
!  Author:
!
!    Alan Kaylor Cline,
!    Department of Computer Science,
!    University of Texas at Austin.
!
!  Modified:
!
!    19 February 2002
!
!  Parameters:
!
!    Input, integer N, number of data points.
!
!    Input, real*8 Y(N), the Y values of the data.
!
!    Input, real*8 DELT, the spacing between the X values of the
!    data.  The actual X values are not needed!
!
!    Work array, real*8 WORK(2*(N-1)).  The actual minimum amount
!    of workspace required is two times the number of integer
!    divisors of N-1.
!
!    Output, real*8 RESULT, the approximation to the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer n
!
  real*8 delt
  real*8 fac
  integer i
  integer j
  integer jbak
  integer jj
  integer k
  real*8 result
  real*8 sum2
  real*8 sum1
  real*8 work(2*(n-1))
  real*8 y(n)
!
!  Determine initial trapezoidal rule
!
  sum1 = ( y(1) + y(n) ) / 2.0D+00
  j = -1
 
  do k = 1, n-1
!
!  Check if K divides N-1
!
    if ( ((n-1)/k)*k == n-1 ) then
!
!  Determine the K-point trapezoidal rule.
!
      sum2 = -sum1
      do i = 1, n, (n-1)/k
        sum2 = sum2 + y(i)
      end do
 
      j = j + 2
      work(j) = delt * sum2 * dble ( ( n - 1 ) / k )
      work(j+1) = dble ( ((n-1)/k)**2 )
!
!  Apply Richardson extrapolation.
!
      if ( k /= 1 ) then
 
        do jj = 3, j, 2
          jbak = j+1-jj
          fac = work(j+1) / ( work(j+1) - work(jbak+1) )
          work(jbak) = work(jbak+2) + fac * ( work(jbak) - work(jbak+2) )
        end do
 
      end if
 
    end if
 
  end do
 
  result = work(1)
 
  return
end
subroutine iratex ( func, a, b, epsin, epsout, result, ind )
!
!***********************************************************************
!
!! IRATEX estimates the integral of a function.
!
!
!  Discussion:
!
!    IRATEX estimates the integral from A to B of F(X), using the
!    trapezoidal rule for several stepsizes H.
!
!    Then a rational function of H*H is fitted to these results, and 
!    an estimate of the integral is computed by extrapolating the 
!    results to a stepsize of zero.
!
!  Reference:
!
!    R Bulirsch and J Stoer,
!    Fehlerabschaetzungen und Extrapolation mit rationaled Funktionen
!      bei Verfahren vom Richardson-Typus,
!    (Error estimates and extrapolation with rational functions
!      in processes of Richardson type),
!    Numerische Mathematik,
!    Volume 6 (1964), pages 413-427.
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC.
!    FUNC is the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!      FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real*8 A, the lower limit of integration.
!
!    Input, real*8 B, the upper limit of integration.
!
!    Input, real*8 EPSIN, the requested relative error tolerance.
!
!    Output, real*8 EPSOUT, an estimate of the error in the integration.
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
!    Output, integer IND, error return flag.
!    IND = 0 if the requested accuracy was not achieved,
!    IND = 1 if the accuracy was achieved.
!
  implicit real*8 (a-h,o-z)
!
  real*8 a
  real*8 arg
  real*8 b
  real*8 ba
  logical bo
  logical bu
  real*8 c
  real*8 d(6)
  real*8 d1
  real*8 ddt
  real*8 den
  real*8 dt(7)
  real*8 e
  real*8 ent
  real*8 epsin
  real*8 epsout
  real*8, external :: func
  real*8 gr
  real*8 hm
  integer i
  integer ind
  integer m
  integer mr
  integer n
  integer np1
  logical odd
  real*8 rnderr
  real*8 result
  real*8 sm
  real*8 t
  real*8 t1
  real*8 t2
  real*8 t2a
  real*8 ta
  real*8 tab
  real*8 tb
  real*8 tnt
  real*8 v
  real*8 w
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  rnderr = epsilon ( 1.0D+00 )
  epsin = max ( epsin, 8.0 * rnderr )
  ind = 0
  n = 2
  np1 = 3
  ba = b-a
  t1 = 0.0D+00
  gr = 0.0D+00
  sm = 0.0D+00
  t2a = 0.5D+00 * ( func ( a ) + func ( b ) )
  t2 = t2a
  tb = abs ( t2a )
  c = t2 * ba

  dt(1) = c
  dt(2:7) = 0.0D+00
 
  odd = .true.
  bu = .false.
 
  do m = 1, 15
 
    bo = m>=7
    hm = ba / real(n)
!
!  N+1 is odd.
!
    if ( odd ) then
 
      do i = 1, n, 2
        arg = a + real(i) * hm
        w = func(arg)
        t2 = t2 + w
        tb = abs ( w ) + tb
      end do
 
      ent = t2
      tab = tb * abs ( hm )
      d(1) = 16.0D+00 / 9.0D+00
      d(3) = 4.0D+00 * d(1)
      d(5) = 4.0D+00 * d(3)
!
!  N+1 is even.
!
    else
 
      do i = 1, n, 6
        w = real(i) * hm
        t1 = t1 + func(a+w)+func(b-w)
      end do
 
      ent = t1+t2a
      t2a = t2
      d(1) = 2.25D+00
      d(3) = 9.0D+00
      d(5) = 36.0D+00
 
    end if
 
    ddt = dt(1)
    t = ent*hm
    dt(1) = t
    ent = t
 
    if ( m < 7 ) then
      mr = m
      w = dble ( n * n )
      d(m) = w
    else
      mr = 6
      d(6) = 64.0D+00
      w = 144.0D+00
    end if
 
    do i = 1, mr
 
      d1 = d(i) * ddt
      den = d1 - ent
      e = ent - ddt
      tnt = ent
      v = 0.0D+00
      ent = 0.0D+00
 
      if ( abs ( den ) >= epsin ) then
        e = e / den
        v = tnt * e
        ent = d1 * e
        t = v + t
        ddt = dt(i+1)
      end if
 
      dt(i+1) = v
 
    end do
 
    ta = c
    c = t
    result = c
    if ( .not. bo ) then
      t = t-v
    end if

    v = t-ta
    t = v+t
    epsout = abs ( v )
 
    if ( ta < t ) then
      d1 = ta
      ta = t
      t = d1
    end if
 
    bo = bo .or. ( ta < gr .and. t > sm )
 
    if ( bu .and. bo .and. epsout < epsin * w * tab ) then
      go to 140
    end if
 
    gr = ta + epsin
    sm = t - epsin
    odd = .not. odd
    n = np1
    np1 = n+1
    bu = bo
    d(2) = 4.0D+00
    d(4) = 16.0D+00
 
  end do
 
  bo = .false.
 
  140 continue
 
  v = rnderr*tab
 
  epsout = max ( epsout, v )

  if (bo) ind = 1

  return
end
real*8 function pi ( )
!
!*******************************************************************************
!
!! PI returns the value of pi.
!
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 PI, the value of pi.
!
  implicit real*8 (a-h,o-z)
!
!  real*8 pi
!
   pi = acos(-1.D0)
!  pi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
subroutine plint ( ftab, xtab, ntab, a, b, result )
!
!***********************************************************************
!
!! PLINT approximates the integral of unequally spaced data.
!
!
!  Discussion:
!
!    The method uses piecewise linear interpolation.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), the function values, FTAB(I) = F(XTAB(I)).
!
!    Input, real*8 XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, integer NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 2.
!
!    Input, real*8 A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real*8 B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer ntab
!
  real*8 a
  real*8 b
  real*8 fa
  real*8 fb
  real*8 ftab(ntab)
  integer i
  integer ihi
  integer ilo
  integer ind
  real*8 result
  real*8 slope
  real*8 syl
  real*8 xtab(ntab)
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( ntab < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLINT - Fatal error!'
    write ( *, '(a,i6)' ) '  NTAB < 2, NTAB = ',ntab
    stop
  end if
 
  do i = 2, ntab
    if ( xtab(i) <= xtab(i-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PLINT - Fatal error!'
      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
      write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
      write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
      stop
    end if
  end do
!
!  If A > B, temporarily switch A and B, and store sign.
!
  if ( a > b ) then
    syl = b
    b = a
    a = syl
    ind = -1
  else
    syl = a
    ind = 1
  end if
!
!  Find ILO and IHI so that A <= XTAB(ILO) <= XTAB(IHI) <= B
!  with the possible exception that A and B may be in the same
!  interval, or completely to the right or left of the XTAB's.
!
  ilo = 1
  ihi = ntab
  do i = 1, ntab
    if ( a <= xtab(i) ) then
      exit
    end if
    ilo = ilo+1
  end do
  
  do i = 1, ntab
    if ( b >= xtab(i) ) then
      exit
    end if
    ihi = ihi-1
  end do
!
!  Treat special cases where A, B lie both to left or both to right
!  of XTAB interval, or inbetween same pair of XTAB's.
!
  if ( ihi == 0 ) then
    slope = (ftab(2)-ftab(1))/(xtab(2)-xtab(1))
    fa = ftab(1) + slope*(a-xtab(1))
    fb = ftab(1) + slope*(b-xtab(1))
    result = 0.5 * (b-a) * (fa+fb)
    go to 110
  else if ( ilo == ntab+1 ) then
    slope = (ftab(ntab)-ftab(ntab-1))/(xtab(ntab)-xtab(ntab-1))
    fa = ftab(ntab-1)+slope*(a-xtab(ntab-1))
    fb = ftab(ntab-1)+slope*(b-xtab(ntab-1))
    result = 0.5 * (b-a) * (fa+fb)
    go to 110
  else if ( ihi+1 == ilo ) then
    slope = (ftab(ilo)-ftab(ihi))/(xtab(ilo)-xtab(ihi))
    fa = ftab(ihi)+slope*(a-xtab(ihi))
    fb = ftab(ihi)+slope*(b-xtab(ihi))
    result = 0.5 * (b-a) * (fa+fb)
    go to 110
  end if
!
!  Carry out approximate integration.  We know that ILO is no greater
!  than IHI-1, but equality is possible; A and B may be on either side
!  of a single XTAB(I).  That's OK, then the loop below won't be executed
!  at all.
!
  result = 0.0D+00
  do i = ilo, ihi-1
    result = result + 0.5 * (xtab(i+1)-xtab(i))*(ftab(i)+ftab(i+1))
  end do
!
!  Add contribution from A-ILO and IHI-B.
!  Still have to watch out if ILO = 1 or IHI=NTAB...
!
  if ( ilo == 1 ) then
    slope = (ftab(2)-ftab(1)) / (xtab(2)-xtab(1))
    fa = ftab(1) + slope*(a-xtab(1))
    result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
  else
    slope = (ftab(ilo)-ftab(ilo-1)) / (xtab(ilo)-xtab(ilo-1))
    fa = ftab(ilo-1) + slope*(a-xtab(ilo-1))
    result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
  end if
 
  if ( ihi == ntab ) then
    slope = (ftab(ntab)-ftab(ntab-1)) / (xtab(ntab)-xtab(ntab-1))
    fb = ftab(ntab-1) + slope*(b-xtab(ntab-1))
    result = result + 0.5*(b-xtab(ntab))*(fb+ftab(ntab))
  else
    slope = (ftab(ihi+1)-ftab(ihi)) / (xtab(ihi+1)-xtab(ihi))
    fb = ftab(ihi) + slope*(b-xtab(ihi))
    result = result + 0.5*(b-xtab(ihi))*(fb+ftab(ihi))
  end if
!
!  Restore original values of A and B, reverse sign of integral
!  because of earlier switch.
!
110   continue
 
  if ( ind /= 1 ) then
    ind = 1
    syl = b
    b = a
    a = syl
    result = -result
  end if
 
  return
end
subroutine qnc79 ( func, a, b, err, result, ierr, k )
!
!***********************************************************************
!
!! QNC79 approximates the integral of F(X) using Newton-Cotes quadrature.
!
!
!  Discussion:
!
!    QNC79 is a general purpose program for evaluation of one
!    dimensional integrals  of user defined functions.  QNC79 will
!    pick its own points for evaluation of the integrand and these
!    will vary from problem to problem.
!
!    Thus QNC79 is not designed to integrate over data sets.
!
!    Moderately smooth integrands will be integrated efficiently
!    and reliably.  For problems with strong singularities,
!    oscillations etc., the user may wish to use more sophisticated
!    routines such as those in QUADPACK.
!
!    One measure of the reliability of QNC79 is the output parameter
!    K, giving the number of integrand evaluations that were needed.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!      FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real*8 A, lower limit of integral.
!
!    Input, real*8 B, upper limit of integral.
!
!    Input, real*8 ERR, is a requested error tolerance.
!    Normally pick a value, 0 .LT. ERR .LT. 1.E-3.
!
!    Output, real*8 RESULT, computed value of the integral.
!    Hopefully RESULT is accurate to within ERR times the
!    integral of ABS ( FUNC(X) ).
!
!    Output, integer IERR, a status code
!     1 RESULT most likely meets requested error tolerance.
!    -1 A and B are too nearly equal to allow normal integration.
!     2 RESULT probably does not meet requested error tolerance.
!
!    Output, integer K, the number of function evaluations
!    actually used to do the integration.
!    A value of K .GT. 1000 indicates a difficult problem.
!    Other programs may be more efficient.
!    QNC79 will gracefully give up if K exceeds 2000.
!
  implicit real*8 (a-h,o-z)
!
  integer, parameter :: kmx = 2000
!
  real*8 a
  real*8 aa(40)
  real*8 ae
  real*8 area
  real*8 b
  real*8 bank
  real*8 blocal
  real*8 c
  real*8 ce
  real*8 ee
  real*8 ef
  real*8 eps
  real*8 err
  real*8 f(13)
  real*8 f1(40)
  real*8 f2(40)
  real*8 f3(40)
  real*8 f4(40)
  real*8 f5(40)
  real*8 f6(40)
  real*8 f7(40)
  real*8, external :: func
  real*8 hh(40)
  integer i
  integer i1mach
  integer, save :: icall = 0
  integer ierr
  integer k
  integer, save :: kml = 7
  integer l
  integer lmn
  integer lmx
  integer lr(40)
  integer nbits
  integer nib
  integer, save :: nlmn = 2
  integer nlmx
  real*8 q13
  real*8 q7
  real*8 q7l
  real*8 q7r(40)
  real*8 d1mach
  real*8 result
  real*8 test
  real*8 tol
  real*8 vl(40)
  real*8 vr
  real*8 w1
  real*8 w2
  real*8 w3
  real*8 w4
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( icall /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Fatal error!'
    write ( *, '(a)' ) '  QNC79 was called recursively!'
    stop
  end if
 
  icall = 1
  w1 = 41.0D+00 / 140.0D+00
  w2  = 216.0 / 140.0D+00
  w3 = 27.0 / 140.0D+00
  w4  = 272.0 / 140.0D+00
  nbits = int ( d1mach(5) * real(i1mach(11)) / 0.30102000 )
  nlmx = min ( 40, (nbits*4)/5 )
  result = 0.0D+00
  ierr = 1
  ce = 0.0D+00
  lmx = nlmx
  lmn = nlmn
  if ( b == 0.0D+00 ) go to 3
  if ( sign ( 1.0D+00, b ) * a <= 0.0D+00 ) go to 3
  c = abs ( 1.0D+00 - a / b )

  if ( c > 0.1d+00 ) then
    go to 3
  end if
 
  if ( c <= 0.0D+00 ) then
    ierr  = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Fatal error!'
    write ( *, '(a)' ) '  A and B are too close.'
    stop
  end if
 
  nib = int ( 0.5D+00 - log(c) / log(2.0D+00) )
  lmx = min ( nlmx, nbits-nib-4)

  if ( lmx < 2 ) then
    go to 32
  end if

  lmn = min(lmn,lmx)
 
3 continue
 
  tol = max ( abs ( err ), 2.0D+00**(5-nbits) )
  if ( err == 0.0D+00 ) tol = sqrt ( epsilon ( tol ) )
  eps = tol
  hh(1) = (b-a) / 12.0D+00
  aa(1) = a
  lr(1) = 1
 
  do i = 1, 11, 2
    f(i) = func(a+real(i-1)*hh(1))
  end do
 
  blocal = b
  f(13) = func(blocal)
  k = 7
  l = 1
  area = 0.0D+00
  q7 = 0.0D+00
  ef = 256.0D+00 / 255.0D+00
  bank = 0.0D+00
!
!  Compute refined estimates, estimate the error, etc.
!
5 continue
 
  do i = 2, 12, 2
    f(i) = func ( aa(l)+real(i-1)*hh(l) )
  end do
 
  k = k+6
!
!  Compute left and right half estimates.
!
  q7l = hh(l)*( ( w1*(f(1)+f(7)) + w2*(f(2)+f(6)) ) &
              + ( w3*(f(3)+f(5)) + w4*f(4) ) )

  q7r(l) = hh(l)*( ( w1*(f(7)+f(13)) + w2*(f(8)+f(12)) ) &
                + ( w3*(f(9)+f(11)) + w4*f(10) ) )
!
!  Update estimate of integral of absolute value.
!
  area = area + ( abs ( q7l ) + abs ( q7r(l) ) - abs ( q7 ) )
!
!  Do not bother to test convergence before minimum refinement level.
!
  if (l<lmn) go to 11
!
!  Estimate the error in new value for whole interval, Q13.
!
  q13 = q7l+q7r(l)
  ee = abs ( q7 - q13 ) * ef
!
!  Compute nominal allowed error.
!
  ae = eps*area
!
!  Borrow from bank account, but not too much.
!
  test = min(ae+0.8d+00*bank, 10.0D+00*ae)
!
!  Don't ask for excessive accuracy.
!
  test = max ( test, tol * abs ( q13 ), 0.00003d+00 * tol * area )
!
!  Now, did this interval pass or not?
!
  if ( ee <= test )go to 8
  go to 10
!
!  Have hit max refinement level - penalize the cumulative error.
!
    7 continue
 
  ce = ce+(q7-q13)
  go to 9
!
!  On good intervals accumulate the theoretical estimate.
!
    8 ce = ce+(q7-q13) / 255.0D+00
!
!  Update the bank account.  Don't go into debt.
!
9 continue

  bank = bank + (ae-ee)
  if ( bank < 0.0D+00 ) bank = 0.0D+00
!
!  Did we just finish a left half or a right half?
!
  if ( lr(l) <= 0 ) go to 15
  go to 20
!
!  Consider the left half of the next deeper level.
!
10 continue
 
  if ( k > kmx ) lmx = min(kml,lmx)
  if ( l >= lmx ) go to 7

11 continue

  l = l+1
  eps = eps * 0.5D+00
  if ( l <= 17 ) ef = ef/sqrt(2.0D+00)
  hh(l) = hh(l-1)*0.5D+00
  lr(l) = -1
  aa(l) = aa(l-1)
  q7 = q7l
  f1(l) = f(7)
  f2(l) = f(8)
  f3(l) = f(9)
  f4(l) = f(10)
  f5(l) = f(11)
  f6(l) = f(12)
  f7(l) = f(13)
  f(13) = f(7)
  f(11) = f(6)
  f(9)  = f(5)
  f(7)  = f(4)
  f(5)  = f(3)
  f(3)  = f(2)
  go to 5
!
!  Proceed to right half at this level
!
   15 vl(l) = q13
   16 q7 = q7r(l-1)
  lr(l) = 1
  aa(l) = aa(l)+12.0D+00 * hh(l)
  f(1)  = f1(l)
  f(3)  = f2(l)
  f(5)  = f3(l)
  f(7)  = f4(l)
  f(9)  = f5(l)
  f(11) = f6(l)
  f(13) = f7(l)
  go to 5
!
!  Left and right halves are done, so go back up a level
!
   20 vr = q13
   22 if ( l <= 1 ) go to 30

  if ( l <= 17 ) then
    ef = ef * sqrt ( 2.0D+00 )
  end if

  eps = eps * 2.0D+00
  l = l-1
 
  if ( lr(l) <= 0 ) then
    vl(l) = vl(l+1)+vr
    go to 16
  else
    vr = vl(l+1)+vr
    go to 22
  end if
 
   30 continue
 
  if ( abs ( ce ) > 2.0D+00 * tol * area ) then
    ierr = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Warning!'
    write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
  end if
 
32    continue
 
  result = vr
 
  icall = 0
 
  return
end
subroutine quad ( func, a, b, abserr, relerr, nleast, nmost, work, &
  result )
!
!***********************************************************************
!
!! QUAD approximates the integral of F(X) by Romberg integration.
!
!
!  Discussion:
!
!    The integration is repeated until convergence of the results, 
!    or the maximum number of steps is taken.  The Romberg method 
!    uses the trapezoidal rule, subinterval halving, and Richardson 
!    extrapolation.
!
!    Convergence is declared if either of the following occurs:
!
!      ABS ( RESULT - INTEGRAL ) < ABSERR
!
!    or
!
!      RESULT = integral from A to B (1+Y(X))*FUNC(X)DX for some
!      function Y with ABS ( Y(X) ) < RELERR+RNDERR  with RNDERR the
!      machine rounding error.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    C F Dunkl,
!    Romberg quadrature to prescribed accuracy,
!    SHARE file number 7090-1481 TYQUAD
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC, the name of the function to be integrated.
!    The user must declare the name an external parameter in the
!    calling program, pass the name of the function in FUNC,
!    and write a function of the form FUNCTION FUNC(X) which
!    evaluates the function at the point X.
!
!    Input, real*8 A, the lower limit of integration.
!
!    Input, real*8 B, the upper limit of integration.
!
!    Input, real*8 ABSERR, the absolute error tolerance.
!
!    Input, real*8 RELERR, the relative error tolerance.
!
!    Input, integer NLEAST, the least number of times the integration
!    is to be carried out before the convergence test is made.
!    A value 3 <= NLEAST <= 15 is suggested.
!
!    Input, integer NMOST, the most number of times the
!    integration is to be carried out.
!
!    Workspace, real*8 WORK(NMOST+1).
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer nmost
!
  real*8 a
  real*8 abserr
  real*8 b
  real*8 f
  real*8 fcna
  real*8 fcnb
  real*8 fcnxi
  real*8, external :: func
  real*8 h
  integer i
  integer j
  integer k
  integer nleast
  integer nx
  real*8 qx1
  real*8 qx2
  real*8 relerr
  real*8 result
  real*8 rnderr
  real*8 sum1
  real*8 sumabs
  real*8 t
  real*8 tabs
  real*8 work(nmost+1)
  real*8 x
!
  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  rnderr = epsilon ( 1.0D+00 )
  qx1 = 0
  qx2 = 0
  h = b-a
  fcna = func(a)
  fcnb = func(b)
  tabs = abs ( h ) * ( abs ( fcna ) + abs ( fcnb ) ) / 2.0D+00
  t = h*(fcna+fcnb) / 2.0D+00
  nx = 1
 
  do i = 1, nmost
 
    h = 0.5*h
 
    sum1 = 0.0D+00
    sumabs = 0.0D+00
    do j = 1, nx
      fcnxi = func(a+h*(2*j-1))
      sumabs = sumabs + abs ( fcnxi )
      sum1 = sum1+fcnxi
    end do
 
    t = 0.5D+00 * t + h * sum1
    tabs = tabs/2.0D+00 + abs ( h ) * sumabs
    work(i) = 2.0D+00*(t+h*sum1)/3.0D+00

    if ( i > 1 ) then
!
!  Construct difference table for Richardson extrapolation.
!
      f = 4.0D+00
      do j = 2, i
        k = i+1-j
        f = f*4.0D+00
        work(k) = work(k+1)+(work(k+1)-work(k) ) / ( f - 1.0D+00 )
      end do
!
!  Perform acceptance check.
!
      if ( i >= nleast ) then
        x = abs ( work(1)-qx2 ) + abs ( qx2 - qx1 )
        if ( x <= 3.0D+00 * tabs * ( abs ( relerr ) + rnderr ) ) go to 10
        if ( x <= 3.0D+00 * abs ( abserr ) ) go to 10
      end if
!
!  Save old result, perform bisection, repeat.
!
      qx1 = qx2
    end if
 
    qx2 = work(1)
    nx = nx*2
 
  end do
!
!  Accept result.
!
10    continue
  result = work(1)
  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP swaps two real*8 values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real*8 X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit real*8 (a-h,o-z)
!
  real*8 x
  real*8 y
  real*8 z
!
  z = x
  x = y
  y = z

  return
end
subroutine rminsp ( func, a, b, epsin, epsout, iop, result )
!
!***********************************************************************
!
!! RMINSP approximates the integral of a function using Romberg integration.
!
!
!  Discussion:
!
!    Both the midpoint and trapezoidal rule are used,
!    the intervals are repeatedly bisected, and Richardson
!    extrapolation is carried out to achieve a high accuracy.
!
!    RMINSP can carry out a cosine-transform of the integral.  The
!    only effect this has is to handle a function F(X) which has
!    singularities near the endpoints.  This transform is based on
!    the fact that
!
!      Integral from -1 to 1  ( F(    X))       DX
!
!    equals
!
!      Integral from  0 to PI ( F(COS(Z))*SIN(Z) )  DZ
!
!    If suitable accuracy is not achieved, the internal variable
!    NUPPER might be increased.  Its current value of 9 corresponds
!    to a maximum of 1024 subintervals and 1025 function evaluations.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    T Havie,
!    Algorithm 257,
!    Communications of the Association for Computing Machinery,
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!       FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real*8 A, lower limit of integration.
!
!    Input, real*8 B, upper limit of integration.
!
!    Input, real*8 EPSIN, requested relative error tolerance.
!
!    Output, real*8 EPSOUT, estimated achieved relative error.
!
!    Input, integer IOP, method switch:
!    1, Use ordinary algorithm.
!    2, Use cosine transformation to decrease effect of
!       singularities near the endpoints.
!
!    Output, real*8 RESULT, the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer, parameter :: nupper = 9
!
  real*8 a
  real*8 acof(11)
  real*8 alf
  real*8 alfnj
  real*8 alfno
  real*8 ar
  real*8 b
  real*8 bcof(nupper+1)
  real*8 bet
  real*8 betnj
  real*8 betno
  real*8 const1
  real*8 const2
  real*8 deltan
  real*8 endpts
  real*8 epsin
  real*8 epsout
  real*8 error
  real*8, parameter :: fac1 = 0.411233516712057d+00
  real*8, parameter :: fac2 = 0.822467033441132d+00
  real*8 factor
  real*8, external :: func
  real*8 gamman
  real*8 hnstep
  integer i
  integer index
  integer iop
  integer iout
  integer j
  integer n
  integer nhalf
  real*8 pi
  real*8 r1
  real*8 r2
  real*8 rn
  real*8 rnderr
  real*8 result
  real*8 rounde
  real*8 tend
  real*8 triarg
  real*8 umid
  real*8 xmin
  real*8 xplus
!
  result = 0.0D+00
 
  if ( a == b ) then
    return
  end if
!
!  Set coefficients in formula for accumulated roundoff error,
!  rounde = rnderr*(r1+r2*n), where r1, r2 are two empirical constants
!  and n is the current number of function values used.
!
  rnderr = epsilon ( 1.0D+00 )
 
  r1 = 1.0D+00
  r2 = 2.0D+00
  if ( iop==2 ) r1 = 50.0D+00
  if ( iop==1 ) r2 = 0.01d+00 * r2
  error = epsin
!
!  Initial calculations.
!
  alf = 0.5D+00 * (b-a)
  bet = 0.5D+00 * (b+a)
  acof(1) = func(a)+func(b)
  bcof(1) = func(bet)
!
!  Modified Romberg algorithm, ordinary case.
!
  if ( iop /= 2 ) then

    hnstep = 2.0D+00
    bcof(1) = hnstep*bcof(1)
    factor = 1.0D+00
!
!  Modified Romberg, cosine transformed case.
!
  else
    hnstep = pi()
    ar = fac1
    endpts = acof(1)
    acof(1) = fac2*acof(1)
    bcof(1) = hnstep*bcof(1)-ar*endpts
    factor = 4.0D+00
    ar = ar/4.0D+00
    triarg = pi() / 4.0D+00
    alfno = -1.0D+00
  end if
 
  hnstep = 0.5D+00 * hnstep
  nhalf = 1
  n = 2
  rn = 2.0D+00
  acof(1) = 0.5D+00 * (acof(1)+bcof(1))
  acof(2) = acof(1)-(acof(1)-bcof(1))/(4.0D+00*factor-1.0D+00)
!
!  End of initial calculation.
!
!  Start actual calculations.
!
  do i = 1, nupper
 
    umid = 0.0D+00
!
!  Modified Romberg algorithm, ordinary case.
!  compute first element in mid-point formula for ordinary case
!
    if ( iop == 1 ) then
 
      alfnj = 0.5D+00*hnstep
 
      do j = 1, nhalf
        xplus = alf*alfnj+bet
        xmin = -alf*alfnj+bet
        umid = umid + func(xplus)+func(xmin)
        alfnj = alfnj+hnstep
      end do
 
      umid = hnstep*umid
!
!  Modified Romberg algorithm, cosine transformed case
!  compute first element in mid-point formula for cosine transformed
!  Romberg algorithm
!
    else if ( iop == 2 ) then
 
      const1 = -sin(triarg)
      const2 = 0.5D+00 * alfno/const1
 
      alfno = const1
      betno = const2
      gamman = 1.0D+00 - 2.0D+00 * alfno**2
      deltan = -2.0D+00 * alfno*betno
 
      do j = 1, nhalf
        alfnj = gamman*const1+deltan*const2
        betnj = gamman*const2-deltan*const1
        xplus = alf*alfnj+bet
        xmin = -alf*alfnj+bet
        umid = umid+betnj*(func(xplus)+func(xmin))
        const1 = alfnj
        const2 = betnj
      end do
 
      umid = hnstep*umid-ar*endpts
      ar = ar / 4.0D+00
 
    end if
!
!  Modified Romberg algorithm, calculate (i+1)-th row in the U table
!
    const1 = 4.0D+00 * factor
    index = i+1
 
    do j = 2, i+1
      tend = umid + ( umid - bcof(j-1) ) / ( const1 - 1.0D+00 )
      bcof(j-1) = umid
      umid = tend
      const1 = 4.0D+00 * const1
    end do
 
    bcof(i+1) = tend
    xplus = const1
!
!  Calculation of (i+1)-th row in the U table is finished
!
!  Test to see if the required accuracy has been obtained.
!
    epsout = 1.0D+00
    iout = 1
 
    do j = 1, index
 
      const1 = 0.5D+00 * ( acof(j) + bcof(j) )
      const2 = 0.5D+00 * abs ( ( acof(j) - bcof(j) ) / const1 )
 
      if ( const2 <= epsout ) then
        epsout = const2
        iout = j
      end if
 
      acof(j) = const1
 
    end do
!
!  Testing on accuracy finished
!
    if (iout==index) iout=iout+1
    acof(index+1) = acof(index)-(acof(index)-bcof(index))/(xplus-1.0)
    rounde = rnderr*(r1+r2*rn)

    epsout = max ( epsout, rounde )
    error = max ( error, rounde )
 
    if ( epsout <= error ) go to 10
 
    nhalf = n
    n = 2 * n
    rn = 2.0D+00 * rn
    hnstep = 0.5D+00 * hnstep

    if ( iop > 1 ) then
      triarg = 0.5 * triarg
    end if
 
  end do
!
!  Accuracy not reached with maximum number of subdivisions
!
  n = nhalf
!
!  Calculation for modified Romberg algorithm finished
!
  10  continue
 
  n = 2*n
  index = index-1
  n = n+1
  j = iout
  if ((j-1)>=index) j = index
  tend = alf * (2.0D+00*acof(j)-bcof(j))
  umid = alf * bcof(j)
  result = alf * acof(iout)

  return
end
subroutine rvec_even ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_EVEN returns N real*8 values, evenly spaced between ALO and AHI.
!
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real*8 A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
!
  implicit real*8 (a-h,o-z)
!
  integer n
!
  real*8 a(n)
  real*8 ahi
  real*8 alo
  integer i
!
  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( dble ( n - i ) * alo + dble ( i - 1 ) * ahi ) / dble ( n - 1 )
    end do

  end if

  return
end
subroutine simp ( func, a, b, eps, result )
!
!***********************************************************************
!
!! SIMP approximates the integral of a function using an adaptive Simpson's rule.
!
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    J N Lyness,
!    Algorithm 379,
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the Association for Computing Machinery,
!    Volume 13 (1970), pages 260-263.
!
!    W M McKeeman and L Tesler,
!    Algorithm 182,
!    Nonrecursive adaptive integration,
!    Communications of the Association for Computing Machinery,
!    Volume 6 (1963), page 315.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!    FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real*8 A, the lower limit of integration.
!
!    Input, real*8 B, the upper limit of integration.
!
!    Input, real*8 EPS, the requested error tolerance.
!
!    Output, real*8 RESULT, the approximation to the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer, parameter :: maxlev = 30
!
  real*8 a
  real*8 a1
  real*8 absar
  real*8 b
  real*8 da
  real*8 dx(maxlev)
  real*8 ep
  real*8 eps
  real*8 epsp(maxlev)
  real*8 est
  real*8 est1
  real*8 est2(maxlev)
  real*8 est3(maxlev)
  real*8 f1
  real*8 f2(maxlev)
  real*8 f3(maxlev)
  real*8 f4(maxlev)
  real*8 fa
  real*8 fb
  real*8 fbp(maxlev)
  real*8 fm
  real*8 fmp(maxlev)
  real*8, external :: func
  integer i
  integer j
  integer l
  integer lvl
  integer nrtr(maxlev)
  real*8 pval(maxlev,3)
  real*8 result
  real*8 sum1
  real*8 sx
  real*8 x2(maxlev)
  real*8 x3(maxlev)
!
  result = 0.0D+00
 
  if ( a == b ) then
    return
  end if
 
  ep = eps
  a1 = a
  nrtr(1:maxlev) = 0
  pval(1:maxlev,1:3) = 0.0D+00
 
  lvl = 0
  absar = 0.0D+00
  est = 0.0D+00
  da = b - a1

  fa = func ( a1 )
  fm = 4.0D+00 * func ( (a1+b) * 0.5D+00 )
  fb = func ( b )
!
!  1 = RECUR
!
   30 continue
 
  lvl = lvl + 1
  dx(lvl) = da / 3.0D+00
  sx = dx(lvl)/6.0D+00
  f1 = 4.0D+00 * func(0.5*dx(lvl)+a1)
  x2(lvl) = a1+dx(lvl)
  f2(lvl) = func(x2(lvl))
  x3(lvl) = x2(lvl)+dx(lvl)
  f3(lvl) = func(x3(lvl))
  epsp(lvl) = ep
  f4(lvl) = 4.0D+00 * func(dx(lvl)*0.5D+00+x3(lvl))
  fmp(lvl) = fm
  est1 = sx*(fa+f1+f2(lvl))
  fbp(lvl) = fb
  est2(lvl) = sx * (f2(lvl)+f3(lvl)+fm)
  est3(lvl) = sx * (f3(lvl)+f4(lvl)+fb)
  sum1 = est1+est2(lvl)+est3(lvl)
  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
    + abs ( est3(lvl) )
  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) go to 40
  if ( lvl >= maxlev ) go to 50
!
!  2 = UP
!
40 continue
 
  if ( lvl > 1 ) then
    lvl = lvl-1
  end if

  l = nrtr(lvl)

  if ( l == 0 ) then
    go to 50
  end if

  pval(lvl,l) = sum1

  if ( l == 1 ) go to 60
  if ( l == 2 ) go to 70
  if ( l == 3 ) go to 80
 
50 continue
 
  nrtr(lvl) = 1
  est = est1
  fm = f1
  fb = f2(lvl)
  ep = epsp(lvl) / 1.7d+00
  da = dx(lvl)
  go to 30
 
60 continue
 
  nrtr(lvl) = 2
  fa = f2(lvl)
  fm = fmp(lvl)
  fb = f3(lvl)
  est = est2(lvl)
  a1 = x2(lvl)
  ep = epsp(lvl) / 1.7d+00
  da = dx(lvl)
  go to 30
 
70 continue
 
  nrtr(lvl) = 3
  fa = f3(lvl)
  fm = f4(lvl)
  fb = fbp(lvl)
  est = est3(lvl)
  a1 = x3(lvl)
  ep = epsp(lvl) / 1.7d+00
  da = dx(lvl)
  go to 30
 
80 continue
 
  sum1 = pval(lvl,1)+pval(lvl,2)+pval(lvl,3)

  if ( lvl > 1 ) then
    go to 40
  end if
 
90 continue
 
  result = sum1
 
  return
end
subroutine simpne ( x, y, num, result )
!
!***********************************************************************
!
!! SIMPNE approximates the integral of unevenly spaced data.
!
!
!  Discussion:
!
!    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
!    to the data and integrates that exactly.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8 X(NUM), contains the X values of the data, in order.
!
!    Input, real*8 Y(NUM), contains the Y values of the data.
!
!    Input, integer NUM, number of data points.  NUM must be at least 3.
!
!    Output, real*8 RESULT.
!    RESULT is the approximate value of the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer num
!
  real*8 del(3)
  real*8 e
  real*8 f
  real*8 feints
  real*8 g(3)
  integer i
  integer n
  real*8 pi(3)
  real*8 result
  real*8 sum1
  real*8 x(num)
  real*8 x1
  real*8 x2
  real*8 x3
  real*8 y(num)
!
  result = 0.0D+00
 
  if ( num <= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIMPNE - Fatal error!'
    write ( *, '(a)' ) '  NUM <= 2.'
    stop
  end if
 
  n = 1
 
  do
 
    x1 = x(n)
    x2 = x(n+1)
    x3 = x(n+2)
    e = x3*x3-x1*x1
    f = x3*x3*x3-x1*x1*x1
    feints = x3-x1
    del(1) = x3-x2
    del(2) = x1-x3
    del(3) = x2-x1
    g(1) = x2+x3
    g(2) = x1+x3
    g(3) = x1+x2
    pi(1) = x2*x3
    pi(2) = x1*x3
    pi(3) = x1*x2
 
    sum1 = 0.0D+00
    do i = 1, 3
      sum1 = sum1 + y(n-1+i)*del(i)*(f/3.0D+00-g(i)*0.5D+00*e+pi(i)*feints)
    end do
    result = result - sum1 / ( del(1) * del(2) * del(3) )
 
    n = n+2

    if ( n + 1 >= num ) then
      exit
    end if

  end do
 
  if ( mod(num,2) /= 0 ) then
    return
  end if

  n = num-2
  x3 = x(num)
  x2 = x(num-1)
  x1 = x(num-2)
  e = x3*x3-x2*x2
  f = x3*x3*x3-x2*x2*x2
  feints = x3-x2
  del(1) = x3-x2
  del(2) = x1-x3
  del(3) = x2-x1
  g(1) = x2+x3
  g(2) = x1+x3
  g(3) = x1+x2
  pi(1) = x2*x3
  pi(2) = x1*x3
  pi(3) = x1*x2
 
  sum1 = 0.0D+00
  do i = 1, 3
    sum1 = sum1 + y(n-1+i) * del(i) * &
      ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
  end do
 
  result = result - sum1 / ( del(1) * del(2) * del(3) )
 
  return
end
subroutine simpsn ( h, y, num, result )
!
!***********************************************************************
!
!! SIMPSN approximates the integral of evenly spaced data.
!
!
!  Discussion:
!
!    Simpson's rule is used.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real*8 H, specifies the increment between the
!    X values.  Note that the actual X values are not needed,
!    just the constant spacing!
!
!    Input, real*8 Y(NUM), the data.
!
!    Input, integer NUM, the number of data points.  NUM must be at least 3.
!
!    Output, real*8 RESULT, the value of the integral
!    from the first to the last point.
!
  implicit real*8 (a-h,o-z)
!
  integer num
!
  real*8 del(3)
  real*8 f
  real*8 g(3)
  real*8 h
  integer i
  integer n
  real*8 pii(3)
  real*8 result
  real*8 sum1
  real*8 y(num)
!
  result = 0.0D+00
 
  if ( num <= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIMPSN - Fatal error!'
    write ( *, '(a,i6)' ) '  NUM < 2, NUM = ', num
    stop
  end if
 
  if ( mod ( num, 2 ) == 0 ) then
    n = num-1
  else
    n = num
  end if
 
  result = y(1) + y(n) + 4.0D+00 * y(n-1)
  do i = 2, n-2, 2
    result = result + 4.0D+00 * y(i) + 2.0D+00 * y(i+1)
  end do
  result = h * result / 3.0D+00
 
  if ( mod(num,2) == 1 ) then
    return
  end if
 
  f = h*h*h
  del(1) = h
  del(2) = -2.0D+00 * h
  del(3) = h
  g(1) = h
  g(2) = 0.0D+00
  g(3) = -h
  pii(1) = 0.0D+00
  pii(2) = -h*h
  pii(3) = 0.0D+00
  n = n-1
 
  sum1 = 0.0D+00
  do i = 1, 3
    sum1 = sum1 + y(n-1+i) * del(i) * &
      ( f / 3.0D+00 - g(i) * 0.5D+00 * h * h + pii(i) * h )
  end do
 
  result = result + 0.5D+00 * sum1 / h**3
 
  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit real*8 (a-h,o-z)
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine wedint ( ftab, h, ntab, result )
!
!***********************************************************************
!
!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!
!
!  Modified:
!
!    30 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), contains the tabulated data values.
!
!    Input, real*8 H, is the spacing between the points at which the data
!    was evaluated.
!
!    Input, integer NTAB, is the number of data points.  (NTAB-1) must be
!    divisible by 6.
!
!    Output, real*8 RESULT, is the approximation to the integral.
!
  implicit real*8 (a-h,o-z)
!
  integer ntab
!
  real*8 ftab(ntab)
  real*8 h
  integer i
  real*8 result
!
  result = 0.0D+00
 
  if ( ntab <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB < 2'
    write ( *, '(a,i6)' ) '  NTAB = ', ntab
    stop
  end if
 
  if ( mod ( ntab, 6 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB must equal 6*N+1 for some N!'
    stop
  end if
 
  do i = 1, ntab-6, 6
    result = result & 
      +           ftab(i) &
      + 5.0D+00 * ftab(i+1) &
      +           ftab(i+2) &
      + 6.0D+00 * ftab(i+3) &
      +           ftab(i+4) &
      + 5.0D+00 * ftab(i+5) &
      +           ftab(i+6)
  end do
 
  result = 3.0D+00 * h * result / 10.0D+00
 
  return
end
