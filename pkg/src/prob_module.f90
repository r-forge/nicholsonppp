include 'mers_twist.f90'

module prob_mod

  use mtmod ! mersenne twister stuff et changement de seed
!!!!!!!!!!!!!!!!!!!!!!!!
!!!module contenant differents codes pour generer des nombre aleatoires
!!!inspire de code precedents
!!!oin utilise le mersenne twister pour generer les nombre aleatoire d'une loi nuniforme
!!!puis on s'inspire de random.f90 pour prndre les lois normales...

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!FONCTION random_normal (avec appel au mersenne twister
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION random_normal()!(u,v_in) 
    
    ! Adapted from the following Fortran 77 code
    !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
    !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
    
    !  The function random_normal() returns a normally distributed pseudo-random
    !  number with zero mean and unit variance.
    
    !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
    !  and J.F. Monahan augmented with quadratic bounding curves.

    !REAL, intent(in) :: u,v_in !genere avec un mersenne twiste
    REAL :: random_normal
    !     Local variables
    REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
    r1 = 0.27597, r2 = 0.27846, half=0.5, x, y, q ,u, v
    !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
    DO
       u=grnd()
       v=grnd()
       v = 1.7156 * (v - half)
       !     Evaluate the quadratic form
       x = u - s
       y = ABS(v) - t
       q = x**2 + y*(a*y - b*x)
       
       !     Accept P if inside inner ellipse
       IF (q < r1) EXIT
       !     Reject P if outside outer ellipse
       IF (q > r2) CYCLE
       !     Reject P if outside acceptance region
       IF (v**2 < -4.0*LOG(u)*u**2) EXIT
    END DO
    !     Return ratio of P's coordinates as the normal deviate
    random_normal = v/u
    RETURN
    
  END FUNCTION random_normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: gaussienne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function norm_pdf(x,mean,var)
    implicit none
    real, intent(in) :: x,mean,var
    real, parameter  :: pi = 3.141592653589793
    real :: norm_pdf

    norm_pdf = exp ( -0.5 * (x-mean) * (x-mean) / var ) / sqrt ( 2.0 * pi * var )

    ! return 
  end function norm_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: log de la gaussienne simplifie (attention a utiliser dans unr apport pour justifier la siplification)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function log_norm_pdf(x,mean,var)
    implicit none
    real, intent(in) :: x,mean,var
    ! real, parameter  :: pi = 3.141592653589793
    real :: log_norm_pdf

    log_norm_pdf = -0.5*log(var) - 0.5 * (x-mean) * (x-mean) / var   ! -0.5*log( 2.0 * pi)

    ! return 
  end function log_norm_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Coefficient C(n,k) utilse pour la binomiale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function binomial_coef (n, k)
    implicit none
    integer, intent(in) :: n,k
    integer :: i,mn,mx
    real :: binomial_coef !sinon probleme d'overflow!!

    mn = min ( k, n-k )
    if ( mn < 0 ) then
       binomial_coef = 0.
    else if ( mn == 0 ) then
       binomial_coef = 1.
    else
       mx = max ( k, n-k )
       binomial_coef = mx + 1.
       do i = 2, mn
          binomial_coef = ( binomial_coef * ( mx + i ) ) / i
       end do
    end if
  end function binomial_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Coefficient GAMMA utilse pour la loi beta...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function gamma_log ( x )

    !*****************************************************************************80
    !
    !! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ).
    !
    !  Discussion:
    !
    !    Computation is based on an algorithm outlined in references 1 and 2.  
    !    The program uses rational functions that theoretically approximate 
    !    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The 
    !    approximation for 12 < X is from Hart et al, while approximations 
    !    for X < 12.0D+00 are similar to those in Cody and Hillstrom, 
    !    but are unpublished.
    !
    !    The accuracy achieved depends on the arithmetic system, the compiler,
    !    intrinsic functions, and proper selection of the machine dependent 
    !    constants.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    16 June 1999
    !
    !  Author: 
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    William Cody, Kenneth Hillstrom, 
    !    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
    !    Mathematics of Computation, 
    !    Volume 21, 1967, pages 198-203.
    !
    !    Kenneth Hillstrom, 
    !    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, 
    !    May 1969.
    ! 
    !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
    !    Charles Mesztenyi, John Rice, Henry Thacher, Christoph Witzgall,
    !    Computer Approximations,
    !    Wiley, 1968.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the Gamma function.  
    !    X must be positive.
    !
    !    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma 
    !    function of X.
    !
    !  Local Parameters:
    !
    !    Local, real ( kind = 8 ) BETA, the radix for the floating-point
    !    representation.
    !
    !    Local, integer MAXEXP, the smallest positive power of BETA that overflows.
    !
    !    Local, real ( kind = 8 ) XBIG, the largest argument for which
    !    LN(GAMMA(X)) is representable in the machine, the solution to the equation
    !      LN(GAMMA(XBIG)) = BETA**MAXEXP.
    !
    !    Local, real ( kind = 8 ) FRTBIG, a rough estimate of the fourth root 
    !    of XBIG.
    !
    !  Approximate values for some important machines are:
    !
    !                            BETA      MAXEXP         XBIG     FRTBIG
    !
    !  CRAY-1        (S.P.)        2        8191       9.62D+2461  3.13D+615
    !  Cyber 180/855 (S.P.)        2        1070       1.72D+319   6.44D+79
    !  IEEE (IBM/XT) (S.P.)        2         128       4.08D+36    1.42D+9
    !  IEEE (IBM/XT) (D.P.)        2        1024       2.55D+305   2.25D+76
    !  IBM 3033      (D.P.)       16          63       4.29D+73    2.56D+18
    !  VAX D-Format  (D.P.)        2         127       2.05D+36    1.20D+9
    !  VAX G-Format  (D.P.)        2        1023       1.28D+305   1.89D+76
    !                    
    implicit none

    real    ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
         -1.910444077728D-03, &
         8.4171387781295D-04, &
         -5.952379913043012D-04, &
         7.93650793500350248D-04, &
         -2.777777777777681622553D-03, &
         8.333333333333333331554247D-02, &
         5.7083835261D-03 /)
    real    ( kind = 8 ) corr
    real    ( kind = 8 ), parameter :: d1 = -5.772156649015328605195174D-01
    real    ( kind = 8 ), parameter :: d2 =  4.227843350984671393993777D-01
    real    ( kind = 8 ), parameter :: d4 =  1.791759469228055000094023D+00
    integer ( kind = 4 ) i
    real    ( kind = 8 ), parameter :: frtbig = 1.42D+09
    real    ( kind = 8 ) gamma_log
    real    ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
         4.945235359296727046734888D+00, &
         2.018112620856775083915565D+02, &
         2.290838373831346393026739D+03, &
         1.131967205903380828685045D+04, &
         2.855724635671635335736389D+04, &
         3.848496228443793359990269D+04, &
         2.637748787624195437963534D+04, &
         7.225813979700288197698961D+03 /)
    real    ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
         4.974607845568932035012064D+00, &
         5.424138599891070494101986D+02, &
         1.550693864978364947665077D+04, &
         1.847932904445632425417223D+05, &
         1.088204769468828767498470D+06, &
         3.338152967987029735917223D+06, &
         5.106661678927352456275255D+06, &
         3.074109054850539556250927D+06 /)
    real    ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
         1.474502166059939948905062D+04, &
         2.426813369486704502836312D+06, &
         1.214755574045093227939592D+08, &
         2.663432449630976949898078D+09, &
         2.940378956634553899906876D+10, &
         1.702665737765398868392998D+11, &
         4.926125793377430887588120D+11, &
         5.606251856223951465078242D+11 /)
    real    ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
    real    ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
         6.748212550303777196073036D+01, &
         1.113332393857199323513008D+03, &
         7.738757056935398733233834D+03, &
         2.763987074403340708898585D+04, &
         5.499310206226157329794414D+04, &
         6.161122180066002127833352D+04, &
         3.635127591501940507276287D+04, &
         8.785536302431013170870835D+03 /)
    real    ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
         1.830328399370592604055942D+02, &
         7.765049321445005871323047D+03, &
         1.331903827966074194402448D+05, &
         1.136705821321969608938755D+06, &
         5.267964117437946917577538D+06, &
         1.346701454311101692290052D+07, &
         1.782736530353274213975932D+07, &
         9.533095591844353613395747D+06 /)
    real    ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
         2.690530175870899333379843D+03, &
         6.393885654300092398984238D+05, &
         4.135599930241388052042842D+07, &
         1.120872109616147941376570D+09, &
         1.488613728678813811542398D+10, &
         1.016803586272438228077304D+11, &
         3.417476345507377132798597D+11, &
         4.463158187419713286462081D+11 /)
    real    ( kind = 8 ) res
    real    ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
    real    ( kind = 8 ) x
    real    ( kind = 8 ), parameter :: xbig = 4.08D+36
    real    ( kind = 8 ) xden
    real    ( kind = 8 ) xm1
    real    ( kind = 8 ) xm2
    real    ( kind = 8 ) xm4
    real    ( kind = 8 ) xnum
    real    ( kind = 8 ) xsq
    !
    !  Return immediately if the argument is out of range.
    !
    if ( x <= 0.0D+00 .or. xbig < x ) then
       gamma_log = huge ( gamma_log )
       return
    end if

    if ( x <= epsilon ( x ) ) then

       res = -log ( x )

    else if ( x <= 1.5D+00 ) then

       if ( x < pnt68 ) then
          corr = - log ( x )
          xm1 = x
       else
          corr = 0.0D+00
          xm1 = ( x - 0.5D+00 ) - 0.5D+00
       end if

       if ( x <= 0.5D+00 .or. pnt68 <= x ) then

          xden = 1.0D+00
          xnum = 0.0D+00

          do i = 1, 8
             xnum = xnum * xm1 + p1(i)
             xden = xden * xm1 + q1(i)
          end do

          res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

       else

          xm2 = ( x - 0.5D+00 ) - 0.5D+00
          xden = 1.0D+00
          xnum = 0.0D+00
          do i = 1, 8
             xnum = xnum * xm2 + p2(i)
             xden = xden * xm2 + q2(i)
          end do

          res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

       end if

    else if ( x <= 4.0D+00 ) then

       xm2 = x - 2.0D+00
       xden = 1.0D+00
       xnum = 0.0D+00
       do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
       end do

       res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

    else if ( x <= 12.0D+00 ) then

       xm4 = x - 4.0D+00
       xden = - 1.0D+00
       xnum = 0.0D+00
       do i = 1, 8
          xnum = xnum * xm4 + p4(i)
          xden = xden * xm4 + q4(i)
       end do

       res = d4 + xm4 * ( xnum / xden )

    else

       res = 0.0D+00

       if ( x <= frtbig ) then

          res = c(7)
          xsq = x * x

          do i = 1, 6
             res = res / xsq + c(i)
          end do

       end if

       res = res / x
       corr = log ( x )
       res = res + sqrtpi - 0.5D+00 * corr
       res = res + x * ( corr - 1.0D+00 )

    end if

    gamma_log = res

    return
  end function gamma_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Coefficient BETA utilse pour la loi beta...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function beta ( a, b )

    !*****************************************************************************80
    !! BETA returns the value of the Beta function.
    !  Discussion:
    !    The Beta function is defined as
    !      BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
    !                = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
    !  Licensing:
    !    This code is distributed under the GNU LGPL license. 
    !  Modified:
    !    10 July 1998
    !  Author:
    !    John Burkardt
    !  Parameters:
    !    Input, real ( kind = 8 ) A, B, the parameters of the function.
    !    0.0D+00 < A,
    !    0.0D+00 < B.
    !    Output, real ( kind = 8 ) BETA, the value of the function.
    !
    implicit none

    real ( kind = 8 ), intent(in) :: a,b
    real :: beta !,gamma_log

    if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'BETA - Fatal error!'
       write ( *, '(a)' ) '  Both A and B must be greater than 0.'
       stop
    end if
    beta = exp ( gamma_log( a ) + gamma_log( b ) - gamma_log( a + b ) )

    return
  end function beta

  subroutine normal_01_cdf ( x, cdf )

    !*****************************************************************************80
    !
    !! NORMAL_01_CDF evaluates the Normal 01 CDF.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 February 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference: 
    !
    !    AG Adams,
    !    Algorithm 39, 
    !    Areas Under the Normal Curve,
    !    Computer Journal, 
    !    Volume 12, pages 197-198, 1969.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the CDF.
    !
    !    Output, real ( kind = 8 ) CDF, the value of the CDF.
    !
    implicit none

    real    ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
    real    ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
    real    ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
    real    ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
    real    ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
    real    ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
    real    ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
    real    ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
    real    ( kind = 8 ), parameter :: b1 = 3.8052D-08
    real    ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
    real    ( kind = 8 ), parameter :: b3 = 3.98064794D-04
    real    ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
    real    ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
    real    ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
    real    ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
    real    ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
    real    ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
    real    ( kind = 8 ), parameter :: b10 = 30.789933034D+00
    real    ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
    real    ( kind = 8 ) cdf
    real    ( kind = 8 ) q
    real    ( kind = 8 ) x
    real    ( kind = 8 ) y
    !
    !  |X| <= 1.28.
    !
    if ( abs ( x ) <= 1.28D+00 ) then

       y = 0.5D+00 * x * x

       q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
            + a6 / ( y + a7 ) ) ) )
       !
       !  1.28 < |X| <= 12.7
       !
    else if ( abs ( x ) <= 12.7D+00 ) then

       y = 0.5D+00 * x * x

       q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
            + b2 / ( abs ( x ) + b3 &
            + b4 / ( abs ( x ) - b5 &
            + b6 / ( abs ( x ) + b7 &
            - b8 / ( abs ( x ) + b9 &
            + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
       !
       !  12.7 < |X|
       !
    else

       q = 0.0D+00

    end if
    !
    !  Take account of negative X.
    !
    if ( x < 0.0D+00 ) then
       cdf = q
    else
       cdf = 1.0D+00 - q
    end if

    return
  end subroutine normal_01_cdf

  subroutine gamma_cdf ( x, a, b, c, cdf )

    !*****************************************************************************80
    !
    !! GAMMA_CDF evaluates the Gamma CDF.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 January 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the PDF.
    !    A <= X
    !
    !    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
    !    0.0D+00 < B, 
    !    0.0D+00 < C.
    !
    !    Output, real ( kind = 8 ) CDF, the value of the CDF.
    !
    implicit none

    real    ( kind = 8 ) a
    real    ( kind = 8 ) b
    real    ( kind = 8 ) c
    real    ( kind = 8 ) cdf
    !real    ( kind = 8 ) gamma_inc
    real    ( kind = 8 ) p2
    real    ( kind = 8 ) x
    real    ( kind = 8 ) x2

    x2 = ( x - a ) / b
    p2 = c

    cdf = gamma_inc ( p2, x2 )

    return
  end subroutine gamma_cdf

  subroutine chi_square_cdf ( x, a, cdf )

    !*****************************************************************************80
    !
    !! CHI_SQUARE_CDF evaluates the Chi squared CDF.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 October 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the value of the random deviate.
    !
    !    Input, real ( kind = 8 ) A, the parameter of the distribution, usually
    !    the number of degrees of freedom.
    !
    !    Output, real ( kind = 8 ) CDF, the value of the CDF.
    !
    implicit none

    real    ( kind = 8 ) a
    real    ( kind = 8 ) a2
    real    ( kind = 8 ) b2
    real    ( kind = 8 ) c2
    real    ( kind = 8 ) cdf
    real    ( kind = 8 ) x
    real    ( kind = 8 ) x2

    x2 = 0.5D+00 * x

    a2 = 0.0D+00
    b2 = 1.0D+00
    c2 = 0.5D+00 * a

    call gamma_cdf ( x2, a2, b2, c2, cdf )

    return
  end subroutine chi_square_cdf

  function gamma_inc ( p, x )

    !*****************************************************************************80
    !
    !! GAMMA_INC computes the incomplete Gamma function.
    !
    !  Discussion:
    !
    !    GAMMA_INC(P,       0) = 0, 
    !    GAMMA_INC(P,Infinity) = 1.
    !
    !    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    01 May 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by B L Shea.
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    BL Shea,
    !    Chi-squared and Incomplete Gamma Integral,
    !    Algorithm AS239,
    !    Applied Statistics,
    !    Volume 37, Number 3, 1988, pages 466-473.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P, the exponent parameter.
    !    0.0D+00 < P.
    !
    !    Input, real ( kind = 8 ) X, the integral limit parameter.
    !    If X is less than or equal to 0, GAMMA_INC is returned as 0.
    !
    !    Output, real ( kind = 8 ) GAMMA_INC, the value of the function.
    !
    implicit none

    real    ( kind = 8 ) a
    real    ( kind = 8 ) arg
    real    ( kind = 8 ) b
    real    ( kind = 8 ) c
    real    ( kind = 8 ) cdf
    real    ( kind = 8 ), parameter :: exp_arg_min = -88.0D+00
    real    ( kind = 8 ) gamma_inc
    !real    ( kind = 8 ) gamma_log
    real    ( kind = 8 ), parameter :: overflow = 1.0D+37
    real    ( kind = 8 ) p
    real    ( kind = 8 ), parameter :: plimit = 1000.0D+00
    real    ( kind = 8 ) pn1
    real    ( kind = 8 ) pn2
    real    ( kind = 8 ) pn3
    real    ( kind = 8 ) pn4
    real    ( kind = 8 ) pn5
    real    ( kind = 8 ) pn6
    real    ( kind = 8 ) rn
    real    ( kind = 8 ), parameter :: tol = 1.0D-07
    real    ( kind = 8 ) x
    real    ( kind = 8 ), parameter :: xbig = 1.0D+08

    gamma_inc = 0.0D+00

    if ( p <= 0.0D+00 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'GAMMA_INC - Fatal error!'
       write ( *, '(a)' ) '  Parameter P <= 0.'
       stop
    end if

    if ( x <= 0.0D+00 ) then
       gamma_inc = 0.0D+00
       return
    end if
    !
    !  Use a normal approximation if PLIMIT < P.
    !
    if ( plimit < p ) then
       pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p ) ** ( 1.0D+00 / 3.0D+00 ) &
            + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )
       call normal_01_cdf ( pn1, cdf )
       gamma_inc = cdf
       return
    end if
    !
    !  Is X extremely large compared to P?
    !
    if ( xbig < x ) then
       gamma_inc = 1.0D+00
       return
    end if
    !
    !  Use Pearson's series expansion.
    !  (P is not large enough to force overflow in the log of Gamma.
    !
    if ( x <= 1.0D+00 .or. x < p ) then

       arg = p * log ( x ) - x - gamma_log ( p + 1.0D+00 )
       c = 1.0D+00
       gamma_inc = 1.0D+00
       a = p

       do

          a = a + 1.0D+00
          c = c * x / a
          gamma_inc = gamma_inc + c

          if ( c <= tol ) then
             exit
          end if

       end do

       arg = arg + log ( gamma_inc )

       if ( exp_arg_min <= arg ) then
          gamma_inc = exp ( arg )
       else
          gamma_inc = 0.0D+00
       end if

    else
       !
       !  Use a continued fraction expansion.
       !
       arg = p * log ( x ) - x - gamma_log ( p )
       a = 1.0D+00 - p
       b = a + x + 1.0D+00
       c = 0.0D+00
       pn1 = 1.0D+00
       pn2 = x
       pn3 = x + 1.0D+00
       pn4 = x * b
       gamma_inc = pn3 / pn4

       do

          a = a + 1.0D+00
          b = b + 2.0D+00
          c = c + 1.0D+00
          pn5 = b * pn3 - a * c * pn1
          pn6 = b * pn4 - a * c * pn2

          if ( 0.0D+00 < abs ( pn6 ) ) then

             rn = pn5 / pn6

             if ( abs ( gamma_inc - rn ) <= min ( tol, tol * rn ) ) then

                arg = arg + log ( gamma_inc )

                if ( exp_arg_min <= arg ) then
                   gamma_inc = 1.0D+00 - exp ( arg )
                else
                   gamma_inc = 1.0D+00
                end if

                return

             end if

             gamma_inc = rn

          end if

          pn1 = pn3
          pn2 = pn4
          pn3 = pn5
          pn4 = pn6
          !
          !  Rescale terms in continued fraction if terms are large.
          !
          if ( overflow <= abs ( pn5 ) ) then
             pn1 = pn1 / overflow
             pn2 = pn2 / overflow
             pn3 = pn3 / overflow
             pn4 = pn4 / overflow
          end if

       end do

    end if

    return
  end function gamma_inc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: BINOMIALE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function binomial_pdf (x,a,b) !x=k;a=n,b=p
    !      PDF(A,B;X) = C(N,X) * B**X * ( 1.0D+00 - B )**( A - X )
    !    Binomial_PDF(1,B;X) = Bernoulli_PDF(B;X).

    implicit none
    real, intent(in) :: b
    integer, intent(in) :: a,x
    real :: cnk
    real :: binomial_pdf

    if ( a < 1 ) then
       binomial_pdf = 0.0D+00
    else if ( x < 0 .or. a < x ) then
       binomial_pdf = 0.0D+00
    else if ( b == 0.0D+00 ) then
       if ( x == 0 ) then
          binomial_pdf = 1.0D+00
       else
          binomial_pdf = 0.0D+00
       end if
    else if ( b == 1.0D+00 ) then
       if ( x == a ) then
          binomial_pdf = 1.0D+00
       else
          binomial_pdf = 0.0D+00
       end if
    else

       cnk= binomial_coef(a , x)
       binomial_pdf = cnk * b**x * ( 1.0D+00 - b )**( a - x )

    end if

  end function binomial_pdf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: BINOMIALE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function beta_pdf ( x, a, b)

    !
    !      PDF(A,B;X) = X**(A-1) * (1-X)**(B-1) / BETA(A,B).
    !
    !    A = B = 1 yields the Uniform distribution on [0,1].
    !    A = B = 1/2 yields the Arcsin distribution.
    !        B = 1 yields the power function distribution.
    !    A = B -> Infinity tends to the Normal distribution.

    !    John Burkardt
    !  Parameters:
    !    Input, real ( kind = 8 ) X, the argument of the PDF.
    !    0.0D+00 <= X <= 1.0.
    !    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
    !    0.0D+00 < A,
    !    0.0D+00 < B.
    !    Output, real ( kind = 8 ) PDF, the value of the PDF.
    !
    implicit none

    real ( kind = 8 ), intent(in) :: x,a,b
    real :: beta_pdf

    if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
       beta_pdf = 0.0D+00
    else
       beta_pdf = x**( a - 1.0D+00 ) * ( 1.0D+00 - x )**( b - 1.0D+00 ) / beta ( a, b )
    end if

    return
  end function beta_pdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION random_beta(aa, bb, first) RESULT(fn_val)
    
    ! Adapted from Fortran 77 code from the book:
    !     Dagpunar, J. 'Principles of random variate generation'
    !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    
    ! FUNCTION GENERATES A RANDOM VARIATE IN [0,1]
    ! FROM A BETA DISTRIBUTION WITH DENSITY
    ! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
    ! USING CHENG'S LOG LOGISTIC METHOD.
    
    !     AA = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
    !     BB = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
    
    REAL, INTENT(IN)    :: aa, bb
    LOGICAL, INTENT(IN) :: first
    REAL                :: fn_val
    
    !     Local variables
    REAL, PARAMETER  :: aln4 = 1.3862944
    REAL             :: a, b, g, r, s, x, y, z , vsmall = TINY(1.0), vlarge = HUGE(1.0),zero = 0.0, one = 1.0, two = 2.0
    REAL, SAVE       :: d, f, h, t, c
    LOGICAL, SAVE    :: swap
    
    IF (aa <= zero .OR. bb <= zero) THEN
       WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE(S)'
       STOP
    END IF
    
    IF (first) THEN                        ! Initialization, if necessary
       a = aa
       b = bb
       swap = b > a
       IF (swap) THEN
          g = b
          b = a
          a = g
       END IF
       d = a/b
       f = a+b
       IF (b > one) THEN
          h = SQRT((two*a*b - f)/(f - two))
          t = one
       ELSE
          h = b
          t = one/(one + (a/(vlarge*b))**b)
       END IF
       c = a+h
    END IF
    
    DO
       CALL RANDOM_NUMBER(r)
       CALL RANDOM_NUMBER(x)
       !  print *,r,x
       s = r*r*x
       IF (r < vsmall .OR. s <= zero) CYCLE
       IF (r < t) THEN
          x = LOG(r/(one - r))/h
          y = d*EXP(x)
          z = c*x + f*LOG((one + d)/(one + y)) - aln4
          IF (s - one > z) THEN
             IF (s - s*z > one) CYCLE
             IF (LOG(s) > z) CYCLE
          END IF
          fn_val = y/(one + y)
       ELSE
          IF (4.0*s > (one + one/d)**f) CYCLE
          fn_val = one
       END IF
       EXIT
    END DO
    
    IF (swap) fn_val = one - fn_val
    RETURN
  END FUNCTION random_beta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION random_binomial2(n, pp, first) RESULT(ival)
    !**********************************************************************
    !     Translated to Fortran 90 by Alan Miller from:
    !                              RANLIB
    !
    !     Library of Fortran Routines for Random Number Generation
    !
    !                      Compiled and Written by:
    !
    !                           Barry W. Brown
    !                            James Lovato
    !
    !               Department of Biomathematics, Box 237
    !               The University of Texas, M.D. Anderson Cancer Center
    !               1515 Holcombe Boulevard
    !               Houston, TX      77030
    !
    ! This work was supported by grant CA-16672 from the National Cancer Institute.
    
    !                    GENerate BINomial random deviate
    
    !                              Function
    
    !     Generates a single random deviate from a binomial
    !     distribution whose number of trials is N and whose
    !     probability of an event in each trial is P.
    
    !                              Arguments
    
    !     N  --> The number of trials in the binomial distribution
    !            from which a random deviate is to be generated.
    !                              INTEGER N
    
    !     P  --> The probability of an event in each trial of the
    !            binomial distribution from which a random deviate
    !            is to be generated.
    !                              REAL P
    
    !     FIRST --> Set FIRST = .TRUE. for the first call to perform initialization
    !               the set FIRST = .FALSE. for further calls using the same pair
    !               of parameter values (N, P).
    !                              LOGICAL FIRST
    
    !     random_binomial2 <-- A random deviate yielding the number of events
    !                from N independent trials, each of which has
    !                a probability of event P.
    !                              INTEGER random_binomial
    
    !                              Method
    
    !     This is algorithm BTPE from:
    
    !         Kachitvichyanukul, V. and Schmeiser, B. W.
    !         Binomial Random Variate Generation.
    !         Communications of the ACM, 31, 2 (February, 1988) 216.
    
    !**********************************************************************
    
    !*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
    
    !     ..
    !     .. Scalar Arguments ..
    REAL, INTENT(IN)    :: pp
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: first
    INTEGER             :: ival
    !     ..
    !     .. Local Scalars ..
    REAL            :: alv, amaxp, f, f1, f2, u, v, w, w2, x, x1, x2, ynorm, z, z2
    REAL, PARAMETER :: zero = 0.0, half = 0.5, one = 1.0
    INTEGER         :: i, ix, ix1, k, mp
    INTEGER, SAVE   :: m
    REAL, SAVE      :: p, q, xnp, ffm, fm, xnpq, p1, xm, xl, xr, c, al, xll,  &
    xlr, p2, p3, p4, qn, r, g
    
    !     ..
    !     .. Executable Statements ..
    
    !*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
    
    IF (first) THEN
       p = MIN(pp, one-pp)
       q = one - p
       xnp = n * p
    END IF
    
    IF (xnp > 30.) THEN
       IF (first) THEN
          ffm = xnp + p
          m = ffm
          fm = m
          xnpq = xnp * q
          p1 = INT(2.195*SQRT(xnpq) - 4.6*q) + half
          xm = fm + half
          xl = xm - p1
          xr = xm + p1
          c = 0.134 + 20.5 / (15.3 + fm)
          al = (ffm-xl) / (ffm - xl*p)
          xll = al * (one + half*al)
          al = (xr - ffm) / (xr*q)
          xlr = al * (one + half*al)
          p2 = p1 * (one + c + c)
          p3 = p2 + c / xll
          p4 = p3 + c / xlr
       END IF
       
       !*****GENERATE VARIATE, Binomial mean at least 30.
       
20       u=grnd() !CALL RANDOM_NUMBER(u)
       u = u * p4
       v=grnd() !CALL RANDOM_NUMBER(v)
       
       !     TRIANGULAR REGION
       
       IF (u <= p1) THEN
          ix = xm - p1 * v + u
          GO TO 110
       END IF
       
       !     PARALLELOGRAM REGION
       
       IF (u <= p2) THEN
          x = xl + (u-p1) / c
          v = v * c + one - ABS(xm-x) / p1
          IF (v > one .OR. v <= zero) GO TO 20
          ix = x
       ELSE
          
          !     LEFT TAIL
          
          IF (u <= p3) THEN
             ix = xl + LOG(v) / xll
             IF (ix < 0) GO TO 20
             v = v * (u-p2) * xll
          ELSE
             
             !     RIGHT TAIL
             
             ix = xr - LOG(v) / xlr
             IF (ix > n) GO TO 20
             v = v * (u-p3) * xlr
          END IF
       END IF
       
       !*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
       
       k = ABS(ix-m)
       IF (k <= 20 .OR. k >= xnpq/2-1) THEN
          
          !     EXPLICIT EVALUATION
          
          f = one
          r = p / q
          g = (n+1) * r
          IF (m < ix) THEN
             mp = m + 1
             DO i = mp, ix
                f = f * (g/i-r)
             END DO
             
          ELSE IF (m > ix) THEN
             ix1 = ix + 1
             DO i = ix1, m
                f = f / (g/i-r)
             END DO
          END IF
          
          IF (v > f) THEN
             GO TO 20
          ELSE
             GO TO 110
          END IF
       END IF
       
       !     SQUEEZING USING UPPER AND LOWER BOUNDS ON LOG(F(X))
       
       amaxp = (k/xnpq) * ((k*(k/3. + .625) + .1666666666666)/xnpq + half)
       ynorm = -k * k / (2.*xnpq)
       alv = LOG(v)
       IF (alv<ynorm - amaxp) GO TO 110
       IF (alv>ynorm + amaxp) GO TO 20
       
       !     STIRLING'S (actually de Moivre's) FORMULA TO MACHINE ACCURACY FOR
       !     THE FINAL ACCEPTANCE/REJECTION TEST
       
       x1 = ix + 1
       f1 = fm + one
       z = n + 1 - fm
       w = n - ix + one
       z2 = z * z
       x2 = x1 * x1
       f2 = f1 * f1
       w2 = w * w
       IF (alv - (xm*LOG(f1/x1) + (n-m+half)*LOG(z/w) + (ix-m)*LOG(w*p/(x1*q)) +    &
          (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320. +               &
          (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320. +                &
          (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320. +               &
          (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.) > zero) THEN
          GO TO 20
       ELSE
          GO TO 110
       END IF
       
    ELSE
       !     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
       IF (first) THEN
          qn = q ** n
          r = p / q
          g = r * (n+1)
       END IF
       
90       ix = 0
       f = qn
       u=grnd() !CALL RANDOM_NUMBER(u)
100      IF (u >= f) THEN
          IF (ix > 110) GO TO 90
          u = u - f
          ix = ix + 1
          f = f * (g/ix - r)
          GO TO 100
       END IF
    END IF
    
110   IF (pp > half) ix = n - ix
    ival = ix
    RETURN
    
  END FUNCTION random_binomial2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module prob_mod




