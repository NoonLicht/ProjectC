#include "MySourceLib.h"

using namespace std;

//####################################################################
struct simp {
 int n;
 double cp;
 double cko;
 double *x;
 int *r;
};
simp sm;

struct simpf{
 int k;
 int *n;
 double *x;
 double *p;
 string ts;
};
simpf smf;

double algdiv ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
//
//  Discussion:
//
//    In this algorithm, DEL(X) is the function defined by
//
//      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI )
//                      + DEL(X).
//
//  Parameters:
//
//    Input, double *A, *B, define the arguments.
//
//    Output, double ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
//
{
  static double algdiv;
  static double c;
  static double c0 =  0.833333333333333e-01;
  static double c1 = -0.277777777760991e-02;
  static double c2 =  0.793650666825390e-03;
  static double c3 = -0.595202931351870e-03;
  static double c4 =  0.837308034031215e-03;
  static double c5 = -0.165322962780713e-02;
  static double d;
  static double h;
  static double s11;
  static double s3;
  static double s5;
  static double s7;
  static double s9;
  static double t;
  static double T1;
  static double u;
  static double v;
  static double w;
  static double x;
  static double x2;

  if ( *b <= *a )
  {
    h = *b / *a;
    c = 1.0e0 / ( 1.0e0 + h );
    x = h / ( 1.0e0 + h );
    d = *a + ( *b - 0.5e0 );
  }
  else
  {
    h = *a / *b;
    c = h / ( 1.0e0 + h );
    x = 1.0e0 / ( 1.0e0 + h );
    d = *b + ( *a - 0.5e0 );
  }
//
//  SET SN = (1 - X**N)/(1 - X)
//
  x2 = x * x;
  s3 = 1.0e0 + ( x + x2 );
  s5 = 1.0e0 + ( x + x2 * s3 );
  s7 = 1.0e0 + ( x + x2 * s5 );
  s9 = 1.0e0 + ( x + x2 * s7 );
  s11 = 1.0e0 + ( x + x2 * s9 );
//
//  SET W = DEL(B) - DEL(A + B)
//
  t = pow ( 1.0e0 / *b, 2.0 );

  w = (((( c5 * s11  * t
         + c4 * s9 ) * t
         + c3 * s7 ) * t
         + c2 * s5 ) * t
         + c1 * s3 ) * t
         + c0;

  w *= ( c / *b );
//
//  Combine the results.
//
  T1 = *a / *b;
  u = d * alnrel ( &T1 );
  v = *a * ( log ( *b ) - 1.0e0 );

  if ( v < u )
  {
    algdiv = w - v - u;
  }
  else
  {
    algdiv = w - u - v;
  }
  return algdiv;
}
//****************************************************************************80

double alnrel ( double *a )

//****************************************************************************80
//
//  Purpose:
//
//    ALNREL evaluates the function ln ( 1 + A ).
//
//  Modified:
//
//    17 November 2006
//
//  Reference:
//
//    Armido DiDinato, Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *A, the argument.
//
//    Output, double ALNREL, the value of ln ( 1 + A ).
//
{
  double alnrel;
  static double p1 = -0.129418923021993e+01;
  static double p2 =  0.405303492862024e+00;
  static double p3 = -0.178874546012214e-01;
  static double q1 = -0.162752256355323e+01;
  static double q2 =  0.747811014037616e+00;
  static double q3 = -0.845104217945565e-01;
  double t;
  double t2;
  double w;
  double x;

  if ( fabs ( *a ) <= 0.375e0 )
  {
    t = *a / ( *a + 2.0e0 );
    t2 = t * t;
    w = (((p3*t2+p2)*t2+p1)*t2+1.0e0)
      / (((q3*t2+q2)*t2+q1)*t2+1.0e0);
    alnrel = 2.0e0 * t * w;
  }
  else
  {
    x = 1.0e0 + *a;
    alnrel = log ( x );
  }
  return alnrel;
}
//****************************************************************************80

double apser ( double *a, double *b, double *x, double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
//
//  Discussion:
//
//    APSER is used only for cases where
//
//      A <= min ( EPS, EPS * B ),
//      B * X <= 1, and
//      X <= 0.5.
//
//  Parameters:
//
//    Input, double *A, *B, *X, the parameters of teh
//    incomplete beta ratio.
//
//    Input, double *EPS, a tolerance.
//
//    Output, double APSER, the computed value of the
//    incomplete beta ratio.
//
{
  static double g = 0.577215664901533e0;
  static double apser,aj,bx,c,j,s,t,tol;

    bx = *b**x;
    t = *x-bx;
    if(*b**eps > 2.e-2) goto S10;
    c = log(*x)+psi(b)+g+t;
    goto S20;
S10:
    c = log(bx)+g+t;
S20:
    tol = 5.0e0**eps*fabs(c);
    j = 1.0e0;
    s = 0.0e0;
S30:
    j = j + 1.0e0;
    t = t * (*x-bx/j);
    aj = t/j;
    s = s + aj;
    if(fabs(aj) > tol) goto S30;
    apser = -(*a*(c+s));
    return apser;
}
//****************************************************************************80

double bcorr ( double *a0, double *b0 )

//****************************************************************************80
//
//  Purpose:
//
//    BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
//
//  Discussion:
//
//    The function DEL(A) is a remainder term that is used in the expression:
//
//      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A )
//        - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
//
//    or, in other words, DEL ( A ) is defined as:
//
//      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A )
//        + A + 0.5 * ln ( 2 * PI ).
//
//  Parameters:
//
//    Input, double *A0, *B0, the arguments.
//    It is assumed that 8 <= A0 and 8 <= B0.
//
//    Output, double *BCORR, the value of the function.
//
{
  static double c0 =  0.833333333333333e-01;
  static double c1 = -0.277777777760991e-02;
  static double c2 =  0.793650666825390e-03;
  static double c3 = -0.595202931351870e-03;
  static double c4 =  0.837308034031215e-03;
  static double c5 = -0.165322962780713e-02;
  static double bcorr,a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;

  a = fifdmin1 ( *a0, *b0 );
  b = fifdmax1 ( *a0, *b0 );
  h = a / b;
  c = h / ( 1.0e0 + h );
  x = 1.0e0 / ( 1.0e0 + h );
  x2 = x * x;
//
//  SET SN = (1 - X**N)/(1 - X)
//
  s3 = 1.0e0 + ( x + x2 );
  s5 = 1.0e0 + ( x + x2 * s3 );
  s7 = 1.0e0 + ( x + x2 * s5 );
  s9 = 1.0e0 + ( x + x2 * s7 );
  s11 = 1.0e0 + ( x + x2 * s9 );
//
//  SET W = DEL(B) - DEL(A + B)
//
  t = pow ( 1.0e0 / b, 2.0 );

  w = (((( c5 * s11  * t + c4
              * s9 ) * t + c3
              * s7 ) * t + c2
              * s5 ) * t + c1
              * s3 ) * t + c0;
  w *= ( c / b );
//
//  COMPUTE  DEL(A) + W
//
  t = pow ( 1.0e0 / a, 2.0 );

  bcorr = ((((( c5 * t + c4 )
                   * t + c3 )
                   * t + c2 )
                   * t + c1 )
                   * t + c0 ) / a + w;
  return bcorr;
}
//****************************************************************************80

double beta ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA evaluates the beta function.
//
//  Modified:
//
//    03 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the arguments of the beta function.
//
//    Output, double BETA, the value of the beta function.
//
{
  return ( exp ( beta_log ( &a, &b ) ) );
}
//****************************************************************************80

double beta_asym ( double *a, double *b, double *lambda, double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that both A and B
//    are greater than or equal to 15.
//
//    Input, double *LAMBDA, the value of ( A + B ) * Y - B.
//    It is assumed that 0 <= LAMBDA.
//
//    Input, double *EPS, the tolerance.
//
{
  static double e0 = 1.12837916709551e0;
  static double e1 = .353553390593274e0;
  static int num = 20;
//
//  NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
//            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
//            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
//     E0 = 2/SQRT(PI)
//     E1 = 2**(-3/2)
//
  static int K3 = 1;
  static double value;
  static double bsum,dsum,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,t0,t1,u,w,w0,z,z0,
    z2,zn,znm1;
  static int i,im1,imj,j,m,mm1,mmj,n,np1;
  static double a0[21],b0[21],c[21],d[21],T1,T2;

    value = 0.0e0;
    if(*a >= *b) goto S10;
    h = *a/ *b;
    r0 = 1.0e0/(1.0e0+h);
    r1 = (*b-*a)/ *b;
    w0 = 1.0e0/sqrt(*a*(1.0e0+h));
    goto S20;
S10:
    h = *b/ *a;
    r0 = 1.0e0/(1.0e0+h);
    r1 = (*b-*a)/ *a;
    w0 = 1.0e0/sqrt(*b*(1.0e0+h));
S20:
    T1 = -(*lambda/ *a);
    T2 = *lambda/ *b;
    f = *a*rlog1(&T1)+*b*rlog1(&T2);
    t = exp(-f);
    if(t == 0.0e0) return value;
    z0 = sqrt(f);
    z = 0.5e0*(z0/e1);
    z2 = f+f;
    a0[0] = 2.0e0/3.0e0*r1;
    c[0] = -(0.5e0*a0[0]);
    d[0] = -c[0];
    j0 = 0.5e0/e0 * error_fc ( &K3, &z0 );
    j1 = e1;
    sum = j0+d[0]*w0*j1;
    s = 1.0e0;
    h2 = h*h;
    hn = 1.0e0;
    w = w0;
    znm1 = z;
    zn = z2;
    for ( n = 2; n <= num; n += 2 )
    {
        hn = h2*hn;
        a0[n-1] = 2.0e0*r0*(1.0e0+h*hn)/((double)n+2.0e0);
        np1 = n+1;
        s += hn;
        a0[np1-1] = 2.0e0*r1*s/((double)n+3.0e0);
        for ( i = n; i <= np1; i++ )
        {
            r = -(0.5e0*((double)i+1.0e0));
            b0[0] = r*a0[0];
            for ( m = 2; m <= i; m++ )
            {
                bsum = 0.0e0;
                mm1 = m-1;
                for ( j = 1; j <= mm1; j++ )
                {
                    mmj = m-j;
                    bsum += (((double)j*r-(double)mmj)*a0[j-1]*b0[mmj-1]);
                }
                b0[m-1] = r*a0[m-1]+bsum/(double)m;
            }
            c[i-1] = b0[i-1]/((double)i+1.0e0);
            dsum = 0.0e0;
            im1 = i-1;
            for ( j = 1; j <= im1; j++ )
            {
                imj = i-j;
                dsum += (d[imj-1]*c[j-1]);
            }
            d[i-1] = -(dsum+c[i-1]);
        }
        j0 = e1*znm1+((double)n-1.0e0)*j0;
        j1 = e1*zn+(double)n*j1;
        znm1 = z2*znm1;
        zn = z2*zn;
        w = w0*w;
        t0 = d[n-1]*w*j0;
        w = w0*w;
        t1 = d[np1-1]*w*j1;
        sum += (t0+t1);
        if(fabs(t0)+fabs(t1) <= *eps*sum) goto S80;
    }
S80:
    u = exp(-bcorr(a,b));
    value = e0*t*u*sum;
    return value;
}
//****************************************************************************80

double beta_frac ( double *a, double *b, double *x, double *y, double *lambda, double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that both A and
//    B are greater than 1.
//
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Input, double *LAMBDA, the value of ( A + B ) * Y - B.
//
//    Input, double *EPS, a tolerance.
//
//    Output, double BETA_FRAC, the value of the continued
//    fraction approximation for IX(A,B).
//
{
  static double bfrac,alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,t,w,yp1;

  bfrac = beta_rcomp ( a, b, x, y );

  if ( bfrac == 0.0e0 )
  {
    return bfrac;
  }

  c = 1.0e0+*lambda;
  c0 = *b/ *a;
  c1 = 1.0e0+1.0e0/ *a;
  yp1 = *y+1.0e0;
  n = 0.0e0;
  p = 1.0e0;
  s = *a+1.0e0;
  an = 0.0e0;
  bn = anp1 = 1.0e0;
  bnp1 = c/c1;
  r = c1/c;
//
//  CONTINUED FRACTION CALCULATION
//
S10:
  n = n + 1.0e0;
  t = n/ *a;
  w = n*(*b-n)**x;
  e = *a/s;
  alpha = p*(p+c0)*e*e*(w**x);
  e = (1.0e0+t)/(c1+t+t);
  beta = n+w/s+e*(c+n*yp1);
  p = 1.0e0+t;
  s += 2.0e0;
//
//  UPDATE AN, BN, ANP1, AND BNP1
//
  t = alpha*an+beta*anp1;
  an = anp1;
  anp1 = t;
  t = alpha*bn+beta*bnp1;
  bn = bnp1;
  bnp1 = t;
  r0 = r;
  r = anp1/bnp1;

  if ( fabs(r-r0) <= (*eps) * r )
  {
    goto S20;
  }
//
//  RESCALE AN, BN, ANP1, AND BNP1
//
  an /= bnp1;
  bn /= bnp1;
  anp1 = r;
  bnp1 = 1.0e0;
  goto S10;
//
//  TERMINATION
//
S20:
  bfrac = bfrac * r;
  return bfrac;
}
//****************************************************************************80

void beta_grat ( double *a, double *b, double *x, double *y, double *w,double *eps,int *ierr )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that 15 <= A
//    and B <= 1, and that B is less than A.
//
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Input/output, double *W, a quantity to which the
//    result of the computation is to be added on output.
//
//    Input, double *EPS, a tolerance.
//
//    Output, int *IERR, an error flag, which is 0 if no error
//    was detected.
//
{
  static double bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,t2,u,v,z;
  static int i,n,nm1;
  static double c[30],d[30],T1;

    bm1 = *b-0.5e0-0.5e0;
    nu = *a+0.5e0*bm1;
    if(*y > 0.375e0) goto S10;
    T1 = -*y;
    lnx = alnrel(&T1);
    goto S20;
S10:
    lnx = log(*x);
S20:
    z = -(nu*lnx);
    if(*b*z == 0.0e0) goto S70;
//
//  COMPUTATION OF THE EXPANSION
//  SET R = EXP(-Z)*Z**B/GAMMA(B)
//
    r = *b*(1.0e0+gam1(b))*exp(*b*log(z));
    r *= (exp(*a*lnx)*exp(0.5e0*bm1*lnx));
    u = algdiv(b,a)+*b*log(nu);
    u = r*exp(-u);
    if(u == 0.0e0) goto S70;
    gamma_rat1 ( b, &z, &r, &p, &q, eps );
    v = 0.25e0*pow(1.0e0/nu,2.0);
    t2 = 0.25e0*lnx*lnx;
    l = *w/u;
    j = q/r;
    sum = j;
    t = cn = 1.0e0;
    n2 = 0.0e0;
    for ( n = 1; n <= 30; n++ )
    {
        bp2n = *b+n2;
        j = (bp2n*(bp2n+1.0e0)*j+(z+bp2n+1.0e0)*t)*v;
        n2 = n2 + 2.0e0;
        t *= t2;
        cn /= (n2*(n2+1.0e0));
        c[n-1] = cn;
        s = 0.0e0;
        if(n == 1) goto S40;
        nm1 = n-1;
        coef = *b-(double)n;
        for ( i = 1; i <= nm1; i++ )
        {
            s = s + (coef*c[i-1]*d[n-i-1]);
            coef = coef + *b;
        }
S40:
        d[n-1] = bm1*cn+s/(double)n;
        dj = d[n-1]*j;
        sum = sum + dj;
        if(sum <= 0.0e0) goto S70;
        if(fabs(dj) <= *eps*(sum+l)) goto S60;
    }
S60:
//
//  ADD THE RESULTS TO W
//
    *ierr = 0;
    *w = *w + (u*sum);
    return;
S70:
//
//  THE EXPANSION CANNOT BE COMPUTED
//
    *ierr = 1;
    return;
}
//****************************************************************************80

void beta_inc ( double *a, double *b, double *x, double *y, double *w, double *w1, int *ierr )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC evaluates the incomplete beta function IX(A,B).
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.
//
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Output, double *W, *W1, the values of IX(A,B) and
//    1-IX(A,B).
//
//    Output, int *IERR, the error flag.
//    0, no error was detected.
//    1, A or B is negative;
//    2, A = B = 0;
//    3, X < 0 or 1 < X;
//    4, Y < 0 or 1 < Y;
//    5, X + Y /= 1;
//    6, X = A = 0;
//    7, Y = B = 0.
//
{
  static int K1 = 1;
  static double a0,b0,eps,lambda,t,x0,y0,z;
  static int ierr1,ind,n;
  static double T2,T3,T4,T5;
//
//  EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
//  NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
//
    eps = dpmpar ( &K1 );
    *w = *w1 = 0.0e0;
    if(*a < 0.0e0 || *b < 0.0e0) goto S270;
    if(*a == 0.0e0 && *b == 0.0e0) goto S280;
    if(*x < 0.0e0 || *x > 1.0e0) goto S290;
    if(*y < 0.0e0 || *y > 1.0e0) goto S300;
    z = *x+*y-0.5e0-0.5e0;
    if(fabs(z) > 3.0e0*eps) goto S310;
    *ierr = 0;
    if(*x == 0.0e0) goto S210;
    if(*y == 0.0e0) goto S230;
    if(*a == 0.0e0) goto S240;
    if(*b == 0.0e0) goto S220;
    eps = fifdmax1(eps,1.e-15);
    if(fifdmax1(*a,*b) < 1.e-3*eps) goto S260;
    ind = 0;
    a0 = *a;
    b0 = *b;
    x0 = *x;
    y0 = *y;
    if(fifdmin1(a0,b0) > 1.0e0) goto S40;
//
//  PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
//
    if(*x <= 0.5e0) goto S10;
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
S10:
    if(b0 < fifdmin1(eps,eps*a0)) goto S90;
    if(a0 < fifdmin1(eps,eps*b0) && b0*x0 <= 1.0e0) goto S100;
    if(fifdmax1(a0,b0) > 1.0e0) goto S20;
    if(a0 >= fifdmin1(0.2e0,b0)) goto S110;
    if(pow(x0,a0) <= 0.9e0) goto S110;
    if(x0 >= 0.3e0) goto S120;
    n = 20;
    goto S140;
S20:
    if(b0 <= 1.0e0) goto S110;
    if(x0 >= 0.3e0) goto S120;
    if(x0 >= 0.1e0) goto S30;
    if(pow(x0*b0,a0) <= 0.7e0) goto S110;
S30:
    if(b0 > 15.0e0) goto S150;
    n = 20;
    goto S140;
S40:
//
//  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
//
    if(*a > *b) goto S50;
    lambda = *a-(*a+*b)**x;
    goto S60;
S50:
    lambda = (*a+*b)**y-*b;
S60:
    if(lambda >= 0.0e0) goto S70;
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
    lambda = fabs(lambda);
S70:
    if(b0 < 40.0e0 && b0*x0 <= 0.7e0) goto S110;
    if(b0 < 40.0e0) goto S160;
    if(a0 > b0) goto S80;
    if(a0 <= 100.0e0) goto S130;
    if(lambda > 0.03e0*a0) goto S130;
    goto S200;
S80:
    if(b0 <= 100.0e0) goto S130;
    if(lambda > 0.03e0*b0) goto S130;
    goto S200;
S90:
//
//  EVALUATION OF THE APPROPRIATE ALGORITHM
//
    *w = fpser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S100:
    *w1 = apser(&a0,&b0,&x0,&eps);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S110:
    *w = beta_pser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S120:
    *w1 = beta_pser(&b0,&a0,&y0,&eps);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S130:
    T2 = 15.0e0*eps;
    *w = beta_frac ( &a0,&b0,&x0,&y0,&lambda,&T2 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S140:
    *w1 = beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
    b0 = b0 + (double)n;
S150:
    T3 = 15.0e0*eps;
    beta_grat (&b0,&a0,&y0,&x0,w1,&T3,&ierr1);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S160:
    n = ( int ) b0;
    b0 -= (double)n;
    if(b0 != 0.0e0) goto S170;
    n -= 1;
    b0 = 1.0e0;
S170:
    *w = beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
    if(x0 > 0.7e0) goto S180;
    *w = *w + beta_pser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S180:
    if(a0 > 15.0e0) goto S190;
    n = 20;
    *w = *w + beta_up ( &a0, &b0, &x0, &y0, &n, &eps );
    a0 = a0 + (double)n;
S190:
    T4 = 15.0e0*eps;
    beta_grat ( &a0, &b0, &x0, &y0, w, &T4, &ierr1 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S200:
    T5 = 100.0e0*eps;
    *w = beta_asym ( &a0, &b0, &lambda, &T5 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S210:
//
//  TERMINATION OF THE PROCEDURE
//
    if(*a == 0.0e0) goto S320;
S220:
    *w = 0.0e0;
    *w1 = 1.0e0;
    return;
S230:
    if(*b == 0.0e0) goto S330;
S240:
    *w = 1.0e0;
    *w1 = 0.0e0;
    return;
S250:
    if(ind == 0) return;
    t = *w;
    *w = *w1;
    *w1 = t;
    return;
S260:
//
//  PROCEDURE FOR A AND B .LT. 1.E-3*EPS
//
    *w = *b/(*a+*b);
    *w1 = *a/(*a+*b);
    return;
S270:
//
//  ERROR RETURN
//
    *ierr = 1;
    return;
S280:
    *ierr = 2;
    return;
S290:
    *ierr = 3;
    return;
S300:
    *ierr = 4;
    return;
S310:
    *ierr = 5;
    return;
S320:
    *ierr = 6;
    return;
S330:
    *ierr = 7;
    return;
}
//****************************************************************************80

void beta_inc_values ( int *n_data, double *a, double *b, double *x,double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0
//      BETA_INC(A,B,1.0) = 1.0
//
//    Note that in Mathematica, the expressions:
//
//      BETA[A,B]   = Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//      BETA[X,A,B] = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//
//    and thus, to evaluate the incomplete Beta function requires:
//
//      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
//
//  Modified:
//
//    09 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, *B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 30

  double a_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,  1.0E+00,
     1.0E+00,  1.0E+00,  1.0E+00,  1.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  5.5E+00, 10.0E+00, 10.0E+00,
    10.0E+00, 10.0E+00, 20.0E+00, 20.0E+00,
    20.0E+00, 20.0E+00, 20.0E+00, 30.0E+00,
    30.0E+00, 40.0E+00 };
  double b_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,  0.5E+00,
     0.5E+00,  0.5E+00,  0.5E+00,  1.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  5.0E+00,  0.5E+00,  5.0E+00,
     5.0E+00, 10.0E+00,  5.0E+00, 10.0E+00,
    10.0E+00, 20.0E+00, 20.0E+00, 10.0E+00,
    10.0E+00, 20.0E+00 };
  double fx_vec[N_MAX] = {
    0.0637686E+00, 0.2048328E+00, 1.0000000E+00, 0.0E+00,
    0.0050126E+00, 0.0513167E+00, 0.2928932E+00, 0.5000000E+00,
    0.028E+00,     0.104E+00,     0.216E+00,     0.352E+00,
    0.500E+00,     0.648E+00,     0.784E+00,     0.896E+00,
    0.972E+00,     0.4361909E+00, 0.1516409E+00, 0.0897827E+00,
    1.0000000E+00, 0.5000000E+00, 0.4598773E+00, 0.2146816E+00,
    0.9507365E+00, 0.5000000E+00, 0.8979414E+00, 0.2241297E+00,
    0.7586405E+00, 0.7001783E+00 };
  double x_vec[N_MAX] = {
    0.01E+00, 0.10E+00, 1.00E+00, 0.0E+00,
    0.01E+00, 0.10E+00, 0.50E+00, 0.50E+00,
    0.1E+00,  0.2E+00,  0.3E+00,  0.4E+00,
    0.5E+00,  0.6E+00,  0.7E+00,  0.8E+00,
    0.9E+00,  0.50E+00, 0.90E+00, 0.50E+00,
    1.00E+00, 0.50E+00, 0.80E+00, 0.60E+00,
    0.80E+00, 0.50E+00, 0.60E+00, 0.70E+00,
    0.80E+00, 0.70E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *b = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

double beta_log ( double *a0, double *b0 )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_LOG evaluates the logarithm of the beta function.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *A0, *B0, the parameters of the function.
//    A0 and B0 should be nonnegative.
//
//    Output, double *BETA_LOG, the value of the logarithm
//    of the Beta function.
//
{
  static double e = .918938533204673e0;
  static double value,a,b,c,h,u,v,w,z;
  static int i,n;
  static double T1;

    a = fifdmin1(*a0,*b0);
    b = fifdmax1(*a0,*b0);
    if(a >= 8.0e0) goto S100;
    if(a >= 1.0e0) goto S20;
//
//  PROCEDURE WHEN A .LT. 1
//
    if(b >= 8.0e0) goto S10;
    T1 = a+b;
    value = gamma_log ( &a )+( gamma_log ( &b )- gamma_log ( &T1 ));
    return value;
S10:
    value = gamma_log ( &a )+algdiv(&a,&b);
    return value;
S20:
//
//  PROCEDURE WHEN 1 .LE. A .LT. 8
//
    if(a > 2.0e0) goto S40;
    if(b > 2.0e0) goto S30;
    value = gamma_log ( &a )+ gamma_log ( &b )-gsumln(&a,&b);
    return value;
S30:
    w = 0.0e0;
    if(b < 8.0e0) goto S60;
    value = gamma_log ( &a )+algdiv(&a,&b);
    return value;
S40:
//
//  REDUCTION OF A WHEN B .LE. 1000
//
    if(b > 1000.0e0) goto S80;
    n = ( int ) ( a - 1.0e0 );
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        a -= 1.0e0;
        h = a/b;
        w *= (h/(1.0e0+h));
    }
    w = log(w);
    if(b < 8.0e0) goto S60;
    value = w+ gamma_log ( &a )+algdiv(&a,&b);
    return value;
S60:
//
//  REDUCTION OF B WHEN B .LT. 8
//
    n = ( int ) ( b - 1.0e0 );
    z = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b -= 1.0e0;
        z *= (b/(a+b));
    }
    value = w+log(z)+( gamma_log ( &a )+( gamma_log ( &b )-gsumln(&a,&b)));
    return value;
S80:
//
//  REDUCTION OF A WHEN B .GT. 1000
//
    n = ( int ) ( a - 1.0e0 );
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        a -= 1.0e0;
        w *= (a/(1.0e0+a/b));
    }
    value = log(w)-(double)n*log(b)+( gamma_log ( &a )+algdiv(&a,&b));
    return value;
S100:
//
//  PROCEDURE WHEN A .GE. 8
//
    w = bcorr(&a,&b);
    h = a/b;
    c = h/(1.0e0+h);
    u = -((a-0.5e0)*log(c));
    v = b*alnrel(&h);
    if(u <= v) goto S110;
    value = -(0.5e0*log(b))+e+w-v-u;
    return value;
S110:
    value = -(0.5e0*log(b))+e+w-u-v;
    return value;
}
//****************************************************************************80

double beta_pser ( double *a, double *b, double *x, double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
//
//  Discussion:
//
//    BETA_PSER is used when B <= 1 or B*X <= 0.7.
//
//  Parameters:
//
//    Input, double *A, *B, the parameters.
//
//    Input, double *X, the point where the function
//    is to be evaluated.
//
//    Input, double *EPS, the tolerance.
//
//    Output, double BETA_PSER, the approximate value of IX(A,B)(X).
//
{
  static double bpser,a0,apb,b0,c,n,sum,t,tol,u,w,z;
  static int i,m;

    bpser = 0.0e0;
    if(*x == 0.0e0) return bpser;
//
//  COMPUTE THE FACTOR X**A/(A*BETA(A,B))
//
    a0 = fifdmin1(*a,*b);
    if(a0 < 1.0e0) goto S10;
    z = *a*log(*x)-beta_log(a,b);
    bpser = exp(z)/ *a;
    goto S100;
S10:
    b0 = fifdmax1(*a,*b);
    if(b0 >= 8.0e0) goto S90;
    if(b0 > 1.0e0) goto S40;
//
//  PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
//
    bpser = pow(*x,*a);
    if(bpser == 0.0e0) return bpser;
    apb = *a+*b;
    if(apb > 1.0e0) goto S20;
    z = 1.0e0+gam1(&apb);
    goto S30;
S20:
    u = *a+*b-1.e0;
    z = (1.0e0+gam1(&u))/apb;
S30:
    c = (1.0e0+gam1(a))*(1.0e0+gam1(b))/z;
    bpser *= (c*(*b/apb));
    goto S100;
S40:
//
//  PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
//
    u = gamma_ln1 ( &a0 );
    m = ( int ) ( b0 - 1.0e0 );
    if(m < 1) goto S60;
    c = 1.0e0;
    for ( i = 1; i <= m; i++ )
    {
        b0 -= 1.0e0;
        c *= (b0/(a0+b0));
    }
    u = log(c)+u;
S60:
    z = *a*log(*x)-u;
    b0 -= 1.0e0;
    apb = a0+b0;
    if(apb > 1.0e0) goto S70;
    t = 1.0e0+gam1(&apb);
    goto S80;
S70:
    u = a0+b0-1.e0;
    t = (1.0e0+gam1(&u))/apb;
S80:
    bpser = exp(z)*(a0/ *a)*(1.0e0+gam1(&b0))/t;
    goto S100;
S90:
//
//  PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
//
    u = gamma_ln1 ( &a0 ) + algdiv ( &a0, &b0 );
    z = *a*log(*x)-u;
    bpser = a0/ *a*exp(z);
S100:
    if(bpser == 0.0e0 || *a <= 0.1e0**eps) return bpser;
//
//  COMPUTE THE SERIES
//
    sum = n = 0.0e0;
    c = 1.0e0;
    tol = *eps/ *a;
S110:
    n = n + 1.0e0;
    c *= ((0.5e0+(0.5e0-*b/n))**x);
    w = c/(*a+n);
    sum = sum + w;
    if(fabs(w) > tol) goto S110;
    bpser *= (1.0e0+*a*sum);
    return bpser;
}
//****************************************************************************80

double beta_rcomp ( double *a, double *b, double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_RCOMP evaluates X**A * Y**B / Beta(A,B).
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the Beta function.
//    A and B should be nonnegative.
//
//    Input, double *X, *Y, define the numerator of the fraction.
//
//    Output, double BETA_RCOMP, the value of X**A * Y**B / Beta(A,B).
//
{
  static double Const = .398942280401433e0;
  static double brcomp,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
  static int i,n;
//
//  CONST = 1/SQRT(2*PI)
//
  static double T1,T2;

    brcomp = 0.0e0;
    if(*x == 0.0e0 || *y == 0.0e0) return brcomp;
    a0 = fifdmin1(*a,*b);
    if(a0 >= 8.0e0) goto S130;
    if(*x > 0.375e0) goto S10;
    lnx = log(*x);
    T1 = -*x;
    lny = alnrel(&T1);
    goto S30;
S10:
    if(*y > 0.375e0) goto S20;
    T2 = -*y;
    lnx = alnrel(&T2);
    lny = log(*y);
    goto S30;
S20:
    lnx = log(*x);
    lny = log(*y);
S30:
    z = *a*lnx+*b*lny;
    if(a0 < 1.0e0) goto S40;
    z -= beta_log(a,b);
    brcomp = exp(z);
    return brcomp;
S40:
//
//  PROCEDURE FOR A .LT. 1 OR B .LT. 1
//
    b0 = fifdmax1(*a,*b);
    if(b0 >= 8.0e0) goto S120;
    if(b0 > 1.0e0) goto S70;
//
//  ALGORITHM FOR B0 .LE. 1
//
    brcomp = exp(z);
    if(brcomp == 0.0e0) return brcomp;
    apb = *a+*b;
    if(apb > 1.0e0) goto S50;
    z = 1.0e0+gam1(&apb);
    goto S60;
S50:
    u = *a+*b-1.e0;
    z = (1.0e0+gam1(&u))/apb;
S60:
    c = (1.0e0+gam1(a))*(1.0e0+gam1(b))/z;
    brcomp = brcomp*(a0*c)/(1.0e0+a0/b0);
    return brcomp;
S70:
//
//  ALGORITHM FOR 1 .LT. B0 .LT. 8
//
    u = gamma_ln1 ( &a0 );
    n = ( int ) ( b0 - 1.0e0 );
    if(n < 1) goto S90;
    c = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b0 -= 1.0e0;
        c *= (b0/(a0+b0));
    }
    u = log(c)+u;
S90:
    z -= u;
    b0 -= 1.0e0;
    apb = a0+b0;
    if(apb > 1.0e0) goto S100;
    t = 1.0e0+gam1(&apb);
    goto S110;
S100:
    u = a0+b0-1.e0;
    t = (1.0e0+gam1(&u))/apb;
S110:
    brcomp = a0*exp(z)*(1.0e0+gam1(&b0))/t;
    return brcomp;
S120:
//
//  ALGORITHM FOR B0 .GE. 8
//
    u = gamma_ln1 ( &a0 ) + algdiv ( &a0, &b0 );
    brcomp = a0*exp(z-u);
    return brcomp;
S130:
//
//  PROCEDURE FOR A .GE. 8 AND B .GE. 8
//
    if(*a > *b) goto S140;
    h = *a/ *b;
    x0 = h/(1.0e0+h);
    y0 = 1.0e0/(1.0e0+h);
    lambda = *a-(*a+*b)**x;
    goto S150;
S140:
    h = *b/ *a;
    x0 = 1.0e0/(1.0e0+h);
    y0 = h/(1.0e0+h);
    lambda = (*a+*b)**y-*b;
S150:
    e = -(lambda/ *a);
    if(fabs(e) > 0.6e0) goto S160;
    u = rlog1(&e);
    goto S170;
S160:
    u = e-log(*x/x0);
S170:
    e = lambda/ *b;
    if(fabs(e) > 0.6e0) goto S180;
    v = rlog1(&e);
    goto S190;
S180:
    v = e-log(*y/y0);
S190:
    z = exp(-(*a*u+*b*v));
    brcomp = Const*sqrt(*b*x0)*z*exp(-bcorr(a,b));
    return brcomp;
}
//****************************************************************************80

double beta_rcomp1 ( int *mu, double *a, double *b, double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_RCOMP1 evaluates exp(MU) * X**A * Y**B / Beta(A,B).
//
//  Parameters:
//
//    Input, int MU, ?
//
//    Input, double A, B, the parameters of the Beta function.
//    A and B should be nonnegative.
//
//    Input, double X, Y, ?
//
//    Output, double BETA_RCOMP1, the value of
//    exp(MU) * X**A * Y**B / Beta(A,B).
//
{
  static double Const = .398942280401433e0;
  static double brcmp1,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
  static int i,n;
//
//     CONST = 1/SQRT(2*PI)
//
  static double T1,T2,T3,T4;

    a0 = fifdmin1(*a,*b);
    if(a0 >= 8.0e0) goto S130;
    if(*x > 0.375e0) goto S10;
    lnx = log(*x);
    T1 = -*x;
    lny = alnrel(&T1);
    goto S30;
S10:
    if(*y > 0.375e0) goto S20;
    T2 = -*y;
    lnx = alnrel(&T2);
    lny = log(*y);
    goto S30;
S20:
    lnx = log(*x);
    lny = log(*y);
S30:
    z = *a*lnx+*b*lny;
    if(a0 < 1.0e0) goto S40;
    z -= beta_log(a,b);
    brcmp1 = esum(mu,&z);
    return brcmp1;
S40:
//
//   PROCEDURE FOR A .LT. 1 OR B .LT. 1
//
    b0 = fifdmax1(*a,*b);
    if(b0 >= 8.0e0) goto S120;
    if(b0 > 1.0e0) goto S70;
//
//  ALGORITHM FOR B0 .LE. 1
//
    brcmp1 = esum(mu,&z);
    if(brcmp1 == 0.0e0) return brcmp1;
    apb = *a+*b;
    if(apb > 1.0e0) goto S50;
    z = 1.0e0+gam1(&apb);
    goto S60;
S50:
    u = *a+*b-1.e0;
    z = (1.0e0+gam1(&u))/apb;
S60:
    c = (1.0e0+gam1(a))*(1.0e0+gam1(b))/z;
    brcmp1 = brcmp1*(a0*c)/(1.0e0+a0/b0);
    return brcmp1;
S70:
//
//  ALGORITHM FOR 1 .LT. B0 .LT. 8
//
    u = gamma_ln1 ( &a0 );
    n = ( int ) ( b0 - 1.0e0 );
    if(n < 1) goto S90;
    c = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b0 -= 1.0e0;
        c *= (b0/(a0+b0));
    }
    u = log(c)+u;
S90:
    z -= u;
    b0 -= 1.0e0;
    apb = a0+b0;
    if(apb > 1.0e0) goto S100;
    t = 1.0e0+gam1(&apb);
    goto S110;
S100:
    u = a0+b0-1.e0;
    t = (1.0e0+gam1(&u))/apb;
S110:
    brcmp1 = a0*esum(mu,&z)*(1.0e0+gam1(&b0))/t;
    return brcmp1;
S120:
//
//  ALGORITHM FOR B0 .GE. 8
//
    u = gamma_ln1 ( &a0 ) + algdiv ( &a0, &b0 );
    T3 = z-u;
    brcmp1 = a0*esum(mu,&T3);
    return brcmp1;
S130:
//
//    PROCEDURE FOR A .GE. 8 AND B .GE. 8
//
    if(*a > *b) goto S140;
    h = *a/ *b;
    x0 = h/(1.0e0+h);
    y0 = 1.0e0/(1.0e0+h);
    lambda = *a-(*a+*b)**x;
    goto S150;
S140:
    h = *b/ *a;
    x0 = 1.0e0/(1.0e0+h);
    y0 = h/(1.0e0+h);
    lambda = (*a+*b)**y-*b;
S150:
    e = -(lambda/ *a);
    if(fabs(e) > 0.6e0) goto S160;
    u = rlog1(&e);
    goto S170;
S160:
    u = e-log(*x/x0);
S170:
    e = lambda/ *b;
    if(fabs(e) > 0.6e0) goto S180;
    v = rlog1(&e);
    goto S190;
S180:
    v = e-log(*y/y0);
S190:
    T4 = -(*a*u+*b*v);
    z = esum(mu,&T4);
    brcmp1 = Const*sqrt(*b*x0)*z*exp(-bcorr(a,b));
    return brcmp1;
}
//****************************************************************************80

double beta_up ( double *a, double *b, double *x, double *y, int *n,double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.
//
//    Input, double *X, *Y, ?
//
//    Input, int *N, ?
//
//    Input, double *EPS, the tolerance.
//
//    Output, double BETA_UP, the value of IX(A,B) - IX(A+N,B).
//
{
  static int K1 = 1;
  static int K2 = 0;
  static double bup,ap1,apb,d,l,r,t,w;
  static int i,k,kp1,mu,nm1;
//
//  OBTAIN THE SCALING FACTOR EXP(-MU) AND
//  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
//
    apb = *a+*b;
    ap1 = *a+1.0e0;
    mu = 0;
    d = 1.0e0;
    if(*n == 1 || *a < 1.0e0) goto S10;
    if(apb < 1.1e0*ap1) goto S10;
    mu = ( int ) fabs ( exparg(&K1) );
    k = ( int ) exparg ( &K2 );
    if(k < mu) mu = k;
    t = mu;
    d = exp(-t);
S10:
    bup = beta_rcomp1 ( &mu, a, b, x, y ) / *a;
    if(*n == 1 || bup == 0.0e0) return bup;
    nm1 = *n-1;
    w = d;
//
//  LET K BE THE INDEX OF THE MAXIMUM TERM
//
    k = 0;
    if(*b <= 1.0e0) goto S50;
    if(*y > 1.e-4) goto S20;
    k = nm1;
    goto S30;
S20:
    r = (*b-1.0e0)**x/ *y-*a;
    if(r < 1.0e0) goto S50;
    t = ( double ) nm1;
    k = nm1;
    if ( r < t ) k = ( int ) r;
S30:
//
//          ADD THE INCREASING TERMS OF THE SERIES
//
    for ( i = 1; i <= k; i++ )
    {
        l = i-1;
        d = (apb+l)/(ap1+l)**x*d;
        w = w + d;
    }
    if(k == nm1) goto S70;
S50:
//
//          ADD THE REMAINING TERMS OF THE SERIES
//
    kp1 = k+1;
    for ( i = kp1; i <= nm1; i++ )
    {
        l = i-1;
        d = (apb+l)/(ap1+l)**x*d;
        w = w + d;
        if(d <= *eps*w) goto S70;
    }
S70:
//
//  TERMINATE THE PROCEDURE
//
    bup *= w;
    return bup;
}
//****************************************************************************80

void binomial_cdf_values ( int *n_data, int *a, double *b, int *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *A, double *B, the parameters of the function.
//
//    Output, int *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 17

  int a_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  4,  4,  4,
     4, 10, 10, 10,
    10, 10, 10, 10,
    10 };
  double b_vec[N_MAX] = {
    0.05E+00, 0.05E+00, 0.05E+00, 0.50E+00,
    0.50E+00, 0.25E+00, 0.25E+00, 0.25E+00,
    0.25E+00, 0.05E+00, 0.10E+00, 0.15E+00,
    0.20E+00, 0.25E+00, 0.30E+00, 0.40E+00,
    0.50E+00 };
  double fx_vec[N_MAX] = {
    0.9025E+00, 0.9975E+00, 1.0000E+00, 0.2500E+00,
    0.7500E+00, 0.3164E+00, 0.7383E+00, 0.9492E+00,
    0.9961E+00, 0.9999E+00, 0.9984E+00, 0.9901E+00,
    0.9672E+00, 0.9219E+00, 0.8497E+00, 0.6331E+00,
    0.3770E+00 };
  int x_vec[N_MAX] = {
     0, 1, 2, 0,
     1, 0, 1, 2,
     3, 4, 4, 4,
     4, 4, 4, 4,
     4 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0.0E+00;
    *x = 0;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

void cdfbet ( int *which, double *p, double *q, double *x, double *y,double *a, double *b, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFBET evaluates the CDF of the Beta Distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the beta distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly by code associated with the reference.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The beta density is proportional to t^(A-1) * (1-t)^(B-1).
//
//  Modified:
//
//    09 June 2004
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which of the next four argument
//    values is to be calculated from the others.
//    1: Calculate P and Q from X, Y, A and B;
//    2: Calculate X and Y from P, Q, A and B;
//    3: Calculate A from P, Q, X, Y and B;
//    4: Calculate B from P, Q, X, Y and A.
//
//    Input/output, double *P, the integral from 0 to X of the
//    chi-square distribution.  Input range: [0, 1].
//
//    Input/output, double *Q, equals 1-P.  Input range: [0, 1].
//
//    Input/output, double *X, the upper limit of integration
//    of the beta density.  If it is an input value, it should lie in
//    the range [0,1].  If it is an output value, it will be searched for
//    in the range [0,1].
//
//    Input/output, double *Y, equal to 1-X.  If it is an input
//    value, it should lie in the range [0,1].  If it is an output value,
//    it will be searched for in the range [0,1].
//
//    Input/output, double *A, the first parameter of the beta
//    density.  If it is an input value, it should lie in the range
//    (0, +infinity).  If it is an output value, it will be searched
//    for in the range [1D-300,1D300].
//
//    Input/output, double *B, the second parameter of the beta
//    density.  If it is an input value, it should lie in the range
//    (0, +infinity).  If it is an output value, it will be searched
//    for in the range [1D-300,1D300].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if X + Y /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300
# define one 1.0e0

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K3 = 1.0e0;
  static double K8 = 0.5e0;
  static double K9 = 5.0e0;
  static double fx,xhi,xlo,cum,ccum,xy,pq;
  static unsigned long qhi,qleft,qporq;
  static double T4,T5,T6,T7,T10,T11,T12,T13,T14,T15;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q < 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S150;
//
//     X
//
    if(!(*x < 0.0e0 || *x > 1.0e0)) goto S140;
    if(!(*x < 0.0e0)) goto S120;
    *bound = 0.0e0;
    goto S130;
S120:
    *bound = 1.0e0;
S130:
    *status = -4;
    return;
S150:
S140:
    if(*which == 2) goto S190;
//
//     Y
//
    if(!(*y < 0.0e0 || *y > 1.0e0)) goto S180;
    if(!(*y < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = -5;
    return;
S190:
S180:
    if(*which == 3) goto S210;
//
//     A
//
    if(!(*a <= 0.0e0)) goto S200;
    *bound = 0.0e0;
    *status = -6;
    return;
S210:
S200:
    if(*which == 4) goto S230;
//
//     B
//
    if(!(*b <= 0.0e0)) goto S220;
    *bound = 0.0e0;
    *status = -7;
    return;
S230:
S220:
    if(*which == 1) goto S270;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S260;
    if(!(pq < 0.0e0)) goto S240;
    *bound = 0.0e0;
    goto S250;
S240:
    *bound = 1.0e0;
S250:
    *status = 3;
    return;
S270:
S260:
    if(*which == 2) goto S310;
//
//     X + Y
//
    xy = *x+*y;
    if(!(fabs(xy-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S300;
    if(!(xy < 0.0e0)) goto S280;
    *bound = 0.0e0;
    goto S290;
S280:
    *bound = 1.0e0;
S290:
    *status = 4;
    return;
S310:
S300:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P and Q
//
        cumbet(x,y,a,b,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating X and Y
//
        T4 = atol;
        T5 = tol;
        dstzr(&K2,&K3,&T4,&T5);
        if(!qporq) goto S340;
        *status = 0;
        dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
        *y = one-*x;
S320:
        if(!(*status == 1)) goto S330;
        cumbet(x,y,a,b,&cum,&ccum);
        fx = cum-*p;
        dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
        *y = one-*x;
        goto S320;
S330:
        goto S370;
S340:
        *status = 0;
        dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
        *x = one-*y;
S350:
        if(!(*status == 1)) goto S360;
        cumbet(x,y,a,b,&cum,&ccum);
        fx = ccum-*q;
        dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
        *x = one-*y;
        goto S350;
S370:
S360:
        if(!(*status == -1)) goto S400;
        if(!qleft) goto S380;
        *status = 1;
        *bound = 0.0e0;
        goto S390;
S380:
        *status = 2;
        *bound = 1.0e0;
S400:
S390:
        ;
    }
    else if(3 == *which) {
//
//     Computing A
//
        *a = 5.0e0;
        T6 = zero;
        T7 = inf;
        T10 = atol;
        T11 = tol;
        dstinv(&T6,&T7,&K8,&K8,&K9,&T10,&T11);
        *status = 0;
        dinvr(status,a,&fx,&qleft,&qhi);
S410:
        if(!(*status == 1)) goto S440;
        cumbet(x,y,a,b,&cum,&ccum);
        if(!qporq) goto S420;
        fx = cum-*p;
        goto S430;
S420:
        fx = ccum-*q;
S430:
        dinvr(status,a,&fx,&qleft,&qhi);
        goto S410;
S440:
        if(!(*status == -1)) goto S470;
        if(!qleft) goto S450;
        *status = 1;
        *bound = zero;
        goto S460;
S450:
        *status = 2;
        *bound = inf;
S470:
S460:
        ;
    }
    else if(4 == *which) {
//
//     Computing B
//
        *b = 5.0e0;
        T12 = zero;
        T13 = inf;
        T14 = atol;
        T15 = tol;
        dstinv(&T12,&T13,&K8,&K8,&K9,&T14,&T15);
        *status = 0;
        dinvr(status,b,&fx,&qleft,&qhi);
S480:
        if(!(*status == 1)) goto S510;
        cumbet(x,y,a,b,&cum,&ccum);
        if(!qporq) goto S490;
        fx = cum-*p;
        goto S500;
S490:
        fx = ccum-*q;
S500:
        dinvr(status,b,&fx,&qleft,&qhi);
        goto S480;
S510:
        if(!(*status == -1)) goto S540;
        if(!qleft) goto S520;
        *status = 1;
        *bound = zero;
        goto S530;
S520:
        *status = 2;
        *bound = inf;
S530:
        ;
    }
S540:
    return;
# undef tol
# undef atol
# undef zero
# undef inf
# undef one
}
//****************************************************************************80

void cdfbin ( int *which, double *p, double *q, double *s, double *xn, double *pr, double *ompr, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFBIN evaluates the CDF of the Binomial distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the binomial distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    P is the probablility of S or fewer successes in XN binomial trials,
//    each trial having an individual probability of success of PR.
//
//  Modified:
//
//    09 June 2004
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.24.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which of argument values is to
//    be calculated from the others.
//    1: Calculate P and Q from S, XN, PR and OMPR;
//    2: Calculate S from P, Q, XN, PR and OMPR;
//    3: Calculate XN from P, Q, S, PR and OMPR;
//    4: Calculate PR and OMPR from P, Q, S and XN.
//
//    Input/output, double *P, the cumulation, from 0 to S,
//    of the binomial distribution.  If P is an input value, it should
//    lie in the range [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *S, the number of successes observed.
//    Whether this is an input or output value, it should lie in the
//    range [0,XN].
//
//    Input/output, double *XN, the number of binomial trials.
//    If this is an input value it should lie in the range: (0, +infinity).
//    If it is an output value it will be searched for in the
//    range [1.0D-300, 1.0D+300].
//
//    Input/output, double *PR, the probability of success in each
//    binomial trial.  Whether this is an input or output value, it should
//    lie in the range: [0,1].
//
//    Input/output, double *OMPR, equal to 1-PR.  Whether this is an
//    input or output value, it should lie in the range [0,1].  Also, it should
//    be the case that PR + OMPR = 1.
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if PR + OMPR /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define atol (1.0e-50)
# define tol (1.0e-8)
# define zero (1.0e-300)
# define inf 1.0e300
# define one 1.0e0

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K3 = 0.5e0;
  static double K4 = 5.0e0;
  static double K11 = 1.0e0;
  static double fx,xhi,xlo,cum,ccum,pq,prompr;
  static unsigned long qhi,qleft,qporq;
  static double T5,T6,T7,T8,T9,T10,T12,T13;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 && *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q < 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 3) goto S130;
//
//     XN
//
    if(!(*xn <= 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -5;
    return;
S130:
S120:
    if(*which == 2) goto S170;
//
//     S
//
    if(!(*s < 0.0e0 || *which != 3 && *s > *xn)) goto S160;
    if(!(*s < 0.0e0)) goto S140;
    *bound = 0.0e0;
    goto S150;
S140:
    *bound = *xn;
S150:
    *status = -4;
    return;
S170:
S160:
    if(*which == 4) goto S210;
//
//     PR
//
    if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S200;
    if(!(*pr < 0.0e0)) goto S180;
    *bound = 0.0e0;
    goto S190;
S180:
    *bound = 1.0e0;
S190:
    *status = -6;
    return;
S210:
S200:
    if(*which == 4) goto S250;
//
//     OMPR
//
    if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S240;
    if(!(*ompr < 0.0e0)) goto S220;
    *bound = 0.0e0;
    goto S230;
S220:
    *bound = 1.0e0;
S230:
    *status = -7;
    return;
S250:
S240:
    if(*which == 1) goto S290;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S280;
    if(!(pq < 0.0e0)) goto S260;
    *bound = 0.0e0;
    goto S270;
S260:
    *bound = 1.0e0;
S270:
    *status = 3;
    return;
S290:
S280:
    if(*which == 4) goto S330;
//
//     PR + OMPR
//
    prompr = *pr+*ompr;
    if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S320;
    if(!(prompr < 0.0e0)) goto S300;
    *bound = 0.0e0;
    goto S310;
S300:
    *bound = 1.0e0;
S310:
    *status = 4;
    return;
S330:
S320:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cumbin(s,xn,pr,ompr,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating S
//
        *s = 5.0e0;
        T5 = atol;
        T6 = tol;
        dstinv(&K2,xn,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        dinvr(status,s,&fx,&qleft,&qhi);
S340:
        if(!(*status == 1)) goto S370;
        cumbin(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S350;
        fx = cum-*p;
        goto S360;
S350:
        fx = ccum-*q;
S360:
        dinvr(status,s,&fx,&qleft,&qhi);
        goto S340;
S370:
        if(!(*status == -1)) goto S400;
        if(!qleft) goto S380;
        *status = 1;
        *bound = 0.0e0;
        goto S390;
S380:
        *status = 2;
        *bound = *xn;
S400:
S390:
        ;
    }
    else if(3 == *which) {
//
//     Calculating XN
//
        *xn = 5.0e0;
        T7 = zero;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        dinvr(status,xn,&fx,&qleft,&qhi);
S410:
        if(!(*status == 1)) goto S440;
        cumbin(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S420;
        fx = cum-*p;
        goto S430;
S420:
        fx = ccum-*q;
S430:
        dinvr(status,xn,&fx,&qleft,&qhi);
        goto S410;
S440:
        if(!(*status == -1)) goto S470;
        if(!qleft) goto S450;
        *status = 1;
        *bound = zero;
        goto S460;
S450:
        *status = 2;
        *bound = inf;
S470:
S460:
        ;
    }
    else if(4 == *which) {
//
//     Calculating PR and OMPR
//
        T12 = atol;
        T13 = tol;
        dstzr(&K2,&K11,&T12,&T13);
        if(!qporq) goto S500;
        *status = 0;
        dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
S480:
        if(!(*status == 1)) goto S490;
        cumbin(s,xn,pr,ompr,&cum,&ccum);
        fx = cum-*p;
        dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
        goto S480;
S490:
        goto S530;
S500:
        *status = 0;
        dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
S510:
        if(!(*status == 1)) goto S520;
        cumbin(s,xn,pr,ompr,&cum,&ccum);
        fx = ccum-*q;
        dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
        goto S510;
S530:
S520:
        if(!(*status == -1)) goto S560;
        if(!qleft) goto S540;
        *status = 1;
        *bound = 0.0e0;
        goto S550;
S540:
        *status = 2;
        *bound = 1.0e0;
S550:
        ;
    }
S560:
    return;
# undef atol
# undef tol
# undef zero
# undef inf
# undef one
}
//****************************************************************************80

void cdfchi ( int *which, double *p, double *q, double *x, double *df,int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFCHI evaluates the CDF of the chi square distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the chi square distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The CDF of the chi square distribution can be evaluated
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ ChiSquareDistribution [ DF ], X ]
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.4.19.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X and DF;
//    2: Calculate X from P, Q and DF;
//    3: Calculate DF from P, Q and X.
//
//    Input/output, double *P, the integral from 0 to X of
//    the chi-square distribution.  If this is an input value, it should
//    lie in the range [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *X, the upper limit of integration
//    of the chi-square distribution.  If this is an input
//    value, it should lie in the range: [0, +infinity).  If it is an output
//    value, it will be searched for in the range: [0,1.0D+300].
//
//    Input/output, double *DF, the degrees of freedom of the
//    chi-square distribution.  If this is an input value, it should lie
//    in the range: (0, +infinity).  If it is an output value, it will be
//    searched for in the range: [ 1.0D-300, 1.0D+300].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +10, an error was returned from CUMGAM.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double fx,cum,ccum,pq,porq;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10,T11;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 3.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     X
//
    if(!(*x < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     DF
//
    if(!(*df <= 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 1) goto S190;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S180;
    if(!(pq < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = 3;
    return;
S190:
S180:
    if(*which == 1) goto S220;
//
//     Select the minimum of P or Q
//
    qporq = *p <= *q;
    if(!qporq) goto S200;
    porq = *p;
    goto S210;
S200:
    porq = *q;
S220:
S210:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P and Q
//
        *status = 0;
        cumchi(x,df,p,q);
        if(porq > 1.5e0) {
            *status = 10;
            return;
        }
    }
    else if(2 == *which) {
//
//     Calculating X
//
        *x = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        dinvr(status,x,&fx,&qleft,&qhi);
S230:
        if(!(*status == 1)) goto S270;
        cumchi(x,df,&cum,&ccum);
        if(!qporq) goto S240;
        fx = cum-*p;
        goto S250;
S240:
        fx = ccum-*q;
S250:
        if(!(fx+porq > 1.5e0)) goto S260;
        *status = 10;
        return;
S260:
        dinvr(status,x,&fx,&qleft,&qhi);
        goto S230;
S270:
        if(!(*status == -1)) goto S300;
        if(!qleft) goto S280;
        *status = 1;
        *bound = 0.0e0;
        goto S290;
S280:
        *status = 2;
        *bound = inf;
S300:
S290:
        ;
    }
    else if(3 == *which) {
//
//  Calculating DF
//
        *df = 5.0e0;
        T8 = zero;
        T9 = inf;
        T10 = atol;
        T11 = tol;
        dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
        *status = 0;
        dinvr(status,df,&fx,&qleft,&qhi);
S310:
        if(!(*status == 1)) goto S350;
        cumchi(x,df,&cum,&ccum);
        if(!qporq) goto S320;
        fx = cum-*p;
        goto S330;
S320:
        fx = ccum-*q;
S330:
        if(!(fx+porq > 1.5e0)) goto S340;
        *status = 10;
        return;
S340:
        dinvr(status,df,&fx,&qleft,&qhi);
        goto S310;
S350:
        if(!(*status == -1)) goto S380;
        if(!qleft) goto S360;
        *status = 1;
        *bound = zero;
        goto S370;
S360:
        *status = 2;
        *bound = inf;
S370:
        ;
    }
S380:
    return;
# undef tol
# undef atol
# undef zero
# undef inf
}
//****************************************************************************80

void cdfchn ( int *which, double *p, double *q, double *x, double *df,double *pnonc, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFCHN evaluates the CDF of the Noncentral Chi-Square.
//
//  Discussion:
//
//    This routine calculates any one parameter of the noncentral chi-square
//    distribution given values for the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The computation time required for this routine is proportional
//    to the noncentrality parameter (PNONC).  Very large values of
//    this parameter can consume immense computer resources.  This is
//    why the search range is bounded by 10,000.
//
//    The CDF of the noncentral chi square distribution can be evaluated
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.25.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, DF and PNONC;
//    2: Calculate X from P, DF and PNONC;
//    3: Calculate DF from P, X and PNONC;
//    4: Calculate PNONC from P, X and DF.
//
//    Input/output, double *P, the integral from 0 to X of
//    the noncentral chi-square distribution.  If this is an input
//    value, it should lie in the range: [0, 1.0-1.0D-16).
//
//    Input/output, double *Q, is generally not used by this
//    subroutine and is only included for similarity with other routines.
//    However, if P is to be computed, then a value will also be computed
//    for Q.
//
//    Input, double *X, the upper limit of integration of the
//    noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be sought in the range: [0,1.0D+300].
//
//    Input/output, double *DF, the number of degrees of freedom
//    of the noncentral chi-square distribution.  If this is an input value,
//    it should lie in the range: (0, +infinity).  If it is an output value,
//    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
//
//    Input/output, double *PNONC, the noncentrality parameter of
//    the noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be searched for in the range: [0,1.0D+4]
//
//    Output, int *STATUS, reports on the calculation.
//    0, if calculation completed correctly;
//    -I, if input parameter number I is out of range;
//    1, if the answer appears to be lower than the lowest search bound;
//    2, if the answer appears to be higher than the greatest search bound.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tent4 1.0e4
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define one (1.0e0-1.0e-16)
# define inf 1.0e300

  static double K1 = 0.0e0;
  static double K3 = 0.5e0;
  static double K4 = 5.0e0;
  static double fx,cum,ccum;
  static unsigned long qhi,qleft;
  static double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > one)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = one;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 2) goto S90;
//
//     X
//
    if(!(*x < 0.0e0)) goto S80;
    *bound = 0.0e0;
    *status = -4;
    return;
S90:
S80:
    if(*which == 3) goto S110;
//
//     DF
//
    if(!(*df <= 0.0e0)) goto S100;
    *bound = 0.0e0;
    *status = -5;
    return;
S110:
S100:
    if(*which == 4) goto S130;
//
//     PNONC
//
    if(!(*pnonc < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -6;
    return;
S130:
S120:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P and Q
//
        cumchn(x,df,pnonc,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating X
//
        *x = 5.0e0;
        T2 = inf;
        T5 = atol;
        T6 = tol;
        dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        dinvr(status,x,&fx,&qleft,&qhi);
S140:
        if(!(*status == 1)) goto S150;
        cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,x,&fx,&qleft,&qhi);
        goto S140;
S150:
        if(!(*status == -1)) goto S180;
        if(!qleft) goto S160;
        *status = 1;
        *bound = 0.0e0;
        goto S170;
S160:
        *status = 2;
        *bound = inf;
S180:
S170:
        ;
    }
    else if(3 == *which) {
//
//     Calculating DF
//
        *df = 5.0e0;
        T7 = zero;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        dinvr(status,df,&fx,&qleft,&qhi);
S190:
        if(!(*status == 1)) goto S200;
        cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,df,&fx,&qleft,&qhi);
        goto S190;
S200:
        if(!(*status == -1)) goto S230;
        if(!qleft) goto S210;
        *status = 1;
        *bound = zero;
        goto S220;
S210:
        *status = 2;
        *bound = inf;
S230:
S220:
        ;
    }
    else if(4 == *which) {
//
//     Calculating PNONC
//
        *pnonc = 5.0e0;
        T11 = tent4;
        T12 = atol;
        T13 = tol;
        dstinv(&K1,&T11,&K3,&K3,&K4,&T12,&T13);
        *status = 0;
        dinvr(status,pnonc,&fx,&qleft,&qhi);
S240:
        if(!(*status == 1)) goto S250;
        cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,pnonc,&fx,&qleft,&qhi);
        goto S240;
S250:
        if(!(*status == -1)) goto S280;
        if(!qleft) goto S260;
        *status = 1;
        *bound = zero;
        goto S270;
S260:
        *status = 2;
        *bound = tent4;
S270:
        ;
    }
S280:
    return;
# undef tent4
# undef tol
# undef atol
# undef zero
# undef one
# undef inf
}
//****************************************************************************80

void cdff ( int *which, double *p, double *q, double *f, double *dfn,double *dfd, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFF evaluates the CDF of the F distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the F distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The value of the cumulative F distribution is not necessarily
//    monotone in either degree of freedom.  There thus may be two
//    values that provide a given CDF value.  This routine assumes
//    monotonicity and will find an arbitrary one of the two values.
//
//  Modified:
//
//    14 April 2007
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.6.2.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, DFN and DFD;
//    2: Calculate F from P, Q, DFN and DFD;
//    3: Calculate DFN from P, Q, F and DFD;
//    4: Calculate DFD from P, Q, F and DFN.
//
//    Input/output, double *P, the integral from 0 to F of
//    the F-density.  If it is an input value, it should lie in the
//    range [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *F, the upper limit of integration
//    of the F-density.  If this is an input value, it should lie in the
//    range [0, +infinity).  If it is an output value, it will be searched
//    for in the range [0,1.0D+300].
//
//    Input/output, double *DFN, the number of degrees of
//    freedom of the numerator sum of squares.  If this is an input value,
//    it should lie in the range: (0, +infinity).  If it is an output value,
//    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
//
//    Input/output, double *DFD, the number of degrees of freedom
//    of the denominator sum of squares.  If this is an input value, it should
//    lie in the range: (0, +infinity).  If it is an output value, it will
//    be searched for in the  range: [ 1.0D-300, 1.0D+300].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double pq,fx,cum,ccum;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15;

  *status = 0;
  *bound = 0.0;
//
//  Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     F
//
    if(!(*f < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     DFN
//
    if(!(*dfn <= 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 4) goto S170;
//
//     DFD
//
    if(!(*dfd <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
    if(*which == 1) goto S210;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S200;
    if(!(pq < 0.0e0)) goto S180;
    *bound = 0.0e0;
    goto S190;
S180:
    *bound = 1.0e0;
S190:
    *status = 3;
    return;
S210:
S200:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cumf(f,dfn,dfd,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating F
//
        *f = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        dinvr(status,f,&fx,&qleft,&qhi);
S220:
        if(!(*status == 1)) goto S250;
        cumf(f,dfn,dfd,&cum,&ccum);
        if(!qporq) goto S230;
        fx = cum-*p;
        goto S240;
S230:
        fx = ccum-*q;
S240:
        dinvr(status,f,&fx,&qleft,&qhi);
        goto S220;
S250:
        if(!(*status == -1)) goto S280;
        if(!qleft) goto S260;
        *status = 1;
        *bound = 0.0e0;
        goto S270;
S260:
        *status = 2;
        *bound = inf;
S280:
S270:
        ;
    }
//
//  Calculate DFN.
//
//  Note that, in the original calculation, the lower bound for DFN was 0.
//  Using DFN = 0 causes an error in CUMF when it calls BETA_INC.
//  The lower bound was set to the more reasonable value of 1.
//  JVB, 14 April 2007.
//
  else if ( 3 == *which )
  {

    T8 = 1.0;
    T9 = inf;
    T10 = atol;
    T11 = tol;
    dstinv ( &T8, &T9, &K4, &K4, &K5, &T10, &T11 );

    *status = 0;
    *dfn = 5.0;
    fx = 0.0;

    dinvr ( status, dfn, &fx, &qleft, &qhi );

    while ( *status == 1 )
    {
      cumf ( f, dfn, dfd, &cum, &ccum );

      if ( *p <= *q )
      {
        fx = cum - *p;
      }
      else
      {
        fx = ccum - *q;
      }
      dinvr ( status, dfn, &fx, &qleft, &qhi );
    }

    if ( *status == -1 )
    {
      if ( qleft )
      {
        *status = 1;
        *bound = 1.0;
      }
      else
      {
        *status = 2;
        *bound = inf;
      }
    }
  }
//
//  Calculate DFD.
//
//  Note that, in the original calculation, the lower bound for DFD was 0.
//  Using DFD = 0 causes an error in CUMF when it calls BETA_INC.
//  The lower bound was set to the more reasonable value of 1.
//  JVB, 14 April 2007.
//
//
  else if ( 4 == *which )
  {

    T12 = 1.0;
    T13 = inf;
    T14 = atol;
    T15 = tol;
    dstinv ( &T12, &T13, &K4, &K4, &K5, &T14, &T15 );

    *status = 0;
    *dfd = 5.0;
    fx = 0.0;
    dinvr ( status, dfd, &fx, &qleft, &qhi );

    while ( *status == 1 )
    {
      cumf ( f, dfn, dfd, &cum, &ccum );

      if ( *p <= *q )
      {
        fx = cum - *p;
      }
      else
      {
        fx = ccum - *q;
      }
      dinvr ( status, dfd, &fx, &qleft, &qhi );
    }

    if ( *status == -1 )
    {
      if ( qleft )
      {
        *status = 1;
        *bound = 1.0;
      }
      else
      {
        *status = 2;
        *bound = inf;
      }
    }
  }

  return;
# undef tol
# undef atol
# undef zero
# undef inf
}
//****************************************************************************80

void cdffnc ( int *which, double *p, double *q, double *f, double *dfn,double *dfd, double *phonc, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFFNC evaluates the CDF of the Noncentral F distribution.
//
//  Discussion:
//
//    This routine originally used 1.0E+300 as the upper bound for the
//    interval in which many of the missing parameters are to be sought.
//    Since the underlying rootfinder routine needs to evaluate the
//    function at this point, it is no surprise that the program was
//    experiencing overflows.  A less extravagant upper bound
//    is being tried for now!
//
//
//    This routine calculates any one parameter of the Noncentral F distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The computation time required for this routine is proportional
//    to the noncentrality parameter PNONC.  Very large values of
//    this parameter can consume immense computer resources.  This is
//    why the search range is bounded by 10,000.
//
//    The value of the cumulative noncentral F distribution is not
//    necessarily monotone in either degree of freedom.  There thus
//    may be two values that provide a given CDF value.  This routine
//    assumes monotonicity and will find an arbitrary one of the two
//    values.
//
//    The CDF of the noncentral F distribution can be evaluated
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ NoncentralFRatioDistribution [ DFN, DFD, PNONC ], X ]
//
//  Modified:
//
//    15 June 2004
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.6.20.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, DFN, DFD and PNONC;
//    2: Calculate F from P, Q, DFN, DFD and PNONC;
//    3: Calculate DFN from P, Q, F, DFD and PNONC;
//    4: Calculate DFD from P, Q, F, DFN and PNONC;
//    5: Calculate PNONC from P, Q, F, DFN and DFD.
//
//    Input/output, double *P, the integral from 0 to F of
//    the noncentral F-density.  If P is an input value it should
//    lie in the range [0,1) (Not including 1!).
//
//    Dummy, double *Q, is not used by this subroutine,
//    and is only included for similarity with the other routines.
//    Its input value is not checked.  If P is to be computed, the
//    Q is set to 1 - P.
//
//    Input/output, double *F, the upper limit of integration
//    of the noncentral F-density.  If this is an input value, it should
//    lie in the range: [0, +infinity).  If it is an output value, it
//    will be searched for in the range: [0,1.0D+30].
//
//    Input/output, double *DFN, the number of degrees of freedom
//    of the numerator sum of squares.  If this is an input value, it should
//    lie in the range: (0, +infinity).  If it is an output value, it will
//    be searched for in the range: [ 1.0, 1.0D+30].
//
//    Input/output, double *DFD, the number of degrees of freedom
//    of the denominator sum of squares.  If this is an input value, it should
//    be in range: (0, +infinity).  If it is an output value, it will be
//    searched for in the range [1.0, 1.0D+30].
//
//    Input/output, double *PNONC, the noncentrality parameter
//    If this is an input value, it should be nonnegative.
//    If it is an output value, it will be searched for in the range: [0,1.0D+4].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tent4 1.0e4
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define one (1.0e0-1.0e-16)
# define inf 1.0e300

  static double K1 = 0.0e0;
  static double K3 = 0.5e0;
  static double K4 = 5.0e0;
  static double fx,cum,ccum;
  static unsigned long qhi,qleft;
  static double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 5)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 5.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > one)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = one;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 2) goto S90;
//
//     F
//
    if(!(*f < 0.0e0)) goto S80;
    *bound = 0.0e0;
    *status = -4;
    return;
S90:
S80:
    if(*which == 3) goto S110;
//
//     DFN
//
    if(!(*dfn <= 0.0e0)) goto S100;
    *bound = 0.0e0;
    *status = -5;
    return;
S110:
S100:
    if(*which == 4) goto S130;
//
//     DFD
//
    if(!(*dfd <= 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -6;
    return;
S130:
S120:
    if(*which == 5) goto S150;
//
//     PHONC
//
    if(!(*phonc < 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -7;
    return;
S150:
S140:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cumfnc(f,dfn,dfd,phonc,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating F
//
        *f = 5.0e0;
        T2 = inf;
        T5 = atol;
        T6 = tol;
        dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        dinvr(status,f,&fx,&qleft,&qhi);
S160:
        if(!(*status == 1)) goto S170;
        cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,f,&fx,&qleft,&qhi);
        goto S160;
S170:
        if(!(*status == -1)) goto S200;
        if(!qleft) goto S180;
        *status = 1;
        *bound = 0.0e0;
        goto S190;
S180:
        *status = 2;
        *bound = inf;
S200:
S190:
        ;
    }
    else if(3 == *which) {
//
//     Calculating DFN
//
        *dfn = 5.0e0;
        T7 = zero;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        dinvr(status,dfn,&fx,&qleft,&qhi);
S210:
        if(!(*status == 1)) goto S220;
        cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,dfn,&fx,&qleft,&qhi);
        goto S210;
S220:
        if(!(*status == -1)) goto S250;
        if(!qleft) goto S230;
        *status = 1;
        *bound = zero;
        goto S240;
S230:
        *status = 2;
        *bound = inf;
S250:
S240:
        ;
    }
    else if(4 == *which) {
//
//     Calculating DFD
//
        *dfd = 5.0e0;
        T11 = zero;
        T12 = inf;
        T13 = atol;
        T14 = tol;
        dstinv(&T11,&T12,&K3,&K3,&K4,&T13,&T14);
        *status = 0;
        dinvr(status,dfd,&fx,&qleft,&qhi);
S260:
        if(!(*status == 1)) goto S270;
        cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,dfd,&fx,&qleft,&qhi);
        goto S260;
S270:
        if(!(*status == -1)) goto S300;
        if(!qleft) goto S280;
        *status = 1;
        *bound = zero;
        goto S290;
S280:
        *status = 2;
        *bound = inf;
S300:
S290:
        ;
    }
    else if(5 == *which) {
//
//     Calculating PHONC
//
        *phonc = 5.0e0;
        T15 = tent4;
        T16 = atol;
        T17 = tol;
        dstinv(&K1,&T15,&K3,&K3,&K4,&T16,&T17);
        *status = 0;
        dinvr(status,phonc,&fx,&qleft,&qhi);
S310:
        if(!(*status == 1)) goto S320;
        cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,phonc,&fx,&qleft,&qhi);
        goto S310;
S320:
        if(!(*status == -1)) goto S350;
        if(!qleft) goto S330;
        *status = 1;
        *bound = 0.0e0;
        goto S340;
S330:
        *status = 2;
        *bound = tent4;
S340:
        ;
    }
S350:
    return;
# undef tent4
# undef tol
# undef atol
# undef zero
# undef one
# undef inf
}
//****************************************************************************80

void cdfgam ( int *which, double *p, double *q, double *x, double *shape,
  double *scale, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFGAM evaluates the CDF of the Gamma Distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the Gamma distribution
//    given the others.
//
//    The cumulative distribution function P is calculated directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The gamma density is proportional to T**(SHAPE - 1) * EXP(- SCALE * T)
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Computation of the incomplete gamma function ratios and their inverse,
//    ACM Transactions on Mathematical Software,
//    Volume 12, 1986, pages 377-393.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, SHAPE and SCALE;
//    2: Calculate X from P, Q, SHAPE and SCALE;
//    3: Calculate SHAPE from P, Q, X and SCALE;
//    4: Calculate SCALE from P, Q, X and SHAPE.
//
//    Input/output, double *P, the integral from 0 to X of the
//    Gamma density.  If this is an input value, it should lie in the
//    range: [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *X, the upper limit of integration of
//    the Gamma density.  If this is an input value, it should lie in the
//    range: [0, +infinity).  If it is an output value, it will lie in
//    the range: [0,1E300].
//
//    Input/output, double *SHAPE, the shape parameter of the
//    Gamma density.  If this is an input value, it should lie in the range:
//    (0, +infinity).  If it is an output value, it will be searched for
//    in the range: [1.0D-300,1.0D+300].
//
//    Input/output, double *SCALE, the scale parameter of the
//    Gamma density.  If this is an input value, it should lie in the range
//    (0, +infinity).  If it is an output value, it will be searched for
//    in the range: (1.0D-300,1.0D+300].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +10, if the Gamma or inverse Gamma routine cannot compute the answer.
//    This usually happens only for X and SHAPE very large (more than 1.0D+10.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e300

  static int K1 = 1;
  static double K5 = 0.5e0;
  static double K6 = 5.0e0;
  static double xx,fx,xscale,cum,ccum,pq,porq;
  static int ierr;
  static unsigned long qhi,qleft,qporq;
  static double T2,T3,T4,T7,T8,T9;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     X
//
    if(!(*x < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     SHAPE
//
    if(!(*shape <= 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 4) goto S170;
//
//     SCALE
//
    if(!(*scale <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
    if(*which == 1) goto S210;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S200;
    if(!(pq < 0.0e0)) goto S180;
    *bound = 0.0e0;
    goto S190;
S180:
    *bound = 1.0e0;
S190:
    *status = 3;
    return;
S210:
S200:
    if(*which == 1) goto S240;
//
//     Select the minimum of P or Q
//
    qporq = *p <= *q;
    if(!qporq) goto S220;
    porq = *p;
    goto S230;
S220:
    porq = *q;
S240:
S230:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        *status = 0;
        xscale = *x**scale;
        cumgam(&xscale,shape,p,q);
        if(porq > 1.5e0) *status = 10;
    }
    else if(2 == *which) {
//
//     Computing X
//
        T2 = -1.0e0;
        gamma_inc_inv ( shape, &xx, &T2, p, q, &ierr );
        if(ierr < 0.0e0) {
            *status = 10;
            return;
        }
        else  {
            *x = xx/ *scale;
            *status = 0;
        }
    }
    else if(3 == *which) {
//
//     Computing SHAPE
//
        *shape = 5.0e0;
        xscale = *x**scale;
        T3 = zero;
        T4 = inf;
        T7 = atol;
        T8 = tol;
        dstinv(&T3,&T4,&K5,&K5,&K6,&T7,&T8);
        *status = 0;
        dinvr(status,shape,&fx,&qleft,&qhi);
S250:
        if(!(*status == 1)) goto S290;
        cumgam(&xscale,shape,&cum,&ccum);
        if(!qporq) goto S260;
        fx = cum-*p;
        goto S270;
S260:
        fx = ccum-*q;
S270:
        if(!(qporq && cum > 1.5e0 || !qporq && ccum > 1.5e0)) goto S280;
        *status = 10;
        return;
S280:
        dinvr(status,shape,&fx,&qleft,&qhi);
        goto S250;
S290:
        if(!(*status == -1)) goto S320;
        if(!qleft) goto S300;
        *status = 1;
        *bound = zero;
        goto S310;
S300:
        *status = 2;
        *bound = inf;
S320:
S310:
        ;
    }
    else if(4 == *which) {
//
//     Computing SCALE
//
        T9 = -1.0e0;
        gamma_inc_inv ( shape, &xx, &T9, p, q, &ierr );
        if(ierr < 0.0e0) {
            *status = 10;
            return;
        }
        else  {
            *scale = xx/ *x;
            *status = 0;
        }
    }
    return;
# undef tol
# undef atol
# undef zero
# undef inf
}
//****************************************************************************80

void cdfnbn ( int *which, double *p, double *q, double *s, double *xn,
  double *pr, double *ompr, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFNBN evaluates the CDF of the Negative Binomial distribution
//
//  Discussion:
//
//    This routine calculates any one parameter of the negative binomial
//    distribution given values for the others.
//
//    The cumulative negative binomial distribution returns the
//    probability that there will be F or fewer failures before the
//    S-th success in binomial trials each of which has probability of
//    success PR.
//
//    The individual term of the negative binomial is the probability of
//    F failures before S successes and is
//    Choose( F, S+F-1 ) * PR^(S) * (1-PR)^F
//
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.26.
//
//  Parameters:
//
//    Input, int WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, S, PR and OMPR;
//    2: Calculate F from P, Q, S, PR and OMPR;
//    3: Calculate S from P, Q, F, PR and OMPR;
//    4: Calculate PR and OMPR from P, Q, F and S.
//
//    Input/output, double P, the cumulation from 0 to F of
//    the negative binomial distribution.  If P is an input value, it
//    should lie in the range [0,1].
//
//    Input/output, double Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double F, the upper limit of cumulation of
//    the binomial distribution.  There are F or fewer failures before
//    the S-th success.  If this is an input value, it may lie in the
//    range [0,+infinity), and if it is an output value, it will be searched
//    for in the range [0,1.0D+300].
//
//    Input/output, double S, the number of successes.
//    If this is an input value, it should lie in the range: [0, +infinity).
//    If it is an output value, it will be searched for in the range:
//    [0, 1.0D+300].
//
//    Input/output, double PR, the probability of success in each
//    binomial trial.  Whether an input or output value, it should lie in the
//    range [0,1].
//
//    Input/output, double OMPR, the value of (1-PR).  Whether an
//    input or output value, it should lie in the range [0,1].
//
//    Output, int STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if PR + OMPR /= 1.
//
//    Output, double BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define inf 1.0e300
# define one 1.0e0

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double K11 = 1.0e0;
  static double fx,xhi,xlo,pq,prompr,cum,ccum;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10,T12,T13;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     S
//
    if(!(*s < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     XN
//
    if(!(*xn < 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 4) goto S190;
//
//     PR
//
    if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S180;
    if(!(*pr < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = -6;
    return;
S190:
S180:
    if(*which == 4) goto S230;
//
//     OMPR
//
    if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S220;
    if(!(*ompr < 0.0e0)) goto S200;
    *bound = 0.0e0;
    goto S210;
S200:
    *bound = 1.0e0;
S210:
    *status = -7;
    return;
S230:
S220:
    if(*which == 1) goto S270;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S260;
    if(!(pq < 0.0e0)) goto S240;
    *bound = 0.0e0;
    goto S250;
S240:
    *bound = 1.0e0;
S250:
    *status = 3;
    return;
S270:
S260:
    if(*which == 4) goto S310;
//
//     PR + OMPR
//
    prompr = *pr+*ompr;
    if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S300;
    if(!(prompr < 0.0e0)) goto S280;
    *bound = 0.0e0;
    goto S290;
S280:
    *bound = 1.0e0;
S290:
    *status = 4;
    return;
S310:
S300:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cumnbn(s,xn,pr,ompr,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating S
//
        *s = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        dinvr(status,s,&fx,&qleft,&qhi);
S320:
        if(!(*status == 1)) goto S350;
        cumnbn(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S330;
        fx = cum-*p;
        goto S340;
S330:
        fx = ccum-*q;
S340:
        dinvr(status,s,&fx,&qleft,&qhi);
        goto S320;
S350:
        if(!(*status == -1)) goto S380;
        if(!qleft) goto S360;
        *status = 1;
        *bound = 0.0e0;
        goto S370;
S360:
        *status = 2;
        *bound = inf;
S380:
S370:
        ;
    }
    else if(3 == *which) {
//
//     Calculating XN
//
        *xn = 5.0e0;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
        *status = 0;
        dinvr(status,xn,&fx,&qleft,&qhi);
S390:
        if(!(*status == 1)) goto S420;
        cumnbn(s,xn,pr,ompr,&cum,&ccum);
        if(!qporq) goto S400;
        fx = cum-*p;
        goto S410;
S400:
        fx = ccum-*q;
S410:
        dinvr(status,xn,&fx,&qleft,&qhi);
        goto S390;
S420:
        if(!(*status == -1)) goto S450;
        if(!qleft) goto S430;
        *status = 1;
        *bound = 0.0e0;
        goto S440;
S430:
        *status = 2;
        *bound = inf;
S450:
S440:
        ;
    }
    else if(4 == *which) {
//
//     Calculating PR and OMPR
//
        T12 = atol;
        T13 = tol;
        dstzr(&K2,&K11,&T12,&T13);
        if(!qporq) goto S480;
        *status = 0;
        dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
S460:
        if(!(*status == 1)) goto S470;
        cumnbn(s,xn,pr,ompr,&cum,&ccum);
        fx = cum-*p;
        dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
        *ompr = one-*pr;
        goto S460;
S470:
        goto S510;
S480:
        *status = 0;
        dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
S490:
        if(!(*status == 1)) goto S500;
        cumnbn(s,xn,pr,ompr,&cum,&ccum);
        fx = ccum-*q;
        dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
        *pr = one-*ompr;
        goto S490;
S510:
S500:
        if(!(*status == -1)) goto S540;
        if(!qleft) goto S520;
        *status = 1;
        *bound = 0.0e0;
        goto S530;
S520:
        *status = 2;
        *bound = 1.0e0;
S530:
        ;
    }
S540:
    return;
# undef tol
# undef atol
# undef inf
# undef one
}
//****************************************************************************80

void cdfnor ( int *which, double *p, double *q, double *x, double *mean,
  double *sd, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFNOR evaluates the CDF of the Normal distribution.
//
//  Discussion:
//
//    A slightly modified version of ANORM from SPECFUN
//    is used to calculate the cumulative standard normal distribution.
//
//    The rational functions from pages 90-95 of Kennedy and Gentle
//    are used as starting values to Newton's Iterations which
//    compute the inverse standard normal.  Therefore no searches are
//    necessary for any parameter.
//
//    For X < -15, the asymptotic expansion for the normal is used  as
//    the starting value in finding the inverse standard normal.
//
//    The normal density is proportional to
//    exp( - 0.5D+00 * (( X - MEAN)/SD)**2)
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.2.12.
//
//    William Cody,
//    Algorithm 715: SPECFUN - A Portable FORTRAN Package of
//      Special Function Routines and Test Drivers,
//    ACM Transactions on Mathematical Software,
//    Volume 19, pages 22-32, 1993.
//
//    Kennedy and Gentle,
//    Statistical Computing,
//    Marcel Dekker, NY, 1980,
//    QA276.4  K46
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, MEAN and SD;
//    2: Calculate X from P, Q, MEAN and SD;
//    3: Calculate MEAN from P, Q, X and SD;
//    4: Calculate SD from P, Q, X and MEAN.
//
//    Input/output, double *P, the integral from -infinity to X
//    of the Normal density.  If this is an input or output value, it will
//    lie in the range [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *X, the upper limit of integration of
//    the Normal density.
//
//    Input/output, double *MEAN, the mean of the Normal density.
//
//    Input/output, double *SD, the standard deviation of the
//    Normal density.  If this is an input value, it should lie in the
//    range (0,+infinity).
//
//    Output, int *STATUS, the status of the calculation.
//    0, if calculation completed correctly;
//    -I, if input parameter number I is out of range;
//    1, if answer appears to be lower than lowest search bound;
//    2, if answer appears to be higher than greatest search bound;
//    3, if P + Q /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
  static int K1 = 1;
  static double z,pq;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    *status = 0;
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p <= 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 1) goto S150;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S140;
    if(!(pq < 0.0e0)) goto S120;
    *bound = 0.0e0;
    goto S130;
S120:
    *bound = 1.0e0;
S130:
    *status = 3;
    return;
S150:
S140:
    if(*which == 4) goto S170;
//
//     SD
//
    if(!(*sd <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Computing P
//
        z = (*x-*mean)/ *sd;
        cumnor(&z,p,q);
    }
    else if(2 == *which) {
//
//     Computing X
//
        z = dinvnr(p,q);
        *x = *sd*z+*mean;
    }
    else if(3 == *which) {
//
//     Computing the MEAN
//
        z = dinvnr(p,q);
        *mean = *x-*sd*z;
    }
    else if(4 == *which) {
//
//     Computing SD
//
        z = dinvnr(p,q);
        *sd = (*x-*mean)/z;
    }
    return;
}
//****************************************************************************80

void cdfpoi ( int *which, double *p, double *q, double *s, double *xlam,
  int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFPOI evaluates the CDF of the Poisson distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the Poisson distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.4.21.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from S and XLAM;
//    2: Calculate A from P, Q and XLAM;
//    3: Calculate XLAM from P, Q and S.
//
//    Input/output, double *P, the cumulation from 0 to S of the
//    Poisson density.  Whether this is an input or output value, it will
//    lie in the range [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *S, the upper limit of cumulation of
//    the Poisson CDF.  If this is an input value, it should lie in
//    the range: [0, +infinity).  If it is an output value, it will be
//    searched for in the range: [0,1.0D+300].
//
//    Input/output, double *XLAM, the mean of the Poisson
//    distribution.  If this is an input value, it should lie in the range
//    [0, +infinity).  If it is an output value, it will be searched for
//    in the range: [0,1E300].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define inf 1.0e300

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double fx,cum,ccum,pq;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 3.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     S
//
    if(!(*s < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     XLAM
//
    if(!(*xlam < 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 1) goto S190;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S180;
    if(!(pq < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = 3;
    return;
S190:
S180:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cumpoi(s,xlam,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating S
//
        *s = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        dinvr(status,s,&fx,&qleft,&qhi);
S200:
        if(!(*status == 1)) goto S230;
        cumpoi(s,xlam,&cum,&ccum);
        if(!qporq) goto S210;
        fx = cum-*p;
        goto S220;
S210:
        fx = ccum-*q;
S220:
        dinvr(status,s,&fx,&qleft,&qhi);
        goto S200;
S230:
        if(!(*status == -1)) goto S260;
        if(!qleft) goto S240;
        *status = 1;
        *bound = 0.0e0;
        goto S250;
S240:
        *status = 2;
        *bound = inf;
S260:
S250:
        ;
    }
    else if(3 == *which) {
//
//     Calculating XLAM
//
        *xlam = 5.0e0;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
        *status = 0;
        dinvr(status,xlam,&fx,&qleft,&qhi);
S270:
        if(!(*status == 1)) goto S300;
        cumpoi(s,xlam,&cum,&ccum);
        if(!qporq) goto S280;
        fx = cum-*p;
        goto S290;
S280:
        fx = ccum-*q;
S290:
        dinvr(status,xlam,&fx,&qleft,&qhi);
        goto S270;
S300:
        if(!(*status == -1)) goto S330;
        if(!qleft) goto S310;
        *status = 1;
        *bound = 0.0e0;
        goto S320;
S310:
        *status = 2;
        *bound = inf;
S320:
        ;
    }
S330:
    return;
# undef tol
# undef atol
# undef inf
}
//****************************************************************************80

void cdft ( int *which, double *p, double *q, double *t, double *df,
  int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFT evaluates the CDF of the T distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the T distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.   The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The original version of this routine allowed the search interval
//    to extend from -1.0E+300 to +1.0E+300, which is fine until you
//    try to evaluate a function at such a point!
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.27.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1 : Calculate P and Q from T and DF;
//    2 : Calculate T from P, Q and DF;
//    3 : Calculate DF from P, Q and T.
//
//    Input/output, double *P, the integral from -infinity to T of
//    the T-density.  Whether an input or output value, this will lie in the
//    range [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *T, the upper limit of integration of
//    the T-density.  If this is an input value, it may have any value.
//    It it is an output value, it will be searched for in the range
//    [ -1.0D+30, 1.0D+30 ].
//
//    Input/output, double *DF, the number of degrees of freedom
//    of the T distribution.  If this is an input value, it should lie
//    in the range: (0 , +infinity).  If it is an output value, it will be
//    searched for in the range: [1, 1.0D+10].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define inf 1.0e30
# define maxdf 1.0e10

  static int K1 = 1;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double fx,cum,ccum,pq;
  static unsigned long qhi,qleft,qporq;
  static double T2,T3,T6,T7,T8,T9,T10,T11;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 3.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p <= 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 3) goto S130;
//
//     DF
//
    if(!(*df <= 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -5;
    return;
S130:
S120:
    if(*which == 1) goto S170;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S160;
    if(!(pq < 0.0e0)) goto S140;
    *bound = 0.0e0;
    goto S150;
S140:
    *bound = 1.0e0;
S150:
    *status = 3;
    return;
S170:
S160:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Computing P and Q
//
        cumt(t,df,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Computing T
//     .. Get initial approximation for T
//
        *t = dt1(p,q,df);
        T2 = -inf;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        dstinv(&T2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        dinvr(status,t,&fx,&qleft,&qhi);
S180:
        if(!(*status == 1)) goto S210;
        cumt(t,df,&cum,&ccum);
        if(!qporq) goto S190;
        fx = cum-*p;
        goto S200;
S190:
        fx = ccum-*q;
S200:
        dinvr(status,t,&fx,&qleft,&qhi);
        goto S180;
S210:
        if(!(*status == -1)) goto S240;
        if(!qleft) goto S220;
        *status = 1;
        *bound = -inf;
        goto S230;
S220:
        *status = 2;
        *bound = inf;
S240:
S230:
        ;
    }
    else if(3 == *which) {
//
//     Computing DF
//
        *df = 5.0e0;
        T8 = zero;
        T9 = maxdf;
        T10 = atol;
        T11 = tol;
        dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
        *status = 0;
        dinvr(status,df,&fx,&qleft,&qhi);
S250:
        if(!(*status == 1)) goto S280;
        cumt(t,df,&cum,&ccum);
        if(!qporq) goto S260;
        fx = cum-*p;
        goto S270;
S260:
        fx = ccum-*q;
S270:
        dinvr(status,df,&fx,&qleft,&qhi);
        goto S250;
S280:
        if(!(*status == -1)) goto S310;
        if(!qleft) goto S290;
        *status = 1;
        *bound = zero;
        goto S300;
S290:
        *status = 2;
        *bound = maxdf;
S300:
        ;
    }
S310:
    return;
# undef tol
# undef atol
# undef zero
# undef inf
# undef maxdf
}
//****************************************************************************80

void chi_noncentral_cdf_values ( int *n_data, double *x, double *lambda,
  int *df, double *cdf )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
//
//  Discussion:
//
//    The CDF of the noncentral chi square distribution can be evaluated
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *LAMBDA, the noncentrality parameter.
//
//    Output, int *DF, the number of degrees of freedom.
//
//    Output, double *CDF, the noncentral chi CDF.
//
{
# define N_MAX 27

  double cdf_vec[N_MAX] = {
    0.839944E+00, 0.695906E+00, 0.535088E+00,
    0.764784E+00, 0.620644E+00, 0.469167E+00,
    0.307088E+00, 0.220382E+00, 0.150025E+00,
    0.307116E-02, 0.176398E-02, 0.981679E-03,
    0.165175E-01, 0.202342E-03, 0.498448E-06,
    0.151325E-01, 0.209041E-02, 0.246502E-03,
    0.263684E-01, 0.185798E-01, 0.130574E-01,
    0.583804E-01, 0.424978E-01, 0.308214E-01,
    0.105788E+00, 0.794084E-01, 0.593201E-01 };
  int df_vec[N_MAX] = {
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
     60,  80, 100,
      1,   2,   3,
     10,  10,  10,
     10,  10,  10,
     10,  10,  10 };
  double lambda_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,
     1.0E+00,  1.0E+00,  1.0E+00,
     5.0E+00,  5.0E+00,  5.0E+00,
    20.0E+00, 20.0E+00, 20.0E+00,
    30.0E+00, 30.0E+00, 30.0E+00,
     5.0E+00,  5.0E+00,  5.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00 };
  double x_vec[N_MAX] = {
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
    60.000E+00, 60.000E+00, 60.000E+00,
     0.050E+00,  0.050E+00,  0.050E+00,
     4.000E+00,  4.000E+00,  4.000E+00,
     5.000E+00,  5.000E+00,  5.000E+00,
     6.000E+00,  6.000E+00,  6.000E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *lambda = 0.0E+00;
    *df = 0;
    *cdf = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *df = df_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void chi_square_cdf_values ( int *n_data, int *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
//
//  Discussion:
//
//    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
//    commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ChiSquareDistribution[DF], X ]
//
//  Modified:
//
//    11 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  int a_vec[N_MAX] = {
     1,  2,  1,  2,
     1,  2,  3,  4,
     1,  2,  3,  4,
     5,  3,  3,  3,
     3,  3, 10, 10,
    10 };
  double fx_vec[N_MAX] = {
    0.0796557E+00, 0.00498752E+00, 0.112463E+00,    0.00995017E+00,
    0.472911E+00,  0.181269E+00,   0.0597575E+00,   0.0175231E+00,
    0.682689E+00,  0.393469E+00,   0.198748E+00,    0.090204E+00,
    0.0374342E+00, 0.427593E+00,   0.608375E+00,    0.738536E+00,
    0.828203E+00,  0.88839E+00,    0.000172116E+00, 0.00365985E+00,
    0.0185759E+00 };
  double x_vec[N_MAX] = {
    0.01E+00, 0.01E+00, 0.02E+00, 0.02E+00,
    0.40E+00, 0.40E+00, 0.40E+00, 0.40E+00,
    1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
    1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
    5.00E+00, 6.00E+00, 1.00E+00, 2.00E+00,
    3.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

void cumbet ( double *x, double *y, double *a, double *b, double *cum,
  double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMBET evaluates the cumulative incomplete beta distribution.
//
//  Discussion:
//
//    This routine calculates the CDF to X of the incomplete beta distribution
//    with parameters A and B.  This is the integral from 0 to x
//    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
//
//  Modified:
//
//    14 March 2006
//
//  Reference:
//
//    A R Didonato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios.
//    ACM Transactions on Mathematical Software,
//    Volume 18, Number 3, September 1992, pages 360-373.
//
//  Parameters:
//
//    Input, double *X, the upper limit of integration.
//
//    Input, double *Y, the value of 1-X.
//
//    Input, double *A, *B, the parameters of the distribution.
//
//    Output, double *CUM, *CCUM, the values of the cumulative
//    density function and complementary cumulative density function.
//
{
  static int ierr;

  if ( *x <= 0.0 )
  {
    *cum = 0.0;
    *ccum = 1.0;
  }
  else if ( *y <= 0.0 )
  {
    *cum = 1.0;
    *ccum = 0.0;
  }
  else
  {
    beta_inc ( a, b, x, y, cum, ccum, &ierr );
  }
  return;
}
//****************************************************************************80

void cumbin ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMBIN evaluates the cumulative binomial distribution.
//
//  Discussion:
//
//    This routine returns the probability of 0 to S successes in XN binomial
//    trials, each of which has a probability of success, PR.
//
//  Modified:
//
//    14 March 2006
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.24.
//
//  Parameters:
//
//    Input, double *S, the upper limit of summation.
//
//    Input, double *XN, the number of trials.
//
//    Input, double *PR, the probability of success in one trial.
//
//    Input, double *OMPR, equals ( 1 - PR ).
//
//    Output, double *CUM, the cumulative binomial distribution.
//
//    Output, double *CCUM, the complement of the cumulative
//    binomial distribution.
//
{
  static double T1,T2;

  if ( *s < *xn )
  {
    T1 = *s + 1.0;
    T2 = *xn - *s;
    cumbet ( pr, ompr, &T1, &T2, ccum, cum );
  }
  else
  {
    *cum = 1.0;
    *ccum = 0.0;
  }
  return;
}
//****************************************************************************80

void cumchi ( double *x, double *df, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMCHI evaluates the cumulative chi-square distribution.
//
//  Parameters:
//
//    Input, double *X, the upper limit of integration.
//
//    Input, double *DF, the degrees of freedom of the
//    chi-square distribution.
//
//    Output, double *CUM, the cumulative chi-square distribution.
//
//    Output, double *CCUM, the complement of the cumulative
//    chi-square distribution.
//
{
  static double a;
  static double xx;

  a = *df * 0.5;
  xx = *x * 0.5;
  cumgam ( &xx, &a, cum, ccum );
  return;
}
//****************************************************************************80

void cumchn ( double *x, double *df, double *pnonc, double *cum,
  double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMCHN evaluates the cumulative noncentral chi-square distribution.
//
//  Discussion:
//
//    Calculates the cumulative noncentral chi-square
//    distribution, i.e., the probability that a random variable
//    which follows the noncentral chi-square distribution, with
//    noncentrality parameter PNONC and continuous degrees of
//    freedom DF, is less than or equal to X.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.4.25.
//
//  Parameters:
//
//    Input, double *X, the upper limit of integration.
//
//    Input, double *DF, the number of degrees of freedom.
//
//    Input, double *PNONC, the noncentrality parameter of
//    the noncentral chi-square distribution.
//
//    Output, double *CUM, *CCUM, the CDF and complementary
//    CDF of the noncentral chi-square distribution.
//
//  Local Parameters:
//
//    Local, double EPS, the convergence criterion.  The sum
//    stops when a term is less than EPS*SUM.
//
//    Local, int NTIRED, the maximum number of terms to be evaluated
//    in each sum.
//
//    Local, bool QCONV, is TRUE if convergence was achieved, that is,
//    the program did not stop on NTIRED criterion.
//
{
# define dg(i) (*df+2.0e0*(double)(i))
# define qsmall(xx) (int)(sum < 1.0e-20 || (xx) < eps*sum)
# define qtired(i) (int)((i) > ntired)

  static double eps = 1.0e-5;
  static int ntired = 1000;
  static double adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt,lfact,pcent,pterm,sum,
    sumadj,term,wt,xnonc;
  static int i,icent,iterb,iterf;
  static double T1,T2,T3;

    if(!(*x <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return;
S10:
    if(!(*pnonc <= 1.0e-10)) goto S20;
//
//     When non-centrality parameter is (essentially) zero,
//     use cumulative chi-square distribution
//
    cumchi(x,df,cum,ccum);
    return;
S20:
    xnonc = *pnonc/2.0e0;
//
//     The following code calculates the weight, chi-square, and
//     adjustment term for the central term in the infinite series.
//     The central term is the one in which the poisson weight is
//     greatest.  The adjustment term is the amount that must
//     be subtracted from the chi-square to move up two degrees
//     of freedom.
//
    icent = fifidint(xnonc);
    if(icent == 0) icent = 1;
    chid2 = *x/2.0e0;
//
//     Calculate central weight term
//
    T1 = (double)(icent+1);
    lfact = gamma_log ( &T1 );
    lcntwt = -xnonc+(double)icent*log(xnonc)-lfact;
    centwt = exp(lcntwt);
//
//     Calculate central chi-square
//
    T2 = dg(icent);
    cumchi(x,&T2,&pcent,ccum);
//
//     Calculate central adjustment term
//
    dfd2 = dg(icent)/2.0e0;
    T3 = 1.0e0+dfd2;
    lfact = gamma_log ( &T3 );
    lcntaj = dfd2*log(chid2)-chid2-lfact;
    centaj = exp(lcntaj);
    sum = centwt*pcent;
//
//     Sum backwards from the central term towards zero.
//     Quit whenever either
//     (1) the zero term is reached, or
//     (2) the term gets small relative to the sum, or
//     (3) More than NTIRED terms are totaled.
//
    iterb = 0;
    sumadj = 0.0e0;
    adj = centaj;
    wt = centwt;
    i = icent;
    goto S40;
S30:
    if(qtired(iterb) || qsmall(term) || i == 0) goto S50;
S40:
    dfd2 = dg(i)/2.0e0;
//
//     Adjust chi-square for two fewer degrees of freedom.
//     The adjusted value ends up in PTERM.
//
    adj = adj*dfd2/chid2;
    sumadj = sumadj + adj;
    pterm = pcent+sumadj;
//
//     Adjust poisson weight for J decreased by one
//
    wt *= ((double)i/xnonc);
    term = wt*pterm;
    sum = sum + term;
    i -= 1;
    iterb = iterb + 1;
    goto S30;
S50:
    iterf = 0;
//
//     Now sum forward from the central term towards infinity.
//     Quit when either
//     (1) the term gets small relative to the sum, or
//     (2) More than NTIRED terms are totaled.
//
    sumadj = adj = centaj;
    wt = centwt;
    i = icent;
    goto S70;
S60:
    if(qtired(iterf) || qsmall(term)) goto S80;
S70:
//
//     Update weights for next higher J
//
    wt *= (xnonc/(double)(i+1));
//
//     Calculate PTERM and add term to sum
//
    pterm = pcent-sumadj;
    term = wt*pterm;
    sum = sum + term;
//
//  Update adjustment term for DF for next iteration
//
    i = i + 1;
    dfd2 = dg(i)/2.0e0;
    adj = adj*chid2/dfd2;
    sumadj = sum + adj;
    iterf = iterf + 1;
    goto S60;
S80:
    *cum = sum;
    *ccum = 0.5e0+(0.5e0-*cum);
    return;
# undef dg
# undef qsmall
# undef qtired
}
//****************************************************************************80

void cumf ( double *f, double *dfn, double *dfd, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMF evaluates the cumulative F distribution.
//
//  Discussion:
//
//    CUMF computes the integral from 0 to F of the F density with DFN
//    numerator and DFD denominator degrees of freedom.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.28.
//
//  Parameters:
//
//    Input, double *F, the upper limit of integration.
//
//    Input, double *DFN, *DFD, the number of degrees of
//    freedom for the numerator and denominator.
//
//    Output, double *CUM, *CCUM, the value of the F CDF and
//    the complementary F CDF.
//
{
# define half 0.5e0
# define done 1.0e0

  static double dsum,prod,xx,yy;
  static int ierr;
  static double T1,T2;

  if(!(*f <= 0.0e0)) goto S10;
  *cum = 0.0e0;
  *ccum = 1.0e0;
  return;
S10:
  prod = *dfn**f;
//
//     XX is such that the incomplete beta with parameters
//     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
//     YY is 1 - XX
//     Calculate the smaller of XX and YY accurately
//
  dsum = *dfd+prod;
  xx = *dfd/dsum;

  if ( xx > half )
  {
    yy = prod/dsum;
    xx = done-yy;
  }
  else
  {
    yy = done-xx;
  }

  T1 = *dfd*half;
  T2 = *dfn*half;
  beta_inc ( &T1, &T2, &xx, &yy, ccum, cum, &ierr );
  return;
# undef half
# undef done
}
//****************************************************************************80

void cumfnc ( double *f, double *dfn, double *dfd, double *pnonc,
  double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMFNC evaluates the cumulative noncentral F distribution.
//
//  Discussion:
//
//    This routine computes the noncentral F distribution with DFN and DFD
//    degrees of freedom and noncentrality parameter PNONC.
//
//    The series is calculated backward and forward from J = LAMBDA/2
//    (this is the term with the largest Poisson weight) until
//    the convergence criterion is met.
//
//    The sum continues until a succeeding term is less than EPS
//    times the sum (or the sum is less than 1.0e-20).  EPS is
//    set to 1.0e-4 in a data statement which can be changed.
//
//
//    The original version of this routine allowed the input values
//    of DFN and DFD to be negative (nonsensical) or zero (which
//    caused numerical overflow.)  I have forced both these values
//    to be at least 1.
//
//  Modified:
//
//    15 June 2004
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
//
//  Parameters:
//
//    Input, double *F, the upper limit of integration.
//
//    Input, double *DFN, *DFD, the number of degrees of freedom
//    in the numerator and denominator.  Both DFN and DFD must be positive,
//    and normally would be integers.  This routine requires that they
//    be no less than 1.
//
//    Input, double *PNONC, the noncentrality parameter.
//
//    Output, double *CUM, *CCUM, the noncentral F CDF and
//    complementary CDF.
//
{
# define qsmall(x) (int)(sum < 1.0e-20 || (x) < eps*sum)
# define half 0.5e0
# define done 1.0e0

  static double eps = 1.0e-4;
  static double dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,
    upterm,xmult,xnonc;
  static int i,icent,ierr;
  static double T1,T2,T3,T4,T5,T6;

    if(!(*f <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return;
S10:
    if(!(*pnonc < 1.0e-10)) goto S20;
//
//  Handle case in which the non-centrality parameter is
//  (essentially) zero.
//
    cumf(f,dfn,dfd,cum,ccum);
    return;
S20:
    xnonc = *pnonc/2.0e0;
//
//  Calculate the central term of the poisson weighting factor.
//
    icent = ( int ) xnonc;
    if(icent == 0) icent = 1;
//
//  Compute central weight term
//
    T1 = (double)(icent+1);
    centwt = exp(-xnonc+(double)icent*log(xnonc)- gamma_log ( &T1 ) );
//
//  Compute central incomplete beta term
//  Assure that minimum of arg to beta and 1 - arg is computed
//  accurately.
//
    prod = *dfn**f;
    dsum = *dfd+prod;
    yy = *dfd/dsum;
    if(yy > half) {
        xx = prod/dsum;
        yy = done-xx;
    }
    else  xx = done-yy;
    T2 = *dfn*half+(double)icent;
    T3 = *dfd*half;
    beta_inc ( &T2, &T3, &xx, &yy, &betdn, &dummy, &ierr );
    adn = *dfn/2.0e0+(double)icent;
    aup = adn;
    b = *dfd/2.0e0;
    betup = betdn;
    sum = centwt*betdn;
//
//  Now sum terms backward from icent until convergence or all done
//
    xmult = centwt;
    i = icent;
    T4 = adn+b;
    T5 = adn+1.0e0;
    dnterm = exp( gamma_log ( &T4 ) - gamma_log ( &T5 )
      - gamma_log ( &b ) + adn * log ( xx ) + b * log(yy));
S30:
    if(qsmall(xmult*betdn) || i <= 0) goto S40;
    xmult *= ((double)i/xnonc);
    i -= 1;
    adn -= 1.0;
    dnterm = (adn+1.0)/((adn+b)*xx)*dnterm;
    betdn += dnterm;
    sum += (xmult*betdn);
    goto S30;
S40:
    i = icent+1;
//
//  Now sum forwards until convergence
//
    xmult = centwt;
    if(aup-1.0+b == 0) upterm = exp(-gamma_log ( &aup )
      - gamma_log ( &b ) + (aup-1.0)*log(xx)+
      b*log(yy));
    else  {
        T6 = aup-1.0+b;
        upterm = exp( gamma_log ( &T6 ) - gamma_log ( &aup )
          - gamma_log ( &b ) + (aup-1.0)*log(xx)+b*
          log(yy));
    }
    goto S60;
S50:
    if(qsmall(xmult*betup)) goto S70;
S60:
    xmult *= (xnonc/(double)i);
    i += 1;
    aup += 1.0;
    upterm = (aup+b-2.0e0)*xx/(aup-1.0)*upterm;
    betup -= upterm;
    sum += (xmult*betup);
    goto S50;
S70:
    *cum = sum;
    *ccum = 0.5e0+(0.5e0-*cum);
    return;
# undef qsmall
# undef half
# undef done
}
//****************************************************************************80

void cumgam ( double *x, double *a, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMGAM evaluates the cumulative incomplete gamma distribution.
//
//  Discussion:
//
//    This routine computes the cumulative distribution function of the
//    incomplete gamma distribution, i.e., the integral from 0 to X of
//
//      (1/GAM(A))*EXP(-T)*T**(A-1) DT
//
//    where GAM(A) is the complete gamma function of A, i.e.,
//
//      GAM(A) = integral from 0 to infinity of EXP(-T)*T**(A-1) DT
//
//  Parameters:
//
//    Input, double *X, the upper limit of integration.
//
//    Input, double *A, the shape parameter of the incomplete
//    Gamma distribution.
//
//    Output, double *CUM, *CCUM, the incomplete Gamma CDF and
//    complementary CDF.
//
{
  static int K1 = 0;

  if(!(*x <= 0.0e0)) goto S10;
  *cum = 0.0e0;
  *ccum = 1.0e0;
  return;
S10:
  gamma_inc ( a, x, cum, ccum, &K1 );
//
//     Call gratio routine
//
    return;
}
//****************************************************************************80

void cumnbn ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMNBN evaluates the cumulative negative binomial distribution.
//
//  Discussion:
//
//    This routine returns the probability that there will be F or
//    fewer failures before there are S successes, with each binomial
//    trial having a probability of success PR.
//
//    Prob(# failures = F | S successes, PR)  =
//                        ( S + F - 1 )
//                        (            ) * PR^S * (1-PR)^F
//                        (      F     )
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.26.
//
//  Parameters:
//
//    Input, double *F, the number of failures.
//
//    Input, double *S, the number of successes.
//
//    Input, double *PR, *OMPR, the probability of success on
//    each binomial trial, and the value of (1-PR).
//
//    Output, double *CUM, *CCUM, the negative binomial CDF,
//    and the complementary CDF.
//
{
  static double T1;

  T1 = *s+1.e0;
  cumbet(pr,ompr,xn,&T1,cum,ccum);
  return;
}
//****************************************************************************80

void cumnor ( double *arg, double *result, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMNOR computes the cumulative normal distribution.
//
//  Discussion:
//
//    This function evaluates the normal distribution function:
//
//                              / x
//                     1       |       -t*t/2
//          P(x) = ----------- |      e       dt
//                 sqrt(2 pi)  |
//                             /-oo
//
//    This transportable program uses rational functions that
//    theoretically approximate the normal distribution function to
//    at least 18 significant decimal digits.  The accuracy achieved
//    depends on the arithmetic system, the compiler, the intrinsic
//    functions, and proper selection of the machine-dependent
//    constants.
//
//  Author:
//
//    William Cody
//    Mathematics and Computer Science Division
//    Argonne National Laboratory
//    Argonne, IL 60439
//
//  Reference:
//
//    William Cody,
//    Rational Chebyshev approximations for the error function,
//    Mathematics of Computation,
//    1969, pages 631-637.
//
//    William Cody,
//    Algorithm 715:
//    SPECFUN - A Portable FORTRAN Package of Special Function Routines
//      and Test Drivers,
//    ACM Transactions on Mathematical Software,
//    Volume 19, 1993, pages 22-32.
//
//  Parameters:
//
//    Input, double *ARG, the upper limit of integration.
//
//    Output, double *CUM, *CCUM, the Normal density CDF and
//    complementary CDF.
//
//  Local Parameters:
//
//    Local, double EPS, the argument below which anorm(x)
//    may be represented by 0.5D+00 and above which  x*x  will not underflow.
//    A conservative value is the largest machine number X
//    such that   1.0D+00 + X = 1.0D+00   to machine precision.
//
{
  static double a[5] = {
    2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
    1.8154981253343561249e04,6.5682337918207449113e-2
  };
  static double b[4] = {
    4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
    4.5507789335026729956e04
  };
  static double c[9] = {
    3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
    5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
    1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
  };
  static double d[8] = {
    2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
    6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
    3.8912003286093271411e04,1.9685429676859990727e04
  };
  static double half = 0.5e0;
  static double p[6] = {
    2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
    1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
  };
  static double one = 1.0e0;
  static double q[5] = {
    1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
    3.78239633202758244e-3,7.29751555083966205e-5
  };
  static double sixten = 1.60e0;
  static double sqrpi = 3.9894228040143267794e-1;
  static double thrsh = 0.66291e0;
  static double root32 = 5.656854248e0;
  static double zero = 0.0e0;
  static int K1 = 1;
  static int K2 = 2;
  static int i;
  static double del,eps,temp,x,xden,xnum,y,xsq,min;
//
//  Machine dependent constants
//
    eps = dpmpar(&K1)*0.5e0;
    min = dpmpar(&K2);
    x = *arg;
    y = fabs(x);
    if(y <= thrsh) {
//
//  Evaluate  anorm  for  |X| <= 0.66291
//
        xsq = zero;
        if(y > eps) xsq = x*x;
        xnum = a[4]*xsq;
        xden = xsq;
        for ( i = 0; i < 3; i++ )
        {
            xnum = (xnum+a[i])*xsq;
            xden = (xden+b[i])*xsq;
        }
        *result = x*(xnum+a[3])/(xden+b[3]);
        temp = *result;
        *result = half+temp;
        *ccum = half-temp;
    }
//
//  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
//
    else if(y <= root32) {
        xnum = c[8]*y;
        xden = y;
        for ( i = 0; i < 7; i++ )
        {
            xnum = (xnum+c[i])*y;
            xden = (xden+d[i])*y;
        }
        *result = (xnum+c[7])/(xden+d[7]);
        xsq = fifdint(y*sixten)/sixten;
        del = (y-xsq)*(y+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
//
//  Evaluate  anorm  for |X| > sqrt(32)
//
    else  {
        *result = zero;
        xsq = one/(x*x);
        xnum = p[5]*xsq;
        xden = xsq;
        for ( i = 0; i < 4; i++ )
        {
            xnum = (xnum+p[i])*xsq;
            xden = (xden+q[i])*xsq;
        }
        *result = xsq*(xnum+p[4])/(xden+q[4]);
        *result = (sqrpi-*result)/y;
        xsq = fifdint(x*sixten)/sixten;
        del = (x-xsq)*(x+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
    if(*result < min) *result = 0.0e0;
//
//  Fix up for negative argument, erf, etc.
//
    if(*ccum < min) *ccum = 0.0e0;
}
//****************************************************************************80

void cumpoi ( double *s, double *xlam, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMPOI evaluates the cumulative Poisson distribution.
//
//  Discussion:
//
//    CUMPOI returns the probability of S or fewer events in a Poisson
//    distribution with mean XLAM.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.4.21.
//
//  Parameters:
//
//    Input, double *S, the upper limit of cumulation of the
//    Poisson density function.
//
//    Input, double *XLAM, the mean of the Poisson distribution.
//
//    Output, double *CUM, *CCUM, the Poisson density CDF and
//    complementary CDF.
//
{
  static double chi,df;

  df = 2.0e0*(*s+1.0e0);
  chi = 2.0e0**xlam;
  cumchi(&chi,&df,ccum,cum);
  return;
}
//****************************************************************************80

void cumt ( double *t, double *df, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMT evaluates the cumulative T distribution.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.5.27.
//
//  Parameters:
//
//    Input, double *T, the upper limit of integration.
//
//    Input, double *DF, the number of degrees of freedom of
//    the T distribution.
//
//    Output, double *CUM, *CCUM, the T distribution CDF and
//    complementary CDF.
//
{
  static double a;
  static double dfptt;
  static double K2 = 0.5e0;
  static double oma;
  static double T1;
  static double tt;
  static double xx;
  static double yy;

  tt = (*t) * (*t);
  dfptt = ( *df ) + tt;
  xx = *df / dfptt;
  yy = tt / dfptt;
  T1 = 0.5e0 * ( *df );
  cumbet ( &xx, &yy, &T1, &K2, &a, &oma );

  if ( *t <= 0.0e0 )
  {
    *cum = 0.5e0 * a;
    *ccum = oma + ( *cum );
  }
  else
  {
    *ccum = 0.5e0 * a;
    *cum = oma + ( *ccum );
  }
  return;
}
//****************************************************************************80

double dbetrm ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    DBETRM computes the Sterling remainder for the complete beta function.
//
//  Discussion:
//
//    Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
//    where Lgamma is the log of the (complete) gamma function
//
//    Let ZZ be approximation obtained if each log gamma is approximated
//    by Sterling's formula, i.e.,
//    Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5D+00 ) * LOG( Z ) - Z
//
//    The Sterling remainder is Log(Beta(A,B)) - ZZ.
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the Beta function.
//
//    Output, double DBETRM, the Sterling remainder.
//
{
  static double dbetrm,T1,T2,T3;
//
//     Try to sum from smallest to largest
//
    T1 = *a+*b;
    dbetrm = -dstrem(&T1);
    T2 = fifdmax1(*a,*b);
    dbetrm += dstrem(&T2);
    T3 = fifdmin1(*a,*b);
    dbetrm += dstrem(&T3);
    return dbetrm;
}
//****************************************************************************80

double dexpm1 ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    DEXPM1 evaluates the function EXP(X) - 1.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *X, the value at which exp(X)-1 is desired.
//
//    Output, double DEXPM1, the value of exp(X)-1.
//
{
  static double p1 = .914041914819518e-09;
  static double p2 = .238082361044469e-01;
  static double q1 = -.499999999085958e+00;
  static double q2 = .107141568980644e+00;
  static double q3 = -.119041179760821e-01;
  static double q4 = .595130811860248e-03;
  static double dexpm1;
  double w;

  if ( fabs(*x) <= 0.15e0 )
  {
    dexpm1 =   *x * ( ( (
        p2   * *x
      + p1 ) * *x
      + 1.0e0 )
      /((((
        q4   * *x
      + q3 ) * *x
      + q2 ) * *x
      + q1 ) * *x
      + 1.0e0 ) );
  }
  else if ( *x <= 0.0e0 )
  {
    w = exp(*x);
    dexpm1 = w-0.5e0-0.5e0;
  }
  else
  {
    w = exp(*x);
    dexpm1 = w*(0.5e0+(0.5e0-1.0e0/w));
  }

  return dexpm1;
}
//****************************************************************************80

double dinvnr ( double *p, double *q )

//****************************************************************************80
//
//  Purpose:
//
//    DINVNR computes the inverse of the normal distribution.
//
//  Discussion:
//
//    Returns X such that CUMNOR(X)  =   P,  i.e., the  integral from -
//    infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
//
//    The rational function on page 95 of Kennedy and Gentle is used as a start
//    value for the Newton method of finding roots.
//
//  Reference:
//
//    Kennedy and Gentle,
//    Statistical Computing,
//    Marcel Dekker, NY, 1980,
//    QA276.4  K46
//
//  Parameters:
//
//    Input, double *P, *Q, the probability, and the complementary
//    probability.
//
//    Output, double DINVNR, the argument X for which the
//    Normal CDF has the value P.
//
{
# define maxit 100
# define eps (1.0e-13)
# define r2pi 0.3989422804014326e0
# define nhalf (-0.5e0)
# define dennor(x) (r2pi*exp(nhalf*(x)*(x)))

  static double dinvnr,strtx,xcur,cum,ccum,pp,dx;
  static int i;
  static unsigned long qporq;

//
//     FIND MINIMUM OF P AND Q
//
    qporq = *p <= *q;
    if(!qporq) goto S10;
    pp = *p;
    goto S20;
S10:
    pp = *q;
S20:
//
//     INITIALIZATION STEP
//
    strtx = stvaln(&pp);
    xcur = strtx;
//
//     NEWTON INTERATIONS
//
    for ( i = 1; i <= maxit; i++ )
    {
        cumnor(&xcur,&cum,&ccum);
        dx = (cum-pp)/dennor(xcur);
        xcur -= dx;
        if(fabs(dx/xcur) < eps) goto S40;
    }
    dinvnr = strtx;
//
//     IF WE GET HERE, NEWTON HAS FAILED
//
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
S40:
//
//     IF WE GET HERE, NEWTON HAS SUCCEDED
//
    dinvnr = xcur;
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
# undef maxit
# undef eps
# undef r2pi
# undef nhalf
# undef dennor
}
//****************************************************************************80

void dinvr ( int *status, double *x, double *fx,
  unsigned long *qleft, unsigned long *qhi )

//****************************************************************************80
//
//  Purpose:
//
//    DINVR bounds the zero of the function and invokes DZROR.
//
//  Discussion:
//
//    This routine seeks to find bounds on a root of the function and
//    invokes ZROR to perform the zero finding.  STINVR must have been
//    called before this routine in order to set its parameters.
//
//  Reference:
//
//    J C P Bus and T J Dekker,
//    Two Efficient Algorithms with Guaranteed Convergence for
//      Finding a Zero of a Function,
//    ACM Transactions on Mathematical Software,
//    Volume 1, Number 4, pages 330-345, 1975.
//
//  Parameters:
//
//    Input/output, integer STATUS.  At the beginning of a zero finding
//    problem, STATUS should be set to 0 and INVR invoked.  The value
//    of parameters other than X will be ignored on this call.
//    If INVR needs the function to be evaluated, it will set STATUS to 1
//    and return.  The value of the function should be set in FX and INVR
//    again called without changing any of its other parameters.
//    If INVR finishes without error, it returns with STATUS 0, and X an
//    approximate root of F(X).
//    If INVR cannot bound the function, it returns a negative STATUS and
//    sets QLEFT and QHI.
//
//    Output, double precision X, the value at which F(X) is to be evaluated.
//
//    Input, double precision FX, the value of F(X) calculated by the user
//    on the previous call, when INVR returned with STATUS = 1.
//
//    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
//    case, QLEFT is TRUE if the stepping search terminated unsucessfully
//    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
//
//    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
//    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
//    if F(X) < Y.
//
{
  E0000(0,status,x,fx,qleft,qhi,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}
//****************************************************************************80

double dlanor ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    DLANOR evaluates the logarithm of the asymptotic Normal CDF.
//
//  Discussion:
//
//    This routine computes the logarithm of the cumulative normal distribution
//    from abs ( x ) to infinity for  5 <= abs ( X ).
//
//    The relative error at X = 5 is about 0.5D-5.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.2.12.
//
//  Parameters:
//
//    Input, double *X, the value at which the Normal CDF is to be
//    evaluated.  It is assumed that 5 <= abs ( X ).
//
//    Output, double DLANOR, the logarithm of the asymptotic
//    Normal CDF.
//
{
# define dlsqpi 0.91893853320467274177e0

  static double coef[12] = {
    -1.0e0,3.0e0,-15.0e0,105.0e0,-945.0e0,10395.0e0,-135135.0e0,2027025.0e0,
    -34459425.0e0,654729075.0e0,-13749310575.e0,316234143225.0e0
  };
  static int K1 = 12;
  static double dlanor,approx,correc,xx,xx2,T2;

  xx = fabs(*x);
  if ( xx < 5.0e0 ) 
  {
    ftnstop(" Argument too small in DLANOR");
  }
  approx = -dlsqpi-0.5e0*xx*xx-log(xx);
  xx2 = xx*xx;
  T2 = 1.0e0/xx2;
  correc = eval_pol ( coef, &K1, &T2 ) / xx2;
  correc = alnrel ( &correc );
  dlanor = approx+correc;
  return dlanor;
# undef dlsqpi
}
//****************************************************************************80

double dpmpar ( int *i )

//****************************************************************************80
//
//  Purpose:
//
//    DPMPAR provides machine constants for double precision arithmetic.
//
//  Discussion:
//
//     DPMPAR PROVIDES THE double PRECISION MACHINE CONSTANTS FOR
//     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
//     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
//     double PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
//     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
//
//        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
//
//        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
//
//        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
//
//     WRITTEN BY
//        ALFRED H. MORRIS, JR.
//        NAVAL SURFACE WARFARE CENTER
//        DAHLGREN VIRGINIA
//
//     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
//     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
//     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
//
{
  static int K1 = 4;
  static int K2 = 8;
  static int K3 = 9;
  static int K4 = 10;
  static double value,b,binv,bm1,one,w,z;
  static int emax,emin,ibeta,m;

    if(*i > 1) goto S10;
    b = ipmpar(&K1);
    m = ipmpar(&K2);
    value = pow(b,(double)(1-m));
    return value;
S10:
    if(*i > 2) goto S20;
    b = ipmpar(&K1);
    emin = ipmpar(&K3);
    one = 1.0;
    binv = one/b;
    w = pow(b,(double)(emin+2));
    value = w*binv*binv*binv;
    return value;
S20:
    ibeta = ipmpar(&K1);
    m = ipmpar(&K2);
    emax = ipmpar(&K4);
    b = ibeta;
    bm1 = ibeta-1;
    one = 1.0;
    z = pow(b,(double)(m-1));
    w = ((z-one)*b+bm1)/(b*z);
    z = pow(b,(double)(emax-2));
    value = w*z*b*b;
    return value;
}
//****************************************************************************80

void dstinv ( double *zsmall, double *zbig, double *zabsst,
  double *zrelst, double *zstpmu, double *zabsto, double *zrelto )

//****************************************************************************80
//
//  Purpose:
//
//    DSTINV seeks a value X such that F(X) = Y.
//
//  Discussion:
//
//      Double Precision - SeT INverse finder - Reverse Communication
//                              Function
//     Concise Description - Given a monotone function F finds X
//     such that F(X) = Y.  Uses Reverse communication -- see invr.
//     This routine sets quantities needed by INVR.
//          More Precise Description of INVR -
//     F must be a monotone function, the results of QMFINV are
//     otherwise undefined.  QINCR must be .TRUE. if F is non-
//     decreasing and .FALSE. if F is non-increasing.
//     QMFINV will return .TRUE. if and only if F(SMALL) and
//     F(BIG) bracket Y, i. e.,
//          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
//          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
//     if QMFINV returns .TRUE., then the X returned satisfies
//     the following condition.  let
//               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
//     then if QINCR is .TRUE.,
//          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
//     and if QINCR is .FALSE.
//          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
//                              Arguments
//     SMALL --> The left endpoint of the interval to be
//          searched for a solution.
//                    SMALL is DOUBLE PRECISION
//     BIG --> The right endpoint of the interval to be
//          searched for a solution.
//                    BIG is DOUBLE PRECISION
//     ABSSTP, RELSTP --> The initial step size in the search
//          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
//                    ABSSTP is DOUBLE PRECISION
//                    RELSTP is DOUBLE PRECISION
//     STPMUL --> When a step doesn't bound the zero, the step
//                size is multiplied by STPMUL and another step
//                taken.  A popular value is 2.0
//                    DOUBLE PRECISION STPMUL
//     ABSTOL, RELTOL --> Two numbers that determine the accuracy
//          of the solution.  See function for a precise definition.
//                    ABSTOL is DOUBLE PRECISION
//                    RELTOL is DOUBLE PRECISION
//                              Method
//     Compares F(X) with Y for the input value of X then uses QINCR
//     to determine whether to step left or right to bound the
//     desired x.  the initial step size is
//          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
//     Iteratively steps right or left until it bounds X.
//     At each step which doesn't bound X, the step size is doubled.
//     The routine is careful never to step beyond SMALL or BIG.  If
//     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
//     after setting QLEFT and QHI.
//     If X is successfully bounded then Algorithm R of the paper
//     'Two Efficient Algorithms with Guaranteed Convergence for
//     Finding a Zero of a Function' by J. C. P. Bus and
//     T. J. Dekker in ACM Transactions on Mathematical
//     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
//     to find the zero of the function F(X)-Y. This is routine
//     QRZERO.
//
{
  E0000(1,NULL,NULL,NULL,NULL,NULL,zabsst,zabsto,zbig,zrelst,zrelto,zsmall,
    zstpmu);
}
//****************************************************************************80

double dstrem ( double *z )

//****************************************************************************80
//
//  Purpose:
//
//    DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
//
//  Discussion:
//
//    This routine returns
//
//      ln ( Gamma ( Z ) ) - Sterling ( Z )
//
//    where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
//
//    Sterling(Z) = ln ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
//
//    If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
//    with values calculated using Maple.
//
//    Otherwise, the difference is computed explicitly.
//
//  Modified:
//
//    14 June 2004
//
//  Parameters:
//
//    Input, double *Z, the value at which the Sterling
//    remainder is to be calculated.  Z must be positive.
//
//    Output, double DSTREM, the Sterling remainder.
//
{
# define hln2pi 0.91893853320467274178e0
# define ncoef 10

  static double coef[ncoef] = {
    0.0e0,0.0833333333333333333333333333333e0,
    -0.00277777777777777777777777777778e0,0.000793650793650793650793650793651e0,
    -0.000595238095238095238095238095238e0,
    0.000841750841750841750841750841751e0,-0.00191752691752691752691752691753e0,
    0.00641025641025641025641025641026e0,-0.0295506535947712418300653594771e0,
    0.179644372368830573164938490016e0
  };
  static int K1 = 10;
  static double dstrem,sterl,T2;
//
//    For information, here are the next 11 coefficients of the
//    remainder term in Sterling's formula
//            -1.39243221690590111642743221691
//            13.4028640441683919944789510007
//            -156.848284626002017306365132452
//            2193.10333333333333333333333333
//            -36108.7712537249893571732652192
//            691472.268851313067108395250776
//            -0.152382215394074161922833649589D8
//            0.382900751391414141414141414141D9
//            -0.108822660357843910890151491655D11
//            0.347320283765002252252252252252D12
//            -0.123696021422692744542517103493D14
//
    if(*z <= 0.0e0) 
    {
      ftnstop ( "Zero or negative argument in DSTREM" );
    }
    if(!(*z > 6.0e0)) goto S10;
    T2 = 1.0e0/pow(*z,2.0);
    dstrem = eval_pol ( coef, &K1, &T2 )**z;
    goto S20;
S10:
    sterl = hln2pi+(*z-0.5e0)*log(*z)-*z;
    dstrem = gamma_log ( z ) - sterl;
S20:
    return dstrem;
# undef hln2pi
# undef ncoef
}
//****************************************************************************80

void dstzr ( double *zxlo, double *zxhi, double *zabstl, double *zreltl )

//****************************************************************************80
//
//  Purpose:
//
//    DSTXR sets quantities needed by the zero finder.
//
//  Discussion:
//
//     Double precision SeT ZeRo finder - Reverse communication version
//                              Function
//     Sets quantities needed by ZROR.  The function of ZROR
//     and the quantities set is given here.
//     Concise Description - Given a function F
//     find XLO such that F(XLO) = 0.
//          More Precise Description -
//     Input condition. F is a double function of a single
//     double argument and XLO and XHI are such that
//          F(XLO)*F(XHI)  .LE.  0.0
//     If the input condition is met, QRZERO returns .TRUE.
//     and output values of XLO and XHI satisfy the following
//          F(XLO)*F(XHI)  .LE. 0.
//          ABS(F(XLO)  .LE. ABS(F(XHI)
//          ABS(XLO-XHI)  .LE. TOL(X)
//     where
//          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
//     If this algorithm does not find XLO and XHI satisfying
//     these conditions then QRZERO returns .FALSE.  This
//     implies that the input condition was not met.
//                              Arguments
//     XLO --> The left endpoint of the interval to be
//           searched for a solution.
//                    XLO is DOUBLE PRECISION
//     XHI --> The right endpoint of the interval to be
//           for a solution.
//                    XHI is DOUBLE PRECISION
//     ABSTOL, RELTOL --> Two numbers that determine the accuracy
//                      of the solution.  See function for a
//                      precise definition.
//                    ABSTOL is DOUBLE PRECISION
//                    RELTOL is DOUBLE PRECISION
//                              Method
//     Algorithm R of the paper 'Two Efficient Algorithms with
//     Guaranteed Convergence for Finding a Zero of a Function'
//     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
//     Mathematical Software, Volume 1, no. 4 page 330
//     (Dec. '75) is employed to find the zero of F(X)-Y.
//
{
  E0001(1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,zabstl,zreltl,zxhi,zxlo);
}
//****************************************************************************80

double dt1 ( double *p, double *q, double *df )

//****************************************************************************80
//
//  Purpose:
//
//    DT1 computes an approximate inverse of the cumulative T distribution.
//
//  Discussion:
//
//    Returns the inverse of the T distribution function, i.e.,
//    the integral from 0 to INVT of the T density is P. This is an
//    initial approximation.
//
//  Parameters:
//
//    Input, double *P, *Q, the value whose inverse from the
//    T distribution CDF is desired, and the value (1-P).
//
//    Input, double *DF, the number of degrees of freedom of the
//    T distribution.
//
//    Output, double DT1, the approximate value of X for which
//    the T density CDF with DF degrees of freedom has value P.
//
{
  static double coef[4][5] = {
    1.0e0,1.0e0,0.0e0,0.0e0,0.0e0,3.0e0,16.0e0,5.0e0,0.0e0,0.0e0,-15.0e0,17.0e0,
    19.0e0,3.0e0,0.0e0,-945.0e0,-1920.0e0,1482.0e0,776.0e0,79.0e0
  };
  static double denom[4] = {
    4.0e0,96.0e0,384.0e0,92160.0e0
  };
  static int ideg[4] = {
    2,3,4,5
  };
  static double dt1,denpow,sum,term,x,xp,xx;
  static int i;

    x = fabs(dinvnr(p,q));
    xx = x*x;
    sum = x;
    denpow = 1.0e0;
    for ( i = 0; i < 4; i++ )
    {
        term = eval_pol ( &coef[i][0], &ideg[i], &xx ) * x;
        denpow *= *df;
        sum += (term/(denpow*denom[i]));
    }
    if(!(*p >= 0.5e0)) goto S20;
    xp = sum;
    goto S30;
S20:
    xp = -sum;
S30:
    dt1 = xp;
    return dt1;
}
//****************************************************************************80

void dzror ( int *status, double *x, double *fx, double *xlo,
  double *xhi, unsigned long *qleft, unsigned long *qhi )

//****************************************************************************80
//
//  Purpose:
//
//    DZROR seeks the zero of a function using reverse communication.
//
//  Discussion:
//
//     Performs the zero finding.  STZROR must have been called before
//     this routine in order to set its parameters.
//
//
//                              Arguments
//
//
//     STATUS <--> At the beginning of a zero finding problem, STATUS
//                 should be set to 0 and ZROR invoked.  (The value
//                 of other parameters will be ignored on this call.)
//
//                 When ZROR needs the function evaluated, it will set
//                 STATUS to 1 and return.  The value of the function
//                 should be set in FX and ZROR again called without
//                 changing any of its other parameters.
//
//                 When ZROR has finished without error, it will return
//                 with STATUS 0.  In that case (XLO,XHI) bound the answe
//
//                 If ZROR finds an error (which implies that F(XLO)-Y an
//                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
//                 this case, XLO and XHI are undefined.
//                         INTEGER STATUS
//
//     X <-- The value of X at which F(X) is to be evaluated.
//                         DOUBLE PRECISION X
//
//     FX --> The value of F(X) calculated when ZROR returns with
//            STATUS = 1.
//                         DOUBLE PRECISION FX
//
//     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
//             inverval in X containing the solution below.
//                         DOUBLE PRECISION XLO
//
//     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
//             inverval in X containing the solution above.
//                         DOUBLE PRECISION XHI
//
//     QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
//                at XLO.  If it is .FALSE. the search terminated
//                unsucessfully at XHI.
//                    QLEFT is LOGICAL
//
//     QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
//              search and .FALSE. if F(X) .LT. Y at the
//              termination of the search.
//                    QHI is LOGICAL
//
//
{
  E0001(0,status,x,fx,xlo,xhi,qleft,qhi,NULL,NULL,NULL,NULL);
}
//****************************************************************************80

void E0000 ( int IENTRY, int *status, double *x, double *fx,
  unsigned long *qleft, unsigned long *qhi, double *zabsst,
  double *zabsto, double *zbig, double *zrelst,
  double *zrelto, double *zsmall, double *zstpmu )

//****************************************************************************80
//
//  Purpose:
//
//    E0000 is a reverse-communication zero bounder.
//
{
# define qxmon(zx,zy,zz) (int)((zx) <= (zy) && (zy) <= (zz))

  static double absstp;
  static double abstol;
  static double big,fbig,fsmall,relstp,reltol,small,step,stpmul,xhi,
    xlb,xlo,xsave,xub,yy;
  static int i99999;
  static unsigned long qbdd,qcond,qdum1,qdum2,qincr,qlim,qok,qup;
    switch(IENTRY){case 0: goto DINVR; case 1: goto DSTINV;}
DINVR:
    if(*status > 0) goto S310;
    qcond = !qxmon(small,*x,big);
    if(qcond) 
    {
      ftnstop(" SMALL, X, BIG not monotone in INVR");
    }
    xsave = *x;
//
//     See that SMALL and BIG bound the zero and set QINCR
//
    *x = small;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 1;
    goto S300;
S10:
    fsmall = *fx;
    *x = big;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 2;
    goto S300;
S20:
    fbig = *fx;
    qincr = fbig > fsmall;
    if(!qincr) goto S50;
    if(fsmall <= 0.0e0) goto S30;
    *status = -1;
    *qleft = *qhi = 1;
    return;
S30:
    if(fbig >= 0.0e0) goto S40;
    *status = -1;
    *qleft = *qhi = 0;
    return;
S40:
    goto S80;
S50:
    if(fsmall >= 0.0e0) goto S60;
    *status = -1;
    *qleft = 1;
    *qhi = 0;
    return;
S60:
    if(fbig <= 0.0e0) goto S70;
    *status = -1;
    *qleft = 0;
    *qhi = 1;
    return;
S80:
S70:
    *x = xsave;
    step = fifdmax1(absstp,relstp*fabs(*x));
//
//      YY = F(X) - Y
//     GET-FUNCTION-VALUE
//
    i99999 = 3;
    goto S300;
S90:
    yy = *fx;
    if(!(yy == 0.0e0)) goto S100;
    *status = 0;
    qok = 1;
    return;
S100:
    qup = qincr && yy < 0.0e0 || !qincr && yy > 0.0e0;
//
//     HANDLE CASE IN WHICH WE MUST STEP HIGHER
//
    if(!qup) goto S170;
    xlb = xsave;
    xub = fifdmin1(xlb+step,big);
    goto S120;
S110:
    if(qcond) goto S150;
S120:
//
//      YY = F(XUB) - Y
//
    *x = xub;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 4;
    goto S300;
S130:
    yy = *fx;
    qbdd = qincr && yy >= 0.0e0 || !qincr && yy <= 0.0e0;
    qlim = xub >= big;
    qcond = qbdd || qlim;
    if(qcond) goto S140;
    step = stpmul*step;
    xlb = xub;
    xub = fifdmin1(xlb+step,big);
S140:
    goto S110;
S150:
    if(!(qlim && !qbdd)) goto S160;
    *status = -1;
    *qleft = 0;
    *qhi = !qincr;
    *x = big;
    return;
S160:
    goto S240;
S170:
//
//     HANDLE CASE IN WHICH WE MUST STEP LOWER
//
    xub = xsave;
    xlb = fifdmax1(xub-step,small);
    goto S190;
S180:
    if(qcond) goto S220;
S190:
//
//      YY = F(XLB) - Y
//
    *x = xlb;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 5;
    goto S300;
S200:
    yy = *fx;
    qbdd = qincr && yy <= 0.0e0 || !qincr && yy >= 0.0e0;
    qlim = xlb <= small;
    qcond = qbdd || qlim;
    if(qcond) goto S210;
    step = stpmul*step;
    xub = xlb;
    xlb = fifdmax1(xub-step,small);
S210:
    goto S180;
S220:
    if(!(qlim && !qbdd)) goto S230;
    *status = -1;
    *qleft = 1;
    *qhi = qincr;
    *x = small;
    return;
S240:
S230:
    dstzr(&xlb,&xub,&abstol,&reltol);
//
//  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
//
    *status = 0;
    goto S260;
S250:
    if(!(*status == 1)) goto S290;
S260:
    dzror ( status, x, fx, &xlo, &xhi, &qdum1, &qdum2 );
    if(!(*status == 1)) goto S280;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 6;
    goto S300;
S280:
S270:
    goto S250;
S290:
    *x = xlo;
    *status = 0;
    return;
DSTINV:
    small = *zsmall;
    big = *zbig;
    absstp = *zabsst;
    relstp = *zrelst;
    stpmul = *zstpmu;
    abstol = *zabsto;
    reltol = *zrelto;
    return;
S300:
//
//     TO GET-FUNCTION-VALUE
//
    *status = 1;
    return;
S310:
    switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S90;case
      4: goto S130;case 5: goto S200;case 6: goto S270;default: break;}
# undef qxmon
}
//****************************************************************************80

void E0001 ( int IENTRY, int *status, double *x, double *fx,
  double *xlo, double *xhi, unsigned long *qleft,
  unsigned long *qhi, double *zabstl, double *zreltl,
  double *zxhi, double *zxlo )

//****************************************************************************80
//
//  Purpose:
//
//    E00001 is a reverse-communication zero finder.
//
{
# define ftol(zx) (0.5e0*fifdmax1(abstol,reltol*fabs((zx))))

  static double a,abstol,b,c,d,fa,fb,fc,fd,fda;
  static double fdb,m,mb,p,q,reltol,tol,w,xxhi,xxlo;
  static int ext,i99999;
  static unsigned long first,qrzero;
    switch(IENTRY){case 0: goto DZROR; case 1: goto DSTZR;}
DZROR:
    if(*status > 0) goto S280;
    *xlo = xxlo;
    *xhi = xxhi;
    b = *x = *xlo;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 1;
    goto S270;
S10:
    fb = *fx;
    *xlo = *xhi;
    a = *x = *xlo;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 2;
    goto S270;
S20:
//
//     Check that F(ZXLO) < 0 < F(ZXHI)  or
//                F(ZXLO) > 0 > F(ZXHI)
//
    if(!(fb < 0.0e0)) goto S40;
    if(!(*fx < 0.0e0)) goto S30;
    *status = -1;
    *qleft = *fx < fb;
    *qhi = 0;
    return;
S40:
S30:
    if(!(fb > 0.0e0)) goto S60;
    if(!(*fx > 0.0e0)) goto S50;
    *status = -1;
    *qleft = *fx > fb;
    *qhi = 1;
    return;
S60:
S50:
    fa = *fx;
    first = 1;
S70:
    c = a;
    fc = fa;
    ext = 0;
S80:
    if(!(fabs(fc) < fabs(fb))) goto S100;
    if(!(c != a)) goto S90;
    d = a;
    fd = fa;
S90:
    a = b;
    fa = fb;
    *xlo = c;
    b = *xlo;
    fb = fc;
    c = a;
    fc = fa;
S100:
    tol = ftol(*xlo);
    m = (c+b)*.5e0;
    mb = m-b;
    if(!(fabs(mb) > tol)) goto S240;
    if(!(ext > 3)) goto S110;
    w = mb;
    goto S190;
S110:
    tol = fifdsign(tol,mb);
    p = (b-a)*fb;
    if(!first) goto S120;
    q = fa-fb;
    first = 0;
    goto S130;
S120:
    fdb = (fd-fb)/(d-b);
    fda = (fd-fa)/(d-a);
    p = fda*p;
    q = fdb*fa-fda*fb;
S130:
    if(!(p < 0.0e0)) goto S140;
    p = -p;
    q = -q;
S140:
    if(ext == 3) p *= 2.0e0;
    if(!(p*1.0e0 == 0.0e0 || p <= q*tol)) goto S150;
    w = tol;
    goto S180;
S150:
    if(!(p < mb*q)) goto S160;
    w = p/q;
    goto S170;
S160:
    w = mb;
S190:
S180:
S170:
    d = a;
    fd = fa;
    a = b;
    fa = fb;
    b += w;
    *xlo = b;
    *x = *xlo;
//
//     GET-FUNCTION-VALUE
//
    i99999 = 3;
    goto S270;
S200:
    fb = *fx;
    if(!(fc*fb >= 0.0e0)) goto S210;
    goto S70;
S210:
    if(!(w == mb)) goto S220;
    ext = 0;
    goto S230;
S220:
    ext += 1;
S230:
    goto S80;
S240:
    *xhi = c;
    qrzero = fc >= 0.0e0 && fb <= 0.0e0 || fc < 0.0e0 && fb >= 0.0e0;
    if(!qrzero) goto S250;
    *status = 0;
    goto S260;
S250:
    *status = -1;
S260:
    return;
DSTZR:
    xxlo = *zxlo;
    xxhi = *zxhi;
    abstol = *zabstl;
    reltol = *zreltl;
    return;
S270:
//
//     TO GET-FUNCTION-VALUE
//
    *status = 1;
    return;
S280:
    switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S200;
      default: break;}
# undef ftol
}
//****************************************************************************80

void erf_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    ERF_VALUES returns some values of the ERF or "error" function.
//
//  Definition:
//
//    ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double fx_vec[N_MAX] = {
    0.0000000000E+00, 0.1124629160E+00, 0.2227025892E+00, 0.3286267595E+00,
    0.4283923550E+00, 0.5204998778E+00, 0.6038560908E+00, 0.6778011938E+00,
    0.7421009647E+00, 0.7969082124E+00, 0.8427007929E+00, 0.8802050696E+00,
    0.9103139782E+00, 0.9340079449E+00, 0.9522851198E+00, 0.9661051465E+00,
    0.9763483833E+00, 0.9837904586E+00, 0.9890905016E+00, 0.9927904292E+00,
    0.9953222650E+00 };
  double x_vec[N_MAX] = {
    0.0E+00, 0.1E+00, 0.2E+00, 0.3E+00,
    0.4E+00, 0.5E+00, 0.6E+00, 0.7E+00,
    0.8E+00, 0.9E+00, 1.0E+00, 1.1E+00,
    1.2E+00, 1.3E+00, 1.4E+00, 1.5E+00,
    1.6E+00, 1.7E+00, 1.8E+00, 1.9E+00,
    2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

double error_f ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    ERROR_F evaluates the error function ERF.
//
//  Parameters:
//
//    Input, double *X, the argument.
//
//    Output, double ERROR_F, the value of the error function at X.
//
{
  static double c = .564189583547756e0;
  static double a[5] = {
    .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
    .479137145607681e-01,.128379167095513e+00
  };
  static double b[3] = {
    .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
  };
  static double p[8] = {
    -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
    4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
    4.51918953711873e+02,3.00459261020162e+02
  };
  static double q[8] = {
    1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
    2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
    7.90950925327898e+02,3.00459260956983e+02
  };
  static double r[5] = {
    2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
    4.65807828718470e+00,2.82094791773523e-01
  };
  static double s[4] = {
    9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
    1.80124575948747e+01
  };
  static double erf1,ax,bot,t,top,x2;

    ax = fabs(*x);
    if(ax > 0.5e0) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
    erf1 = *x*(top/bot);
    return erf1;
S10:
    if(ax > 4.0e0) goto S20;
    top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
      7];
    bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
      7];
    erf1 = 0.5e0+(0.5e0-exp(-(*x**x))*top/bot);
    if(*x < 0.0e0) erf1 = -erf1;
    return erf1;
S20:
    if(ax >= 5.8e0) goto S30;
    x2 = *x**x;
    t = 1.0e0/x2;
    top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
    bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
    erf1 = (c-top/(x2*bot))/ax;
    erf1 = 0.5e0+(0.5e0-exp(-x2)*erf1);
    if(*x < 0.0e0) erf1 = -erf1;
    return erf1;
S30:
    erf1 = fifdsign(1.0e0,*x);
    return erf1;
}
//****************************************************************************80

double error_fc ( int *ind, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    ERROR_FC evaluates the complementary error function ERFC.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, int *IND, chooses the scaling.
//    If IND is nonzero, then the value returned has been multiplied by
//    EXP(X*X).
//
//    Input, double *X, the argument of the function.
//
//    Output, double ERROR_FC, the value of the complementary
//    error function.
//
{
  static double c = .564189583547756e0;
  static double a[5] = {
    .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
    .479137145607681e-01,.128379167095513e+00
  };
  static double b[3] = {
    .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
  };
  static double p[8] = {
    -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
    4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
    4.51918953711873e+02,3.00459261020162e+02
  };
  static double q[8] = {
    1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
    2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
    7.90950925327898e+02,3.00459260956983e+02
  };
  static double r[5] = {
    2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
    4.65807828718470e+00,2.82094791773523e-01
  };
  static double s[4] = {
    9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
    1.80124575948747e+01
  };
  static int K1 = 1;
  static double erfc1,ax,bot,e,t,top,w;

//
//                     ABS(X) .LE. 0.5
//
    ax = fabs(*x);
    if(ax > 0.5e0) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
    erfc1 = 0.5e0+(0.5e0-*x*(top/bot));
    if(*ind != 0) erfc1 = exp(t)*erfc1;
    return erfc1;
S10:
//
//                  0.5 .LT. ABS(X) .LE. 4
//
    if(ax > 4.0e0) goto S20;
    top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
      7];
    bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
      7];
    erfc1 = top/bot;
    goto S40;
S20:
//
//                      ABS(X) .GT. 4
//
    if(*x <= -5.6e0) goto S60;
    if(*ind != 0) goto S30;
    if(*x > 100.0e0) goto S70;
    if(*x**x > -exparg(&K1)) goto S70;
S30:
    t = pow(1.0e0/ *x,2.0);
    top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
    bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
    erfc1 = (c-t*top/bot)/ax;
S40:
//
//                      FINAL ASSEMBLY
//
    if(*ind == 0) goto S50;
    if(*x < 0.0e0) erfc1 = 2.0e0*exp(*x**x)-erfc1;
    return erfc1;
S50:
    w = *x**x;
    t = w;
    e = w-t;
    erfc1 = (0.5e0+(0.5e0-e))*exp(-t)*erfc1;
    if(*x < 0.0e0) erfc1 = 2.0e0-erfc1;
    return erfc1;
S60:
//
//             LIMIT VALUE FOR LARGE NEGATIVE X
//
    erfc1 = 2.0e0;
    if(*ind != 0) erfc1 = 2.0e0*exp(*x**x);
    return erfc1;
S70:
//
//             LIMIT VALUE FOR LARGE POSITIVE X
//                       WHEN IND = 0
//
    erfc1 = 0.0e0;
    return erfc1;
}
//****************************************************************************80

double esum ( int *mu, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    ESUM evaluates exp ( MU + X ).
//
//  Parameters:
//
//    Input, int *MU, part of the argument.
//
//    Input, double *X, part of the argument.
//
//    Output, double ESUM, the value of exp ( MU + X ).
//
{
  static double esum,w;

    if(*x > 0.0e0) goto S10;
    if(*mu < 0) goto S20;
    w = (double)*mu+*x;
    if(w > 0.0e0) goto S20;
    esum = exp(w);
    return esum;
S10:
    if(*mu > 0) goto S20;
    w = (double)*mu+*x;
    if(w < 0.0e0) goto S20;
    esum = exp(w);
    return esum;
S20:
    w = *mu;
    esum = exp(w)*exp(*x);
    return esum;
}
//****************************************************************************80

double eval_pol ( double a[], int *n, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    EVAL_POL evaluates a polynomial at X.
//
//  Discussion:
//
//    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
//
//  Modified:
//
//    15 December 1999
//
//  Parameters:
//
//    Input, double precision A(0:N), coefficients of the polynomial.
//
//    Input, int *N, length of A.
//
//    Input, double *X, the point at which the polynomial
//    is to be evaluated.
//
//    Output, double EVAL_POL, the value of the polynomial at X.
//
{
  static double devlpl,term;
  static int i;

  term = a[*n-1];
  for ( i = *n-1-1; i >= 0; i-- )
  {
    term = a[i]+term**x;
  }

  devlpl = term;
  return devlpl;
}
//****************************************************************************80

double exparg ( int *l )

//****************************************************************************80
//
//  Purpose:
//
//    EXPARG returns the largest or smallest legal argument for EXP.
//
//  Discussion:
//
//    Only an approximate limit for the argument of EXP is desired.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, int *L, indicates which limit is desired.
//    If L = 0, then the largest positive argument for EXP is desired.
//    Otherwise, the largest negative argument for EXP for which the
//    result is nonzero is desired.
//
//    Output, double EXPARG, the desired value.
//
{
  static int K1 = 4;
  static int K2 = 9;
  static int K3 = 10;
  static double exparg,lnb;
  static int b,m;

    b = ipmpar(&K1);
    if(b != 2) goto S10;
    lnb = .69314718055995e0;
    goto S40;
S10:
    if(b != 8) goto S20;
    lnb = 2.0794415416798e0;
    goto S40;
S20:
    if(b != 16) goto S30;
    lnb = 2.7725887222398e0;
    goto S40;
S30:
    lnb = log((double)b);
S40:
    if(*l == 0) goto S50;
    m = ipmpar(&K2)-1;
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
S50:
    m = ipmpar(&K3);
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
}
//****************************************************************************80

void f_cdf_values ( int *n_data, int *a, int *b, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    F_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by
//    commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[FRatioDistribution[ DFN, DFD ], X ]
//
//  Modified:
//
//    11 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *A, int *B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 20

  int a_vec[N_MAX] = {
    1, 1, 5, 1,
    2, 4, 1, 6,
    8, 1, 3, 6,
    1, 1, 1, 1,
    2, 3, 4, 5 };
  int b_vec[N_MAX] = {
     1,  5,  1,  5,
    10, 20,  5,  6,
    16,  5, 10, 12,
     5,  5,  5,  5,
     5,  5,  5,  5 };
  double fx_vec[N_MAX] = {
    0.500000E+00, 0.499971E+00, 0.499603E+00, 0.749699E+00,
    0.750466E+00, 0.751416E+00, 0.899987E+00, 0.899713E+00,
    0.900285E+00, 0.950025E+00, 0.950057E+00, 0.950193E+00,
    0.975013E+00, 0.990002E+00, 0.994998E+00, 0.999000E+00,
    0.568799E+00, 0.535145E+00, 0.514343E+00, 0.500000E+00 };
  double x_vec[N_MAX] = {
    1.00E+00,  0.528E+00, 1.89E+00,  1.69E+00,
    1.60E+00,  1.47E+00,  4.06E+00,  3.05E+00,
    2.09E+00,  6.61E+00,  3.71E+00,  3.00E+00,
   10.01E+00, 16.26E+00, 22.78E+00, 47.18E+00,
    1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

void f_noncentral_cdf_values ( int *n_data, int *a, int *b, double *lambda,
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
//    in Mathematica by commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *A, int *B, double *LAMBDA, the
//    parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 22

  int a_vec[N_MAX] = {
     1,  1,  1,  1,
     1,  1,  1,  1,
     1,  1,  2,  2,
     3,  3,  4,  4,
     5,  5,  6,  6,
     8, 16 };
  int b_vec[N_MAX] = {
     1,  5,  5,  5,
     5,  5,  5,  5,
     5,  5,  5, 10,
     5,  5,  5,  5,
     1,  5,  6, 12,
    16,  8 };
  double fx_vec[N_MAX] = {
    0.500000E+00, 0.636783E+00, 0.584092E+00, 0.323443E+00,
    0.450119E+00, 0.607888E+00, 0.705928E+00, 0.772178E+00,
    0.819105E+00, 0.317035E+00, 0.432722E+00, 0.450270E+00,
    0.426188E+00, 0.337744E+00, 0.422911E+00, 0.692767E+00,
    0.363217E+00, 0.421005E+00, 0.426667E+00, 0.446402E+00,
    0.844589E+00, 0.816368E+00 };
  double lambda_vec[N_MAX] = {
    0.00E+00,  0.000E+00, 0.25E+00,  1.00E+00,
    1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
    0.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  1.00E+00 };
  double x_vec[N_MAX] = {
    1.00E+00,  1.00E+00, 1.00E+00,  0.50E+00,
    1.00E+00,  2.00E+00, 3.00E+00,  4.00E+00,
    5.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
    1.00E+00,  1.00E+00, 1.00E+00,  2.00E+00,
    1.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
    2.00E+00,  2.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0;
    *lambda = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double fifdint ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    FIFDINT truncates a double number to an integer.
//
//  Parameters:
//
// a     -     number to be truncated
{
  return (double) ((int) a);
}
//****************************************************************************80

double fifdmax1 ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FIFDMAX1 returns the maximum of two numbers a and b
//
//  Parameters:
//
//  a     -      first number
//  b     -      second number
//
{
  if ( a < b )
  {
    return b;
  }
  else
  {
    return a;
  }
}
//****************************************************************************80

double fifdmin1 ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FIFDMIN1 returns the minimum of two numbers.
//
//  Parameters:
//
//  a     -     first number
//  b     -     second number
//
{
  if (a < b) return a;
  else return b;
}
//****************************************************************************80

double fifdsign ( double mag, double sign )

//****************************************************************************80
//
//  Purpose:
//
//    FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
//
//  Parameters:
//
//  mag     -     magnitude
//  sign    -     sign to be transfered
//
{
  if (mag < 0) mag = -mag;
  if (sign < 0) mag = -mag;
  return mag;

}
//****************************************************************************80

long fifidint ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    FIFIDINT truncates a double number to a long integer
//
//  Parameters:
//
//  a - number to be truncated
//
{
  if ( a < 1.0 )
  {
    return (long) 0;
  }
  else
  {
    return ( long ) a;
  }
}
//****************************************************************************80

long fifmod ( long a, long b )

//****************************************************************************80
//
//  Purpose:
//
//    FIFMOD returns the modulo of a and b
//
//  Parameters:
//
//  a - numerator
//  b - denominator
//
{
  return ( a % b );
}
//****************************************************************************80

double fpser ( double *a, double *b, double *x, double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    FPSER evaluates IX(A,B)(X) for very small B.
//
//  Discussion:
//
//    This routine is appropriate for use when
//
//      B < min ( EPS, EPS * A )
//
//    and
//
//      X <= 0.5.
//
//  Parameters:
//
//    Input, double *A, *B, parameters of the function.
//
//    Input, double *X, the point at which the function is to
//    be evaluated.
//
//    Input, double *EPS, a tolerance.
//
//    Output, double FPSER, the value of IX(A,B)(X).
//
{
  static int K1 = 1;
  static double fpser,an,c,s,t,tol;

    fpser = 1.0e0;
    if(*a <= 1.e-3**eps) goto S10;
    fpser = 0.0e0;
    t = *a*log(*x);
    if(t < exparg(&K1)) return fpser;
    fpser = exp(t);
S10:
//
//                NOTE THAT 1/B(A,B) = B
//
    fpser = *b/ *a*fpser;
    tol = *eps/ *a;
    an = *a+1.0e0;
    t = *x;
    s = t/an;
S20:
    an += 1.0e0;
    t = *x*t;
    c = t/an;
    s += c;
    if(fabs(c) > tol) goto S20;
    fpser *= (1.0e0+*a*s);
    return fpser;
}
//****************************************************************************80

void ftnstop ( string msg )

//****************************************************************************80
//
//  Purpose:
//
//    FTNSTOP prints a message to standard error and then exits.
//
//  Parameters:
//
//    Input, string MSG, the message to be printed.
//
{
  //cerr << msg << "\n";

  exit ( 0 );
}
//****************************************************************************80

double gam1 ( double *a )

//****************************************************************************80
//
//  Purpose:
//
//    GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5D+00 <= A <= 1.5
//
//  Parameters:
//
//    Input, double *A, forms the argument of the Gamma function.
//
//    Output, double GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
//
{
  static double s1 = .273076135303957e+00;
  static double s2 = .559398236957378e-01;
  static double p[7] = {
    .577215664901533e+00,-.409078193005776e+00,-.230975380857675e+00,
    .597275330452234e-01,.766968181649490e-02,-.514889771323592e-02,
    .589597428611429e-03
  };
  static double q[5] = {
    .100000000000000e+01,.427569613095214e+00,.158451672430138e+00,
    .261132021441447e-01,.423244297896961e-02
  };
  static double r[9] = {
    -.422784335098468e+00,-.771330383816272e+00,-.244757765222226e+00,
    .118378989872749e+00,.930357293360349e-03,-.118290993445146e-01,
    .223047661158249e-02,.266505979058923e-03,-.132674909766242e-03
  };
  static double gam1,bot,d,t,top,w,T1;

    t = *a;
    d = *a-0.5e0;
    if(d > 0.0e0) t = d-0.5e0;
    T1 = t;
    if(T1 < 0) goto S40;
    else if(T1 == 0) goto S10;
    else  goto S20;
S10:
    gam1 = 0.0e0;
    return gam1;
S20:
    top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
    bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1.0e0;
    w = top/bot;
    if(d > 0.0e0) goto S30;
    gam1 = *a*w;
    return gam1;
S30:
    gam1 = t/ *a*(w-0.5e0-0.5e0);
    return gam1;
S40:
    top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+r[1])*t+
      r[0];
    bot = (s2*t+s1)*t+1.0e0;
    w = top/bot;
    if(d > 0.0e0) goto S50;
    gam1 = *a*(w+0.5e0+0.5e0);
    return gam1;
S50:
    gam1 = t*w/ *a;
    return gam1;
}
//****************************************************************************80

void gamma_inc ( double *a, double *x, double *ans, double *qans, int *ind )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
//
//  Discussion:
//
//    This is certified spaghetti code.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double *A, *X, the arguments of the incomplete
//    gamma ratio.  A and X must be nonnegative.  A and X cannot
//    both be zero.
//
//    Output, double *ANS, *QANS.  On normal output,
//    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
//    A or X is negative, or both are 0, or when the answer is
//    computationally indeterminate because A is extremely large
//    and X is very close to A.
//
//    Input, int *IND, indicates the accuracy request:
//    0, as much accuracy as possible.
//    1, to within 1 unit of the 6-th significant digit,
//    otherwise, to within 1 unit of the 3rd significant digit.
//
{
  static double alog10 = 2.30258509299405e0;
  static double d10 = -.185185185185185e-02;
  static double d20 = .413359788359788e-02;
  static double d30 = .649434156378601e-03;
  static double d40 = -.861888290916712e-03;
  static double d50 = -.336798553366358e-03;
  static double d60 = .531307936463992e-03;
  static double d70 = .344367606892378e-03;
  static double rt2pin = .398942280401433e0;
  static double rtpi = 1.77245385090552e0;
  static double third = .333333333333333e0;
  static double acc0[3] = {
    5.e-15,5.e-7,5.e-4
  };
  static double big[3] = {
    20.0e0,14.0e0,10.0e0
  };
  static double d0[13] = {
    .833333333333333e-01,-.148148148148148e-01,.115740740740741e-02,
    .352733686067019e-03,-.178755144032922e-03,.391926317852244e-04,
    -.218544851067999e-05,-.185406221071516e-05,.829671134095309e-06,
    -.176659527368261e-06,.670785354340150e-08,.102618097842403e-07,
    -.438203601845335e-08
  };
  static double d1[12] = {
    -.347222222222222e-02,.264550264550265e-02,-.990226337448560e-03,
    .205761316872428e-03,-.401877572016461e-06,-.180985503344900e-04,
    .764916091608111e-05,-.161209008945634e-05,.464712780280743e-08,
    .137863344691572e-06,-.575254560351770e-07,.119516285997781e-07
  };
  static double d2[10] = {
    -.268132716049383e-02,.771604938271605e-03,.200938786008230e-05,
    -.107366532263652e-03,.529234488291201e-04,-.127606351886187e-04,
    .342357873409614e-07,.137219573090629e-05,-.629899213838006e-06,
    .142806142060642e-06
  };
  static double d3[8] = {
    .229472093621399e-03,-.469189494395256e-03,.267720632062839e-03,
    -.756180167188398e-04,-.239650511386730e-06,.110826541153473e-04,
    -.567495282699160e-05,.142309007324359e-05
  };
  static double d4[6] = {
    .784039221720067e-03,-.299072480303190e-03,-.146384525788434e-05,
    .664149821546512e-04,-.396836504717943e-04,.113757269706784e-04
  };
  static double d5[4] = {
    -.697281375836586e-04,.277275324495939e-03,-.199325705161888e-03,
    .679778047793721e-04
  };
  static double d6[2] = {
    -.592166437353694e-03,.270878209671804e-03
  };
  static double e00[3] = {
    .25e-3,.25e-1,.14e0
  };
  static double x00[3] = {
    31.0e0,17.0e0,9.7e0
  };
  static int K1 = 1;
  static int K2 = 0;
  static double a2n,a2nm1,acc,am0,amn,an,an0,apn,b2n,b2nm1,c,c0,c1,c2,c3,c4,c5,c6,
    cma,e,e0,g,h,j,l,r,rta,rtx,s,sum,t,t1,tol,twoa,u,w,x0,y,z;
  static int i,iop,m,max,n;
  static double wk[20],T3;
  static int T4,T5;
  static double T6,T7;

//
//  E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
//  NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
//
    e = dpmpar(&K1);
    if(*a < 0.0e0 || *x < 0.0e0) goto S430;
    if(*a == 0.0e0 && *x == 0.0e0) goto S430;
    if(*a**x == 0.0e0) goto S420;
    iop = *ind+1;
    if(iop != 1 && iop != 2) iop = 3;
    acc = fifdmax1(acc0[iop-1],e);
    e0 = e00[iop-1];
    x0 = x00[iop-1];
//
//  SELECT THE APPROPRIATE ALGORITHM
//
    if(*a >= 1.0e0) goto S10;
    if(*a == 0.5e0) goto S390;
    if(*x < 1.1e0) goto S160;
    t1 = *a*log(*x)-*x;
    u = *a*exp(t1);
    if(u == 0.0e0) goto S380;
    r = u*(1.0e0+gam1(a));
    goto S250;
S10:
    if(*a >= big[iop-1]) goto S30;
    if(*a > *x || *x >= x0) goto S20;
    twoa = *a+*a;
    m = fifidint(twoa);
    if(twoa != (double)m) goto S20;
    i = m/2;
    if(*a == (double)i) goto S210;
    goto S220;
S20:
    t1 = *a*log(*x)-*x;
    r = exp(t1)/ gamma_x(a);
    goto S40;
S30:
    l = *x/ *a;
    if(l == 0.0e0) goto S370;
    s = 0.5e0+(0.5e0-l);
    z = rlog(&l);
    if(z >= 700.0e0/ *a) goto S410;
    y = *a*z;
    rta = sqrt(*a);
    if(fabs(s) <= e0/rta) goto S330;
    if(fabs(s) <= 0.4e0) goto S270;
    t = pow(1.0e0/ *a,2.0);
    t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
    t1 -= y;
    r = rt2pin*rta*exp(t1);
S40:
    if(r == 0.0e0) goto S420;
    if(*x <= fifdmax1(*a,alog10)) goto S50;
    if(*x < x0) goto S250;
    goto S100;
S50:
//
//  TAYLOR SERIES FOR P/R
//
    apn = *a+1.0e0;
    t = *x/apn;
    wk[0] = t;
    for ( n = 2; n <= 20; n++ )
    {
        apn += 1.0e0;
        t *= (*x/apn);
        if(t <= 1.e-3) goto S70;
        wk[n-1] = t;
    }
    n = 20;
S70:
    sum = t;
    tol = 0.5e0*acc;
S80:
    apn += 1.0e0;
    t *= (*x/apn);
    sum += t;
    if(t > tol) goto S80;
    max = n-1;
    for ( m = 1; m <= max; m++ )
    {
        n -= 1;
        sum += wk[n-1];
    }
    *ans = r/ *a*(1.0e0+sum);
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S100:
//
//  ASYMPTOTIC EXPANSION
//
    amn = *a-1.0e0;
    t = amn/ *x;
    wk[0] = t;
    for ( n = 2; n <= 20; n++ )
    {
        amn -= 1.0e0;
        t *= (amn/ *x);
        if(fabs(t) <= 1.e-3) goto S120;
        wk[n-1] = t;
    }
    n = 20;
S120:
    sum = t;
S130:
    if(fabs(t) <= acc) goto S140;
    amn -= 1.0e0;
    t *= (amn/ *x);
    sum += t;
    goto S130;
S140:
    max = n-1;
    for ( m = 1; m <= max; m++ )
    {
        n -= 1;
        sum += wk[n-1];
    }
    *qans = r/ *x*(1.0e0+sum);
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S160:
//
//  TAYLOR SERIES FOR P(A,X)/X**A
//
    an = 3.0e0;
    c = *x;
    sum = *x/(*a+3.0e0);
    tol = 3.0e0*acc/(*a+1.0e0);
S170:
    an += 1.0e0;
    c = -(c*(*x/an));
    t = c/(*a+an);
    sum += t;
    if(fabs(t) > tol) goto S170;
    j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
    z = *a*log(*x);
    h = gam1(a);
    g = 1.0e0+h;
    if(*x < 0.25e0) goto S180;
    if(*a < *x/2.59e0) goto S200;
    goto S190;
S180:
    if(z > -.13394e0) goto S200;
S190:
    w = exp(z);
    *ans = w*g*(0.5e0+(0.5e0-j));
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S200:
    l = rexp(&z);
    w = 0.5e0+(0.5e0+l);
    *qans = (w*j-l)*g-h;
    if(*qans < 0.0e0) goto S380;
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S210:
//
//  FINITE SUMS FOR Q WHEN A .GE. 1 AND 2*A IS AN INTEGER
//
    sum = exp(-*x);
    t = sum;
    n = 1;
    c = 0.0e0;
    goto S230;
S220:
    rtx = sqrt(*x);
    sum = error_fc ( &K2, &rtx );
    t = exp(-*x)/(rtpi*rtx);
    n = 0;
    c = -0.5e0;
S230:
    if(n == i) goto S240;
    n += 1;
    c += 1.0e0;
    t = *x*t/c;
    sum += t;
    goto S230;
S240:
    *qans = sum;
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S250:
//
//  CONTINUED FRACTION EXPANSION
//
    tol = fifdmax1(5.0e0*e,acc);
    a2nm1 = a2n = 1.0e0;
    b2nm1 = *x;
    b2n = *x+(1.0e0-*a);
    c = 1.0e0;
S260:
    a2nm1 = *x*a2n+c*a2nm1;
    b2nm1 = *x*b2n+c*b2nm1;
    am0 = a2nm1/b2nm1;
    c += 1.0e0;
    cma = c-*a;
    a2n = a2nm1+cma*a2n;
    b2n = b2nm1+cma*b2n;
    an0 = a2n/b2n;
    if(fabs(an0-am0) >= tol*an0) goto S260;
    *qans = r*an0;
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S270:
//
//  GENERAL TEMME EXPANSION
//
    if(fabs(s) <= 2.0e0*e && *a*e*e > 3.28e-3) goto S430;
    c = exp(-y);
    T3 = sqrt(y);
    w = 0.5e0 * error_fc ( &K1, &T3 );
    u = 1.0e0/ *a;
    z = sqrt(z+z);
    if(l < 1.0e0) z = -z;
    T4 = iop-2;
    if(T4 < 0) goto S280;
    else if(T4 == 0) goto S290;
    else  goto S300;
S280:
    if(fabs(s) <= 1.e-3) goto S340;
    c0 = ((((((((((((d0[12]*z+d0[11])*z+d0[10])*z+d0[9])*z+d0[8])*z+d0[7])*z+d0[
      6])*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
    c1 = (((((((((((d1[11]*z+d1[10])*z+d1[9])*z+d1[8])*z+d1[7])*z+d1[6])*z+d1[5]
      )*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
    c2 = (((((((((d2[9]*z+d2[8])*z+d2[7])*z+d2[6])*z+d2[5])*z+d2[4])*z+d2[3])*z+
      d2[2])*z+d2[1])*z+d2[0])*z+d20;
    c3 = (((((((d3[7]*z+d3[6])*z+d3[5])*z+d3[4])*z+d3[3])*z+d3[2])*z+d3[1])*z+
      d3[0])*z+d30;
    c4 = (((((d4[5]*z+d4[4])*z+d4[3])*z+d4[2])*z+d4[1])*z+d4[0])*z+d40;
    c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z+d50;
    c6 = (d6[1]*z+d6[0])*z+d60;
    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
    goto S310;
S290:
    c0 = (((((d0[5]*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
    c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
    c2 = d2[0]*z+d20;
    t = (c2*u+c1)*u+c0;
    goto S310;
S300:
    t = ((d0[2]*z+d0[1])*z+d0[0])*z-third;
S310:
    if(l < 1.0e0) goto S320;
    *qans = c*(w+rt2pin*t/rta);
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S320:
    *ans = c*(w-rt2pin*t/rta);
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S330:
//
//  TEMME EXPANSION FOR L = 1
//
    if(*a*e*e > 3.28e-3) goto S430;
    c = 0.5e0+(0.5e0-y);
    w = (0.5e0-sqrt(y)*(0.5e0+(0.5e0-y/3.0e0))/rtpi)/c;
    u = 1.0e0/ *a;
    z = sqrt(z+z);
    if(l < 1.0e0) z = -z;
    T5 = iop-2;
    if(T5 < 0) goto S340;
    else if(T5 == 0) goto S350;
    else  goto S360;
S340:
    c0 = ((((((d0[6]*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-
      third;
    c1 = (((((d1[5]*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
    c2 = ((((d2[4]*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z+d20;
    c3 = (((d3[3]*z+d3[2])*z+d3[1])*z+d3[0])*z+d30;
    c4 = (d4[1]*z+d4[0])*z+d40;
    c5 = (d5[1]*z+d5[0])*z+d50;
    c6 = d6[0]*z+d60;
    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
    goto S310;
S350:
    c0 = (d0[1]*z+d0[0])*z-third;
    c1 = d1[0]*z+d10;
    t = (d20*u+c1)*u+c0;
    goto S310;
S360:
    t = d0[0]*z-third;
    goto S310;
S370:
//
//  SPECIAL CASES
//
    *ans = 0.0e0;
    *qans = 1.0e0;
    return;
S380:
    *ans = 1.0e0;
    *qans = 0.0e0;
    return;
S390:
    if(*x >= 0.25e0) goto S400;
    T6 = sqrt(*x);
    *ans = error_f ( &T6 );
    *qans = 0.5e0+(0.5e0-*ans);
    return;
S400:
    T7 = sqrt(*x);
    *qans = error_fc ( &K2, &T7 );
    *ans = 0.5e0+(0.5e0-*qans);
    return;
S410:
    if(fabs(s) <= 2.0e0*e) goto S430;
S420:
    if(*x <= *a) goto S370;
    goto S380;
S430:
//
//  ERROR RETURN
//
    *ans = 2.0e0;
    return;
}
//****************************************************************************80

void gamma_inc_inv ( double *a, double *x, double *x0, double *p, double *q,
  int *ierr )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
//
//  Discussion:
//
//    The routine is given positive A, and nonnegative P and Q where P + Q = 1.
//    The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.
//    Schroder iteration is employed.  The routine attempts to compute X
//    to 10 significant digits if this is possible for the particular computer
//    arithmetic being used.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double *A, the parameter in the incomplete gamma
//    ratio.  A must be positive.
//
//    Output, double *X, the computed point for which the
//    incomplete gamma functions have the values P and Q.
//
//    Input, double *X0, an optional initial approximation
//    for the solution X.  If the user does not want to supply an
//    initial approximation, then X0 should be set to 0, or a negative
//    value.
//
//    Input, double *P, *Q, the values of the incomplete gamma
//    functions, for which the corresponding argument is desired.
//
//    Output, int *IERR, error flag.
//    0, the solution was obtained. Iteration was not used.
//    0 < K, The solution was obtained. IERR iterations were performed.
//    -2, A <= 0
//    -3, No solution was obtained. The ratio Q/A is too large.
//    -4, P + Q /= 1
//    -6, 20 iterations were performed. The most recent value obtained
//        for X is given.  This cannot occur if X0 <= 0.
//    -7, Iteration failed. No value is given for X.
//        This may occur when X is approximately 0.
//    -8, A value for X has been obtained, but the routine is not certain
//        of its accuracy.  Iteration cannot be performed in this
//        case. If X0 <= 0, this can occur only when P or Q is
//        approximately 0. If X0 is positive then this can occur when A is
//        exceedingly close to X and A is extremely large (say A .GE. 1.E20).
//
{
  static double a0 = 3.31125922108741e0;
  static double a1 = 11.6616720288968e0;
  static double a2 = 4.28342155967104e0;
  static double a3 = .213623493715853e0;
  static double b1 = 6.61053765625462e0;
  static double b2 = 6.40691597760039e0;
  static double b3 = 1.27364489782223e0;
  static double b4 = .036117081018842e0;
  static double c = .577215664901533e0;
  static double ln10 = 2.302585e0;
  static double tol = 1.e-5;
  static double amin[2] = {
    500.0e0,100.0e0
  };
  static double bmin[2] = {
    1.e-28,1.e-13
  };
  static double dmin[2] = {
    1.e-06,1.e-04
  };
  static double emin[2] = {
    2.e-03,6.e-03
  };
  static double eps0[2] = {
    1.e-10,1.e-08
  };
  static int K1 = 1;
  static int K2 = 2;
  static int K3 = 3;
  static int K8 = 0;
  static double am1,amax,ap1,ap2,ap3,apn,b,c1,c2,c3,c4,c5,d,e,e2,eps,g,h,pn,qg,qn,
    r,rta,s,s2,sum,t,u,w,xmax,xmin,xn,y,z;
  static int iop;
  static double T4,T5,T6,T7,T9;

//
//  E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
//            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
//            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
//            LARGEST POSITIVE NUMBER.
//
    e = dpmpar(&K1);
    xmin = dpmpar(&K2);
    xmax = dpmpar(&K3);
    *x = 0.0e0;
    if(*a <= 0.0e0) goto S300;
    t = *p+*q-1.e0;
    if(fabs(t) > e) goto S320;
    *ierr = 0;
    if(*p == 0.0e0) return;
    if(*q == 0.0e0) goto S270;
    if(*a == 1.0e0) goto S280;
    e2 = 2.0e0*e;
    amax = 0.4e-10/(e*e);
    iop = 1;
    if(e > 1.e-10) iop = 2;
    eps = eps0[iop-1];
    xn = *x0;
    if(*x0 > 0.0e0) goto S160;
//
//        SELECTION OF THE INITIAL APPROXIMATION XN OF X
//                       WHEN A .LT. 1
//
    if(*a > 1.0e0) goto S80;
    T4 = *a+1.0e0;
    g = gamma_x(&T4);
    qg = *q*g;
    if(qg == 0.0e0) goto S360;
    b = qg/ *a;
    if(qg > 0.6e0**a) goto S40;
    if(*a >= 0.30e0 || b < 0.35e0) goto S10;
    t = exp(-(b+c));
    u = t*exp(t);
    xn = t*exp(u);
    goto S160;
S10:
    if(b >= 0.45e0) goto S40;
    if(b == 0.0e0) goto S360;
    y = -log(b);
    s = 0.5e0+(0.5e0-*a);
    z = log(y);
    t = y-s*z;
    if(b < 0.15e0) goto S20;
    xn = y-s*log(t)-log(1.0e0+s/(t+1.0e0));
    goto S220;
S20:
    if(b <= 0.01e0) goto S30;
    u = ((t+2.0e0*(3.0e0-*a))*t+(2.0e0-*a)*(3.0e0-*a))/((t+(5.0e0-*a))*t+2.0e0);
    xn = y-s*log(t)-log(u);
    goto S220;
S30:
    c1 = -(s*z);
    c2 = -(s*(1.0e0+c1));
    c3 = s*((0.5e0*c1+(2.0e0-*a))*c1+(2.5e0-1.5e0**a));
    c4 = -(s*(((c1/3.0e0+(2.5e0-1.5e0**a))*c1+((*a-6.0e0)**a+7.0e0))*c1+(
      (11.0e0**a-46.0)**a+47.0e0)/6.0e0));
    c5 = -(s*((((-(c1/4.0e0)+(11.0e0**a-17.0e0)/6.0e0)*c1+((-(3.0e0**a)+13.0e0)*
      *a-13.0e0))*c1+0.5e0*(((2.0e0**a-25.0e0)**a+72.0e0)**a-61.0e0))*c1+((
      (25.0e0**a-195.0e0)**a+477.0e0)**a-379.0e0)/12.0e0));
    xn = (((c5/y+c4)/y+c3)/y+c2)/y+c1+y;
    if(*a > 1.0e0) goto S220;
    if(b > bmin[iop-1]) goto S220;
    *x = xn;
    return;
S40:
    if(b**q > 1.e-8) goto S50;
    xn = exp(-(*q/ *a+c));
    goto S70;
S50:
    if(*p <= 0.9e0) goto S60;
    T5 = -*q;
    xn = exp((alnrel(&T5)+ gamma_ln1 ( a ) ) / *a );
    goto S70;
S60:
    xn = exp(log(*p*g)/ *a);
S70:
    if(xn == 0.0e0) goto S310;
    t = 0.5e0+(0.5e0-xn/(*a+1.0e0));
    xn /= t;
    goto S160;
S80:
//
//        SELECTION OF THE INITIAL APPROXIMATION XN OF X
//                       WHEN A .GT. 1
//
    if(*q <= 0.5e0) goto S90;
    w = log(*p);
    goto S100;
S90:
    w = log(*q);
S100:
    t = sqrt(-(2.0e0*w));
    s = t-(((a3*t+a2)*t+a1)*t+a0)/((((b4*t+b3)*t+b2)*t+b1)*t+1.0e0);
    if(*q > 0.5e0) s = -s;
    rta = sqrt(*a);
    s2 = s*s;
    xn = *a+s*rta+(s2-1.0e0)/3.0e0+s*(s2-7.0e0)/(36.0e0*rta)-((3.0e0*s2+7.0e0)*
      s2-16.0e0)/(810.0e0**a)+s*((9.0e0*s2+256.0e0)*s2-433.0e0)/(38880.0e0**a*
      rta);
    xn = fifdmax1(xn,0.0e0);
    if(*a < amin[iop-1]) goto S110;
    *x = xn;
    d = 0.5e0+(0.5e0-*x/ *a);
    if(fabs(d) <= dmin[iop-1]) return;
S110:
    if(*p <= 0.5e0) goto S130;
    if(xn < 3.0e0**a) goto S220;
    y = -(w+ gamma_log ( a ) );
    d = fifdmax1(2.0e0,*a*(*a-1.0e0));
    if(y < ln10*d) goto S120;
    s = 1.0e0-*a;
    z = log(y);
    goto S30;
S120:
    t = *a-1.0e0;
    T6 = -(t/(xn+1.0e0));
    xn = y+t*log(xn)-alnrel(&T6);
    T7 = -(t/(xn+1.0e0));
    xn = y+t*log(xn)-alnrel(&T7);
    goto S220;
S130:
    ap1 = *a+1.0e0;
    if(xn > 0.70e0*ap1) goto S170;
    w += gamma_log ( &ap1 );
    if(xn > 0.15e0*ap1) goto S140;
    ap2 = *a+2.0e0;
    ap3 = *a+3.0e0;
    *x = exp((w+*x)/ *a);
    *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
    *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
    *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2*(1.0e0+*x/ap3))))/ *a);
    xn = *x;
    if(xn > 1.e-2*ap1) goto S140;
    if(xn <= emin[iop-1]*ap1) return;
    goto S170;
S140:
    apn = ap1;
    t = xn/apn;
    sum = 1.0e0+t;
S150:
    apn += 1.0e0;
    t *= (xn/apn);
    sum += t;
    if(t > 1.e-4) goto S150;
    t = w-log(sum);
    xn = exp((xn+t)/ *a);
    xn *= (1.0e0-(*a*log(xn)-xn-t)/(*a-xn));
    goto S170;
S160:
//
//                 SCHRODER ITERATION USING P
//
    if(*p > 0.5e0) goto S220;
S170:
    if(*p <= 1.e10*xmin) goto S350;
    am1 = *a-0.5e0-0.5e0;
S180:
    if(*a <= amax) goto S190;
    d = 0.5e0+(0.5e0-xn/ *a);
    if(fabs(d) <= e2) goto S350;
S190:
    if(*ierr >= 20) goto S330;
    *ierr += 1;
    gamma_inc ( a, &xn, &pn, &qn, &K8 );
    if(pn == 0.0e0 || qn == 0.0e0) goto S350;
    r = rcomp(a,&xn);
    if(r == 0.0e0) goto S350;
    t = (pn-*p)/r;
    w = 0.5e0*(am1-xn);
    if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S200;
    *x = xn*(1.0e0-t);
    if(*x <= 0.0e0) goto S340;
    d = fabs(t);
    goto S210;
S200:
    h = t*(1.0e0+w*t);
    *x = xn*(1.0e0-h);
    if(*x <= 0.0e0) goto S340;
    if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return;
    d = fabs(h);
S210:
    xn = *x;
    if(d > tol) goto S180;
    if(d <= eps) return;
    if(fabs(*p-pn) <= tol**p) return;
    goto S180;
S220:
//
//                 SCHRODER ITERATION USING Q
//
    if(*q <= 1.e10*xmin) goto S350;
    am1 = *a-0.5e0-0.5e0;
S230:
    if(*a <= amax) goto S240;
    d = 0.5e0+(0.5e0-xn/ *a);
    if(fabs(d) <= e2) goto S350;
S240:
    if(*ierr >= 20) goto S330;
    *ierr += 1;
    gamma_inc ( a, &xn, &pn, &qn, &K8 );
    if(pn == 0.0e0 || qn == 0.0e0) goto S350;
    r = rcomp(a,&xn);
    if(r == 0.0e0) goto S350;
    t = (*q-qn)/r;
    w = 0.5e0*(am1-xn);
    if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S250;
    *x = xn*(1.0e0-t);
    if(*x <= 0.0e0) goto S340;
    d = fabs(t);
    goto S260;
S250:
    h = t*(1.0e0+w*t);
    *x = xn*(1.0e0-h);
    if(*x <= 0.0e0) goto S340;
    if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return;
    d = fabs(h);
S260:
    xn = *x;
    if(d > tol) goto S230;
    if(d <= eps) return;
    if(fabs(*q-qn) <= tol**q) return;
    goto S230;
S270:
//
//                       SPECIAL CASES
//
    *x = xmax;
    return;
S280:
    if(*q < 0.9e0) goto S290;
    T9 = -*p;
    *x = -alnrel(&T9);
    return;
S290:
    *x = -log(*q);
    return;
S300:
//
//                       ERROR RETURN
//
    *ierr = -2;
    return;
S310:
    *ierr = -3;
    return;
S320:
    *ierr = -4;
    return;
S330:
    *ierr = -6;
    return;
S340:
    *ierr = -7;
    return;
S350:
    *x = xn;
    *ierr = -8;
    return;
S360:
    *x = xmax;
    *ierr = -8;
    return;
}
//****************************************************************************80

void gamma_inc_values ( int *n_data, double *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function P(A,X) is defined as:
//
//      PN(A,X) = 1/GAMMA(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= PN(A,X) <= 1
//
//    and
//
//      PN(A,INFINITY) = 1.0
//
//    Mathematica can compute this value as
//
//      1 - GammaRegularized[A,X]
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 20

  double a_vec[N_MAX] = {
    0.1E+00,  0.1E+00,  0.1E+00,  0.5E+00,
    0.5E+00,  0.5E+00,  1.0E+00,  1.0E+00,
    1.0E+00,  1.1E+00,  1.1E+00,  1.1E+00,
    2.0E+00,  2.0E+00,  2.0E+00,  6.0E+00,
    6.0E+00, 11.0E+00, 26.0E+00, 41.0E+00 };
  double fx_vec[N_MAX] = {
    0.7420263E+00, 0.9119753E+00, 0.9898955E+00, 0.2931279E+00,
    0.7656418E+00, 0.9921661E+00, 0.0951626E+00, 0.6321206E+00,
    0.9932621E+00, 0.0757471E+00, 0.6076457E+00, 0.9933425E+00,
    0.0091054E+00, 0.4130643E+00, 0.9931450E+00, 0.0387318E+00,
    0.9825937E+00, 0.9404267E+00, 0.4863866E+00, 0.7359709E+00 };
  double x_vec[N_MAX] = {
    3.1622777E-02, 3.1622777E-01, 1.5811388E+00, 7.0710678E-02,
    7.0710678E-01, 3.5355339E+00, 0.1000000E+00, 1.0000000E+00,
    5.0000000E+00, 1.0488088E-01, 1.0488088E+00, 5.2440442E+00,
    1.4142136E-01, 1.4142136E+00, 7.0710678E+00, 2.4494897E+00,
    1.2247449E+01, 1.6583124E+01, 2.5495098E+01, 4.4821870E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

double gamma_ln1 ( double *a )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
//
//  Parameters:
//
//    Input, double *A, defines the argument of the function.
//
//    Output, double GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
//
{
  static double p0 = .577215664901533e+00;
  static double p1 = .844203922187225e+00;
  static double p2 = -.168860593646662e+00;
  static double p3 = -.780427615533591e+00;
  static double p4 = -.402055799310489e+00;
  static double p5 = -.673562214325671e-01;
  static double p6 = -.271935708322958e-02;
  static double q1 = .288743195473681e+01;
  static double q2 = .312755088914843e+01;
  static double q3 = .156875193295039e+01;
  static double q4 = .361951990101499e+00;
  static double q5 = .325038868253937e-01;
  static double q6 = .667465618796164e-03;
  static double r0 = .422784335098467e+00;
  static double r1 = .848044614534529e+00;
  static double r2 = .565221050691933e+00;
  static double r3 = .156513060486551e+00;
  static double r4 = .170502484022650e-01;
  static double r5 = .497958207639485e-03;
  static double s1 = .124313399877507e+01;
  static double s2 = .548042109832463e+00;
  static double s3 = .101552187439830e+00;
  static double s4 = .713309612391000e-02;
  static double s5 = .116165475989616e-03;
  static double gamln1,w,x;

    if(*a >= 0.6e0) goto S10;
    w = ((((((p6**a+p5)**a+p4)**a+p3)**a+p2)**a+p1)**a+p0)/((((((q6**a+q5)**a+
      q4)**a+q3)**a+q2)**a+q1)**a+1.0e0);
    gamln1 = -(*a*w);
    return gamln1;
S10:
    x = *a-0.5e0-0.5e0;
    w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/(((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x
      +1.0e0);
    gamln1 = x*w;
    return gamln1;
}
//****************************************************************************80

double gamma_log ( double *a )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *A, the argument of the function.
//    A should be positive.
//
//    Output, double GAMMA_LOG, the value of ln ( Gamma ( A ) ).
//
{
  static double c0 = .833333333333333e-01;
  static double c1 = -.277777777760991e-02;
  static double c2 = .793650666825390e-03;
  static double c3 = -.595202931351870e-03;
  static double c4 = .837308034031215e-03;
  static double c5 = -.165322962780713e-02;
  static double d = .418938533204673e0;
  static double gamln,t,w;
  static int i,n;
  static double T1;

    if(*a > 0.8e0) goto S10;
    gamln = gamma_ln1 ( a ) - log ( *a );
    return gamln;
S10:
    if(*a > 2.25e0) goto S20;
    t = *a-0.5e0-0.5e0;
    gamln = gamma_ln1 ( &t );
    return gamln;
S20:
    if(*a >= 10.0e0) goto S40;
    n = ( int ) ( *a - 1.25e0 );
    t = *a;
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        t -= 1.0e0;
        w = t*w;
    }
    T1 = t-1.0e0;
    gamln = gamma_ln1 ( &T1 ) + log ( w );
    return gamln;
S40:
    t = pow(1.0e0/ *a,2.0);
    w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ *a;
    gamln = d+w+(*a-0.5e0)*(log(*a)-1.0e0);
    return gamln;
}
//****************************************************************************80

void gamma_rat1 ( double *a, double *x, double *r, double *p, double *q,
  double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
//
//  Parameters:
//
//    Input, double *A, *X, the parameters of the functions.
//    It is assumed that A <= 1.
//
//    Input, double *R, the value exp(-X) * X**A / Gamma(A).
//
//    Output, double *P, *Q, the values of P(A,X) and Q(A,X).
//
//    Input, double *EPS, the tolerance.
//
{
  static int K2 = 0;
  static double a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,t,tol,w,z,T1,T3;

    if(*a**x == 0.0e0) goto S120;
    if(*a == 0.5e0) goto S100;
    if(*x < 1.1e0) goto S10;
    goto S60;
S10:
//
//             TAYLOR SERIES FOR P(A,X)/X**A
//
    an = 3.0e0;
    c = *x;
    sum = *x/(*a+3.0e0);
    tol = 0.1e0**eps/(*a+1.0e0);
S20:
    an += 1.0e0;
    c = -(c*(*x/an));
    t = c/(*a+an);
    sum += t;
    if(fabs(t) > tol) goto S20;
    j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
    z = *a*log(*x);
    h = gam1(a);
    g = 1.0e0+h;
    if(*x < 0.25e0) goto S30;
    if(*a < *x/2.59e0) goto S50;
    goto S40;
S30:
    if(z > -.13394e0) goto S50;
S40:
    w = exp(z);
    *p = w*g*(0.5e0+(0.5e0-j));
    *q = 0.5e0+(0.5e0-*p);
    return;
S50:
    l = rexp(&z);
    w = 0.5e0+(0.5e0+l);
    *q = (w*j-l)*g-h;
    if(*q < 0.0e0) goto S90;
    *p = 0.5e0+(0.5e0-*q);
    return;
S60:
//
//              CONTINUED FRACTION EXPANSION
//
    a2nm1 = a2n = 1.0e0;
    b2nm1 = *x;
    b2n = *x+(1.0e0-*a);
    c = 1.0e0;
S70:
    a2nm1 = *x*a2n+c*a2nm1;
    b2nm1 = *x*b2n+c*b2nm1;
    am0 = a2nm1/b2nm1;
    c += 1.0e0;
    cma = c-*a;
    a2n = a2nm1+cma*a2n;
    b2n = b2nm1+cma*b2n;
    an0 = a2n/b2n;
    if(fabs(an0-am0) >= *eps*an0) goto S70;
    *q = *r*an0;
    *p = 0.5e0+(0.5e0-*q);
    return;
S80:
//
//                SPECIAL CASES
//
    *p = 0.0e0;
    *q = 1.0e0;
    return;
S90:
    *p = 1.0e0;
    *q = 0.0e0;
    return;
S100:
    if(*x >= 0.25e0) goto S110;
    T1 = sqrt(*x);
    *p = error_f ( &T1 );
    *q = 0.5e0+(0.5e0-*p);
    return;
S110:
    T3 = sqrt(*x);
    *q = error_fc ( &K2, &T3 );
    *p = 0.5e0+(0.5e0-*q);
    return;
S120:
    if(*x <= *a) goto S80;
    goto S90;
}
//****************************************************************************80

void gamma_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_VALUES returns some values of the Gamma function.
//
//  Definition:
//
//    GAMMA(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) EXP(-T) dT
//
//  Recursion:
//
//    GAMMA(X+1) = X*GAMMA(X)
//
//  Restrictions:
//
//    0 < X ( a software restriction).
//
//  Special values:
//
//    GAMMA(0.5) = sqrt(PI)
//
//    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 18

  double fx_vec[N_MAX] = {
    4.590845E+00,     2.218160E+00,     1.489192E+00,     1.164230E+00,
    1.0000000000E+00, 0.9513507699E+00, 0.9181687424E+00, 0.8974706963E+00,
    0.8872638175E+00, 0.8862269255E+00, 0.8935153493E+00, 0.9086387329E+00,
    0.9313837710E+00, 0.9617658319E+00, 1.0000000000E+00, 3.6288000E+05,
    1.2164510E+17,    8.8417620E+30 };
  double x_vec[N_MAX] = {
    0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00,
    1.0E+00,  1.1E+00,  1.2E+00,  1.3E+00,
    1.4E+00,  1.5E+00,  1.6E+00,  1.7E+00,
    1.8E+00,  1.9E+00,  2.0E+00, 10.0E+00,
   20.0E+00, 30.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

double gamma_x ( double *a )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_X evaluates the gamma function.
//
//  Discussion:
//
//    This routine was renamed from "GAMMA" to avoid a conflict with the
//    C/C++ math library routine.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double *A, the argument of the Gamma function.
//
//    Output, double GAMMA_X, the value of the Gamma function.
//
{
  static double d = .41893853320467274178e0;
  static double pi = 3.1415926535898e0;
  static double r1 = .820756370353826e-03;
  static double r2 = -.595156336428591e-03;
  static double r3 = .793650663183693e-03;
  static double r4 = -.277777777770481e-02;
  static double r5 = .833333333333333e-01;
  static double p[7] = {
    .539637273585445e-03,.261939260042690e-02,.204493667594920e-01,
    .730981088720487e-01,.279648642639792e+00,.553413866010467e+00,1.0e0
  };
  static double q[7] = {
    -.832979206704073e-03,.470059485860584e-02,.225211131035340e-01,
    -.170458969313360e+00,-.567902761974940e-01,.113062953091122e+01,1.0e0
  };
  static int K2 = 3;
  static int K3 = 0;
  static double Xgamm,bot,g,lnx,s,t,top,w,x,z;
  static int i,j,m,n,T1;

    Xgamm = 0.0e0;
    x = *a;
    if(fabs(*a) >= 15.0e0) goto S110;
//
//            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
//
    t = 1.0e0;
    m = fifidint(*a)-1;
//
//     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
//
    T1 = m;
    if(T1 < 0) goto S40;
    else if(T1 == 0) goto S30;
    else  goto S10;
S10:
    for ( j = 1; j <= m; j++ )
    {
        x -= 1.0e0;
        t = x*t;
    }
S30:
    x -= 1.0e0;
    goto S80;
S40:
//
//     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
//
    t = *a;
    if(*a > 0.0e0) goto S70;
    m = -m-1;
    if(m == 0) goto S60;
    for ( j = 1; j <= m; j++ )
    {
        x += 1.0e0;
        t = x*t;
    }
S60:
    x += (0.5e0+0.5e0);
    t = x*t;
    if(t == 0.0e0) return Xgamm;
S70:
//
//     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
//     CODE MAY BE OMITTED IF DESIRED.
//
    if(fabs(t) >= 1.e-30) goto S80;
    if(fabs(t)*dpmpar(&K2) <= 1.0001e0) return Xgamm;
    Xgamm = 1.0e0/t;
    return Xgamm;
S80:
//
//     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
//
    top = p[0];
    bot = q[0];
    for ( i = 1; i < 7; i++ )
    {
        top = p[i]+x*top;
        bot = q[i]+x*bot;
    }
    Xgamm = top/bot;
//
//     TERMINATION
//
    if(*a < 1.0e0) goto S100;
    Xgamm *= t;
    return Xgamm;
S100:
    Xgamm /= t;
    return Xgamm;
S110:
//
//  EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
//
    if(fabs(*a) >= 1.e3) return Xgamm;
    if(*a > 0.0e0) goto S120;
    x = -*a;
    n = ( int ) x;
    t = x-(double)n;
    if(t > 0.9e0) t = 1.0e0-t;
    s = sin(pi*t)/pi;
    if(fifmod(n,2) == 0) s = -s;
    if(s == 0.0e0) return Xgamm;
S120:
//
//     COMPUTE THE MODIFIED ASYMPTOTIC SUM
//
    t = 1.0e0/(x*x);
    g = ((((r1*t+r2)*t+r3)*t+r4)*t+r5)/x;
//
//     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
//     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
//
    lnx = log(x);
//
//  FINAL ASSEMBLY
//
    z = x;
    g = d+g+(z-0.5e0)*(lnx-1.e0);
    w = g;
    t = g-w;
    if(w > 0.99999e0*exparg(&K3)) return Xgamm;
    Xgamm = exp(w)*(1.0e0+t);
    if(*a < 0.0e0) Xgamm = 1.0e0/(Xgamm*s)/x;
    return Xgamm;
}
//****************************************************************************80

double gsumln ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    GSUMLN evaluates the function ln(Gamma(A + B)).
//
//  Discussion:
//
//    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
//
//  Parameters:
//
//    Input, double *A, *B, values whose sum is the argument of
//    the Gamma function.
//
//    Output, double GSUMLN, the value of ln(Gamma(A+B)).
//
{
  static double gsumln,x,T1,T2;

    x = *a+*b-2.e0;
    if(x > 0.25e0) goto S10;
    T1 = 1.0e0+x;
    gsumln = gamma_ln1 ( &T1 );
    return gsumln;
S10:
    if(x > 1.25e0) goto S20;
    gsumln = gamma_ln1 ( &x ) + alnrel ( &x );
    return gsumln;
S20:
    T2 = x-1.0e0;
    gsumln = gamma_ln1 ( &T2 ) + log ( x * ( 1.0e0 + x ) );
    return gsumln;
}
//****************************************************************************80

int ipmpar ( int *i )

//****************************************************************************80
//
//  Purpose:
//
//    IPMPAR returns integer machine constants.
//
//  Discussion:
//
//    Input arguments 1 through 3 are queries about integer arithmetic.
//    We assume integers are represented in the N-digit, base-A form
//
//      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
//
//    where 0 <= X(0:N-1) < A.
//
//    Then:
//
//      IPMPAR(1) = A, the base of integer arithmetic;
//      IPMPAR(2) = N, the number of base A digits;
//      IPMPAR(3) = A**N - 1, the largest magnitude.
//
//    It is assumed that the single and double precision floating
//    point arithmetics have the same base, say B, and that the
//    nonzero numbers are represented in the form
//
//      sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
//
//    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
//    EMIN <= E <= EMAX.
//
//    Input argument 4 is a query about the base of real arithmetic:
//
//      IPMPAR(4) = B, the base of single and double precision arithmetic.
//
//    Input arguments 5 through 7 are queries about single precision
//    floating point arithmetic:
//
//     IPMPAR(5) = M, the number of base B digits for single precision.
//     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
//     IPMPAR(7) = EMAX, the largest exponent E for single precision.
//
//    Input arguments 8 through 10 are queries about double precision
//    floating point arithmetic:
//
//     IPMPAR(8) = M, the number of base B digits for double precision.
//     IPMPAR(9) = EMIN, the smallest exponent E for double precision.
//     IPMPAR(10) = EMAX, the largest exponent E for double precision.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, and Norman Schryer,
//    Algorithm 528,
//    Framework for a Portable FORTRAN Subroutine Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, 1978, pages 176-188.
//
//  Parameters:
//
//    Input, int *I, the index of the desired constant.
//
//    Output, int IPMPAR, the value of the desired constant.
//
{
  static int imach[11];
  static int ipmpar;
//     MACHINE CONSTANTS FOR AMDAHL MACHINES.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 16;
//   imach[5] = 6;
//   imach[6] = -64;
//   imach[7] = 63;
//   imach[8] = 14;
//   imach[9] = -64;
//   imach[10] = 63;
//
//     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
//       PC 7300, AND AT&T 6300.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -125;
//   imach[7] = 128;
//   imach[8] = 53;
//   imach[9] = -1021;
//   imach[10] = 1024;
//
//     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
//
//   imach[1] = 2;
//   imach[2] = 33;
//   imach[3] = 8589934591;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -256;
//   imach[7] = 255;
//   imach[8] = 60;
//   imach[9] = -256;
//   imach[10] = 255;
//
//     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
//
//   imach[1] = 2;
//   imach[2] = 39;
//   imach[3] = 549755813887;
//   imach[4] = 8;
//   imach[5] = 13;
//   imach[6] = -50;
//   imach[7] = 76;
//   imach[8] = 26;
//   imach[9] = -50;
//   imach[10] = 76;
//
//     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
//
//   imach[1] = 2;
//   imach[2] = 39;
//   imach[3] = 549755813887;
//   imach[4] = 8;
//   imach[5] = 13;
//   imach[6] = -50;
//   imach[7] = 76;
//   imach[8] = 26;
//   imach[9] = -32754;
//   imach[10] = 32780;
//
//     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
//       60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
//       ARITHMETIC (NOS OPERATING SYSTEM).
//
//   imach[1] = 2;
//   imach[2] = 48;
//   imach[3] = 281474976710655;
//   imach[4] = 2;
//   imach[5] = 48;
//   imach[6] = -974;
//   imach[7] = 1070;
//   imach[8] = 95;
//   imach[9] = -926;
//   imach[10] = 1070;
//
//     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
//       ARITHMETIC (NOS/VE OPERATING SYSTEM).
//
//   imach[1] = 2;
//   imach[2] = 63;
//   imach[3] = 9223372036854775807;
//   imach[4] = 2;
//   imach[5] = 48;
//   imach[6] = -4096;
//   imach[7] = 4095;
//   imach[8] = 96;
//   imach[9] = -4096;
//   imach[10] = 4095;
//
//     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
//
//   imach[1] = 2;
//   imach[2] = 63;
//   imach[3] = 9223372036854775807;
//   imach[4] = 2;
//   imach[5] = 47;
//   imach[6] = -8189;
//   imach[7] = 8190;
//   imach[8] = 94;
//   imach[9] = -8099;
//   imach[10] = 8190;
//
//     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
//
//   imach[1] = 2;
//   imach[2] = 15;
//   imach[3] = 32767;
//   imach[4] = 16;
//   imach[5] = 6;
//   imach[6] = -64;
//   imach[7] = 63;
//   imach[8] = 14;
//   imach[9] = -64;
//   imach[10] = 63;
//
//     MACHINE CONSTANTS FOR THE HARRIS 220.
//
//   imach[1] = 2;
//   imach[2] = 23;
//   imach[3] = 8388607;
//   imach[4] = 2;
//   imach[5] = 23;
//   imach[6] = -127;
//   imach[7] = 127;
//   imach[8] = 38;
//   imach[9] = -127;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
//       AND DPS 8/70 SERIES.
//
//   imach[1] = 2;
//   imach[2] = 35;
//   imach[3] = 34359738367;
//   imach[4] = 2;
//   imach[5] = 27;
//   imach[6] = -127;
//   imach[7] = 127;
//   imach[8] = 63;
//   imach[9] = -127;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE HP 2100
//       3 WORD DOUBLE PRECISION OPTION WITH FTN4
//
//   imach[1] = 2;
//   imach[2] = 15;
//   imach[3] = 32767;
//   imach[4] = 2;
//   imach[5] = 23;
//   imach[6] = -128;
//   imach[7] = 127;
//   imach[8] = 39;
//   imach[9] = -128;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE HP 2100
//       4 WORD DOUBLE PRECISION OPTION WITH FTN4
//
//   imach[1] = 2;
//   imach[2] = 15;
//   imach[3] = 32767;
//   imach[4] = 2;
//   imach[5] = 23;
//   imach[6] = -128;
//   imach[7] = 127;
//   imach[8] = 55;
//   imach[9] = -128;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE HP 9000.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -126;
//   imach[7] = 128;
//   imach[8] = 53;
//   imach[9] = -1021;
//   imach[10] = 1024;
//
//     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
//       THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
//       5/7/9 AND THE SEL SYSTEMS 85/86.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 16;
//   imach[5] = 6;
//   imach[6] = -64;
//   imach[7] = 63;
//   imach[8] = 14;
//   imach[9] = -64;
//   imach[10] = 63;
//
//     MACHINE CONSTANTS FOR THE IBM PC.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -125;
//   imach[7] = 128;
//   imach[8] = 53;
//   imach[9] = -1021;
//   imach[10] = 1024;
//
//     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
//       MACFORTRAN II.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -125;
//   imach[7] = 128;
//   imach[8] = 53;
//   imach[9] = -1021;
//   imach[10] = 1024;
//
//     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -127;
//   imach[7] = 127;
//   imach[8] = 56;
//   imach[9] = -127;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
//
//   imach[1] = 2;
//   imach[2] = 35;
//   imach[3] = 34359738367;
//   imach[4] = 2;
//   imach[5] = 27;
//   imach[6] = -128;
//   imach[7] = 127;
//   imach[8] = 54;
//   imach[9] = -101;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
//
//   imach[1] = 2;
//   imach[2] = 35;
//   imach[3] = 34359738367;
//   imach[4] = 2;
//   imach[5] = 27;
//   imach[6] = -128;
//   imach[7] = 127;
//   imach[8] = 62;
//   imach[9] = -128;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
//       32-BIT INTEGER ARITHMETIC.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -127;
//   imach[7] = 127;
//   imach[8] = 56;
//   imach[9] = -127;
//   imach[10] = 127;
//
//     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -125;
//   imach[7] = 128;
//   imach[8] = 53;
//   imach[9] = -1021;
//   imach[10] = 1024;
//
//     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
//       SERIES (MIPS R3000 PROCESSOR).
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -125;
//   imach[7] = 128;
//   imach[8] = 53;
//   imach[9] = -1021;
//   imach[10] = 1024;
//
//     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
//       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
//       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;

//     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
//
//   imach[1] = 2;
//   imach[2] = 35;
//   imach[3] = 34359738367;
//   imach[4] = 2;
//   imach[5] = 27;
//   imach[6] = -128;
//   imach[7] = 127;
//   imach[8] = 60;
//   imach[9] = -1024;
//   imach[10] = 1023;
//
//     MACHINE CONSTANTS FOR THE VAX 11/780.
//
//   imach[1] = 2;
//   imach[2] = 31;
//   imach[3] = 2147483647;
//   imach[4] = 2;
//   imach[5] = 24;
//   imach[6] = -127;
//   imach[7] = 127;
//   imach[8] = 56;
//   imach[9] = -127;
//   imach[10] = 127;
//
    ipmpar = imach[*i];
    return ipmpar;
}
//****************************************************************************80

void negative_binomial_cdf_values ( int *n_data, int *f, int *s, double *p,
  double *cdf )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
//
//  Discussion:
//
//    Assume that a coin has a probability P of coming up heads on
//    any one trial.  Suppose that we plan to flip the coin until we
//    achieve a total of S heads.  If we let F represent the number of
//    tails that occur in this process, then the value of F satisfies
//    a negative binomial PDF:
//
//      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
//
//    The negative binomial CDF is the probability that there are F or
//    fewer failures upon the attainment of the S-th success.  Thus,
//
//      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
//
//  Modified:
//
//    07 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    F C Powell,
//    Statistical Tables for Sociology, Biology and Physical Sciences,
//    Cambridge University Press, 1982.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *F, the maximum number of failures.
//
//    Output, int *S, the number of successes.
//
//    Output, double *P, the probability of a success on one trial.
//
//    Output, double *CDF, the probability of at most F failures before the
//    S-th success.
//
{
# define N_MAX 27

  double cdf_vec[N_MAX] = {
    0.6367, 0.3633, 0.1445,
    0.5000, 0.2266, 0.0625,
    0.3438, 0.1094, 0.0156,
    0.1792, 0.0410, 0.0041,
    0.0705, 0.0109, 0.0007,
    0.9862, 0.9150, 0.7472,
    0.8499, 0.5497, 0.2662,
    0.6513, 0.2639, 0.0702,
    1.0000, 0.0199, 0.0001 };
  int f_vec[N_MAX] = {
     4,  3,  2,
     3,  2,  1,
     2,  1,  0,
     2,  1,  0,
     2,  1,  0,
    11, 10,  9,
    17, 16, 15,
     9,  8,  7,
     2,  1,  0 };
  double p_vec[N_MAX] = {
    0.50, 0.50, 0.50,
    0.50, 0.50, 0.50,
    0.50, 0.50, 0.50,
    0.40, 0.40, 0.40,
    0.30, 0.30, 0.30,
    0.30, 0.30, 0.30,
    0.10, 0.10, 0.10,
    0.10, 0.10, 0.10,
    0.01, 0.01, 0.01 };
  int s_vec[N_MAX] = {
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    0, 1, 2 };

  if (*n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *f = 0;
    *s = 0;
    *p = 0.0E+00;
    *cdf = 0.0E+00;
  }
  else
  {
    *f = f_vec[*n_data-1];
    *s = s_vec[*n_data-1];
    *p = p_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void normal_cdf_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF_VALUES returns some values of the Normal CDF.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output double *FX, the value of the function.
//
{
# define N_MAX 13

  double fx_vec[N_MAX] = {
    0.500000000000000E+00, 0.539827837277029E+00, 0.579259709439103E+00,
    0.617911422188953E+00, 0.655421741610324E+00, 0.691462461274013E+00,
    0.725746882249927E+00, 0.758036347776927E+00, 0.788144601416604E+00,
    0.815939874653241E+00, 0.841344746068543E+00, 0.933192798731142E+00,
    0.977249868051821E+00 };
  double x_vec[N_MAX] = {
    0.00E+00, 0.10E+00, 0.20E+00,
    0.30E+00, 0.40E+00, 0.50E+00,
    0.60E+00, 0.70E+00, 0.80E+00,
    0.90E+00, 1.00E+00, 1.50E+00,
    2.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void poisson_cdf_values ( int *n_data, double *a, int *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF_VALUES returns some values of the Poisson CDF.
//
//  Discussion:
//
//    CDF(X)(A) is the probability of at most X successes in unit time,
//    given that the expected mean number of successes is A.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 653-658.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, the parameter of the function.
//
//    Output, int *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double a_vec[N_MAX] = {
    0.02E+00, 0.10E+00, 0.10E+00, 0.50E+00,
    0.50E+00, 0.50E+00, 1.00E+00, 1.00E+00,
    1.00E+00, 1.00E+00, 2.00E+00, 2.00E+00,
    2.00E+00, 2.00E+00, 5.00E+00, 5.00E+00,
    5.00E+00, 5.00E+00, 5.00E+00, 5.00E+00,
    5.00E+00 };
  double fx_vec[N_MAX] = {
    0.980E+00, 0.905E+00, 0.995E+00, 0.607E+00,
    0.910E+00, 0.986E+00, 0.368E+00, 0.736E+00,
    0.920E+00, 0.981E+00, 0.135E+00, 0.406E+00,
    0.677E+00, 0.857E+00, 0.007E+00, 0.040E+00,
    0.125E+00, 0.265E+00, 0.441E+00, 0.616E+00,
    0.762E+00 };
  int x_vec[N_MAX] = {
     0, 0, 1, 0,
     1, 2, 0, 1,
     2, 3, 0, 1,
     2, 3, 0, 1,
     2, 3, 4, 5,
     6 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *x = 0;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

double psi ( double *xx )

//****************************************************************************80
//
//  Purpose:
//
//    PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
//
//  Discussion:
//
//    The main computation involves evaluation of rational Chebyshev
//    approximations.  PSI was written at Argonne National Laboratory
//    for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
//
//  Reference:
//
//    William Cody, Strecok and Thacher,
//    Chebyshev Approximations for the Psi Function,
//    Mathematics of Computation,
//    Volume 27, 1973, pages 123-127.
//
//  Parameters:
//
//    Input, double *XX, the argument of the psi function.
//
//    Output, double PSI, the value of the psi function.  PSI
//    is assigned the value 0 when the psi function is undefined.
//
{
  static double dx0 = 1.461632144968362341262659542325721325e0;
  static double piov4 = .785398163397448e0;
  static double p1[7] = {
    .895385022981970e-02,.477762828042627e+01,.142441585084029e+03,
    .118645200713425e+04,.363351846806499e+04,.413810161269013e+04,
    .130560269827897e+04
  };
  static double p2[4] = {
    -.212940445131011e+01,-.701677227766759e+01,-.448616543918019e+01,
    -.648157123766197e+00
  };
  static double q1[6] = {
    .448452573429826e+02,.520752771467162e+03,.221000799247830e+04,
    .364127349079381e+04,.190831076596300e+04,.691091682714533e-05
  };
  static double q2[4] = {
    .322703493791143e+02,.892920700481861e+02,.546117738103215e+02,
    .777788548522962e+01
  };
  static int K1 = 3;
  static int K2 = 1;
  static double psi,aug,den,sgn,upper,w,x,xmax1,xmx0,xsmall,z;
  static int i,m,n,nq;
//
//     MACHINE DEPENDENT CONSTANTS ...
//        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
//                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
//                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
//                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
//                 PSI MAY BE REPRESENTED AS ALOG(X).
//        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
//                 MAY BE REPRESENTED BY 1/X.
//
    xmax1 = ipmpar(&K1);
    xmax1 = fifdmin1(xmax1,1.0e0/dpmpar(&K2));
    xsmall = 1.e-9;
    x = *xx;
    aug = 0.0e0;
    if(x >= 0.5e0) goto S50;
//
//     X .LT. 0.5,  USE REFLECTION FORMULA
//     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
//
    if(fabs(x) > xsmall) goto S10;
    if(x == 0.0e0) goto S100;
//
//     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
//     FOR  PI*COTAN(PI*X)
//
    aug = -(1.0e0/x);
    goto S40;
S10:
//
//     REDUCTION OF ARGUMENT FOR COTAN
//
    w = -x;
    sgn = piov4;
    if(w > 0.0e0) goto S20;
    w = -w;
    sgn = -sgn;
S20:
//
//     MAKE AN ERROR EXIT IF X .LE. -XMAX1
//
    if(w >= xmax1) goto S100;
    nq = fifidint(w);
    w -= (double)nq;
    nq = fifidint(w*4.0e0);
    w = 4.0e0*(w-(double)nq*.25e0);
//
//     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
//     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
//     QUADRANT AND DETERMINE SIGN
//
    n = nq/2;
    if(n+n != nq) w = 1.0e0-w;
    z = piov4*w;
    m = n/2;
    if(m+m != n) sgn = -sgn;
//
//     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
//
    n = (nq+1)/2;
    m = n/2;
    m += m;
    if(m != n) goto S30;
//
//     CHECK FOR SINGULARITY
//
    if(z == 0.0e0) goto S100;
//
//     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
//     SIN/COS AS A SUBSTITUTE FOR TAN
//
    aug = sgn*(cos(z)/sin(z)*4.0e0);
    goto S40;
S30:
    aug = sgn*(sin(z)/cos(z)*4.0e0);
S40:
    x = 1.0e0-x;
S50:
    if(x > 3.0e0) goto S70;
//
//     0.5 .LE. X .LE. 3.0
//
    den = x;
    upper = p1[0]*x;
    for ( i = 1; i <= 5; i++ )
    {
        den = (den+q1[i-1])*x;
        upper = (upper+p1[i+1-1])*x;
    }
    den = (upper+p1[6])/(den+q1[5]);
    xmx0 = x-dx0;
    psi = den*xmx0+aug;
    return psi;
S70:
//
//     IF X .GE. XMAX1, PSI = LN(X)
//
    if(x >= xmax1) goto S90;
//
//     3.0 .LT. X .LT. XMAX1
//
    w = 1.0e0/(x*x);
    den = w;
    upper = p2[0]*w;
    for ( i = 1; i <= 3; i++ )
    {
        den = (den+q2[i-1])*w;
        upper = (upper+p2[i+1-1])*w;
    }
    aug = upper/(den+q2[3])-0.5e0/x+aug;
S90:
    psi = aug+log(x);
    return psi;
S100:
//
//     ERROR RETURN
//
    psi = 0.0e0;
    return psi;
}
//****************************************************************************80

void psi_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    PSI_VALUES returns some values of the Psi or Digamma function.
//
//  Discussion:
//
//    PSI(X) = d LN ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
//
//    PSI(1) = - Euler's constant.
//
//    PSI(X+1) = PSI(X) + 1 / X.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 11

  double fx_vec[N_MAX] = {
    -0.5772156649E+00, -0.4237549404E+00, -0.2890398966E+00,
    -0.1691908889E+00, -0.0613845446E+00, -0.0364899740E+00,
     0.1260474528E+00,  0.2085478749E+00,  0.2849914333E+00,
     0.3561841612E+00,  0.4227843351E+00 };
  double x_vec[N_MAX] = {
    1.0E+00,  1.1E+00,  1.2E+00,
    1.3E+00,  1.4E+00,  1.5E+00,
    1.6E+00,  1.7E+00,  1.8E+00,
    1.9E+00,  2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

double rcomp ( double *a, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    RCOMP evaluates exp(-X) * X**A / Gamma(A).
//
//  Parameters:
//
//    Input, double *A, *X, arguments of the quantity to be computed.
//
//    Output, double RCOMP, the value of exp(-X) * X**A / Gamma(A).
//
//  Local parameters:
//
//    RT2PIN = 1/SQRT(2*PI)
//
{
  static double rt2pin = .398942280401433e0;
  static double rcomp,t,t1,u;
    rcomp = 0.0e0;
    if(*a >= 20.0e0) goto S20;
    t = *a*log(*x)-*x;
    if(*a >= 1.0e0) goto S10;
    rcomp = *a*exp(t)*(1.0e0+gam1(a));
    return rcomp;
S10:
    rcomp = exp(t)/ gamma_x(a);
    return rcomp;
S20:
    u = *x/ *a;
    if(u == 0.0e0) return rcomp;
    t = pow(1.0e0/ *a,2.0);
    t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
    t1 -= (*a*rlog(&u));
    rcomp = rt2pin*sqrt(*a)*exp(t1);
    return rcomp;
}
//****************************************************************************80

double rexp ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    REXP evaluates the function EXP(X) - 1.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, double *X, the argument of the function.
//
//    Output, double REXP, the value of EXP(X)-1.
//
{
  static double p1 = .914041914819518e-09;
  static double p2 = .238082361044469e-01;
  static double q1 = -.499999999085958e+00;
  static double q2 = .107141568980644e+00;
  static double q3 = -.119041179760821e-01;
  static double q4 = .595130811860248e-03;
  static double rexp,w;

    if(fabs(*x) > 0.15e0) goto S10;
    rexp = *x*(((p2**x+p1)**x+1.0e0)/((((q4**x+q3)**x+q2)**x+q1)**x+1.0e0));
    return rexp;
S10:
    w = exp(*x);
    if(*x > 0.0e0) goto S20;
    rexp = w-0.5e0-0.5e0;
    return rexp;
S20:
    rexp = w*(0.5e0+(0.5e0-1.0e0/w));
    return rexp;
}
//****************************************************************************80

double rlog ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    RLOG computes  X - 1 - LN(X).
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, double *X, the argument of the function.
//
//    Output, double RLOG, the value of the function.
//
{
  static double a = .566749439387324e-01;
  static double b = .456512608815524e-01;
  static double p0 = .333333333333333e+00;
  static double p1 = -.224696413112536e+00;
  static double p2 = .620886815375787e-02;
  static double q1 = -.127408923933623e+01;
  static double q2 = .354508718369557e+00;
  static double rlog,r,t,u,w,w1;

    if(*x < 0.61e0 || *x > 1.57e0) goto S40;
    if(*x < 0.82e0) goto S10;
    if(*x > 1.18e0) goto S20;
//
//              ARGUMENT REDUCTION
//
    u = *x-0.5e0-0.5e0;
    w1 = 0.0e0;
    goto S30;
S10:
    u = *x-0.7e0;
    u /= 0.7e0;
    w1 = a-u*0.3e0;
    goto S30;
S20:
    u = 0.75e0**x-1.e0;
    w1 = b+u/3.0e0;
S30:
//
//               SERIES EXPANSION
//
    r = u/(u+2.0e0);
    t = r*r;
    w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
    rlog = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
    return rlog;
S40:
    r = *x-0.5e0-0.5e0;
    rlog = r-log(*x);
    return rlog;
}
//****************************************************************************80

double rlog1 ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    RLOG1 evaluates the function X - ln ( 1 + X ).
//
//  Parameters:
//
//    Input, double *X, the argument.
//
//    Output, double RLOG1, the value of X - ln ( 1 + X ).
//
{
  static double a = .566749439387324e-01;
  static double b = .456512608815524e-01;
  static double p0 = .333333333333333e+00;
  static double p1 = -.224696413112536e+00;
  static double p2 = .620886815375787e-02;
  static double q1 = -.127408923933623e+01;
  static double q2 = .354508718369557e+00;
  static double rlog1,h,r,t,w,w1;

    if(*x < -0.39e0 || *x > 0.57e0) goto S40;
    if(*x < -0.18e0) goto S10;
    if(*x > 0.18e0) goto S20;
//
//              ARGUMENT REDUCTION
//
    h = *x;
    w1 = 0.0e0;
    goto S30;
S10:
    h = *x+0.3e0;
    h /= 0.7e0;
    w1 = a-h*0.3e0;
    goto S30;
S20:
    h = 0.75e0**x-0.25e0;
    w1 = b+h/3.0e0;
S30:
//
//               SERIES EXPANSION
//
    r = h/(h+2.0e0);
    t = r*r;
    w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
    rlog1 = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
    return rlog1;
S40:
    w = *x+0.5e0+0.5e0;
    rlog1 = *x-log(w);
    return rlog1;
}
//****************************************************************************80

void student_cdf_values ( int *n_data, int *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CDF_VALUES returns some values of the Student CDF.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 13

  int a_vec[N_MAX] = {
    1, 2, 3, 4,
    5, 2, 5, 2,
    5, 2, 3, 4,
    5 };
  double fx_vec[N_MAX] = {
    0.60E+00, 0.60E+00, 0.60E+00, 0.60E+00,
    0.60E+00, 0.75E+00, 0.75E+00, 0.95E+00,
    0.95E+00, 0.99E+00, 0.99E+00, 0.99E+00,
    0.99E+00 };
  double x_vec[N_MAX] = {
    0.325E+00, 0.289E+00, 0.277E+00, 0.271E+00,
    0.267E+00, 0.816E+00, 0.727E+00, 2.920E+00,
    2.015E+00, 6.965E+00, 4.541E+00, 3.747E+00,
    3.365E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double stvaln ( double *p )

//****************************************************************************80
//
//  Purpose:
//
//    STVALN provides starting values for the inverse of the normal distribution.
//
//  Discussion:
//
//    The routine returns X such that
//      P = CUMNOR(X),
//    that is,
//      P = Integral from -infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU.
//
//  Reference:
//
//    Kennedy and Gentle,
//    Statistical Computing,
//    Marcel Dekker, NY, 1980, page 95,
//    QA276.4  K46
//
//  Parameters:
//
//    Input, double *P, the probability whose normal deviate
//    is sought.
//
//    Output, double STVALN, the normal deviate whose probability
//    is P.
//
{
  static double xden[5] = {
    0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
    0.38560700634e-2
  };
  static double xnum[5] = {
    -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
    -0.453642210148e-4
  };
  static int K1 = 5;
  static double stvaln,sign,y,z;

    if(!(*p <= 0.5e0)) goto S10;
    sign = -1.0e0;
    z = *p;
    goto S20;
S10:
    sign = 1.0e0;
    z = 1.0e0-*p;
S20:
    y = sqrt(-(2.0e0*log(z)));
    stvaln = y+ eval_pol ( xnum, &K1, &y ) / eval_pol ( xden, &K1, &y );
    stvaln = sign*stvaln;
    return stvaln;
}

//****************************************************************************80

double alnorm ( double x, bool upper )

//****************************************************************************80
//
//  Purpose:
//
//    ALNORM computes the cumulative density of the standard normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by David Hill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Hill,
//    Algorithm AS 66:
//    The Normal Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 424-427.
//
//  Parameters:
//
//    Input, double X, is one endpoint of the semi-infinite interval
//    over which the integration takes place.
//
//    Input, bool UPPER, determines whether the upper or lower
//    interval is to be integrated:
//    .TRUE.  => integrate from X to + Infinity;
//    .FALSE. => integrate from - Infinity to X.
//
//    Output, double ALNORM, the integral of the standard normal
//    distribution over the desired interval.
//
{
  double a1 = 5.75885480458;
  double a2 = 2.62433121679;
  double a3 = 5.92885724438;
  double b1 = -29.8213557807;
  double b2 = 48.6959930692;
  double c1 = -0.000000038052;
  double c2 = 0.000398064794;
  double c3 = -0.151679116635;
  double c4 = 4.8385912808;
  double c5 = 0.742380924027;
  double c6 = 3.99019417011;
  double con = 1.28;
  double d1 = 1.00000615302;
  double d2 = 1.98615381364;
  double d3 = 5.29330324926;
  double d4 = -15.1508972451;
  double d5 = 30.789933034;
  double ltone = 7.0;
  double p = 0.398942280444;
  double q = 0.39990348504;
  double r = 0.398942280385;
  bool up;
  double utzero = 18.66;
  double value;
  double y;
  double z;

  up = upper;
  z = x;

  if ( z < 0.0 )
  {
    up = !up;
    z = - z;
  }

  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  y = 0.5 * z * z;

  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y 
      / ( y + a1 + b1 
      / ( y + a2 + b2 
      / ( y + a3 ))));
  }
  else
  {
    value = r * exp ( - y ) 
      / ( z + c1 + d1 
      / ( z + c2 + d2 
      / ( z + c3 + d3 
      / ( z + c4 + d4 
      / ( z + c5 + d5 
      / ( z + c6 ))))));
  }

  if ( !up )
  {
    value = 1.0 - value;
  }

  return value;
}
//****************************************************************************80

double prncst ( double &st, int &idf, double &d,int &ifault)

//****************************************************************************80
//
//  Purpose:
//
//    PRNCST computes the lower tail of noncentral T distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by BE Cooper.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    BE Cooper,
//    Algorithm AS 5:
//    The Integral of the Non-Central T-Distribution,
//    Applied Statistics,
//    Volume 17, Number 2, 1968, page 193.
//
//  Parameters:
//
//    Input, double ST, the argument.
//
//    Input, int IDF, the number of degrees of freedom.
//
//    Input, double D, the noncentrality parameter.
//
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double PRNCST, the value of the lower tail of
//    the noncentral T distribution.
//
//  Local Parameters:
//
//    Local, double G1, 1.0 / sqrt(2.0 * pi)
//
//    Local, double G2, 1.0 / (2.0 * pi)
//
//    Local, double G3, sqrt(2.0 * pi)
//
{
  double a;
  double ak;
  double b;
  double da;
  double drb;
  double emin = 12.5;
  double f;
  double fk;
  double fkm1;
  double fmkm1;
  double fmkm2;
  double g1 = 0.3989422804;
  double g2 = 0.1591549431;
  double g3 = 2.5066282746;
  int ioe;
  int k;
  double rb;
  double sum;
  double value;

  f = ( double ) ( idf );
//
//  For very large IDF, use the normal approximation.
//
  if ( 100 < idf )
  {
    ifault = 1;

    a = sqrt ( 0.5 * f ) * exp ( lgamma ( 0.5 * ( f - 1.0 ) ) 
      - lgamma ( 0.5 * f ) ) * d;

    value = alnorm ( ( st - a ) / sqrt ( f * ( 1.0 + d * d ) 
      / ( f - 2.0 ) - a * a ), false );
    return value;
  }

  ifault = 0;
  ioe = ( idf % 2 );
  a = st / sqrt ( f );
  b = f / ( f + st * st );
  rb = sqrt ( b );
  da = d * a;
  drb = d * rb;

  if ( idf == 1 )
  {
    value = alnorm ( drb, true ) + 2.0 * tfn ( drb, a );
    return value;
  }

  sum = 0.0;

  if ( fabs ( drb ) < emin )
  {
    fmkm2 = a * rb * exp ( - 0.5 * drb * drb ) 
    * alnorm ( a * drb, false ) * g1;
  }
  else
  {
    fmkm2 = 0.0;
  }

  fmkm1 = b * da * fmkm2;

  if ( fabs ( d ) < emin )
  {
    fmkm1 = fmkm1 + b * a * g2 * exp ( - 0.5 * d * d );
  }

  if ( ioe == 0 )
  {
    sum = fmkm2;
  }
  else
  {
    sum = fmkm1;
  }

  ak = 1.0;
  fk = 2.0;

  for ( k = 2; k <= idf - 2; k = k + 2 )
  {
    fkm1 = fk - 1.0;
    fmkm2 = b * ( da * ak * fmkm1 + fmkm2 ) * fkm1 / fk;
    ak = 1.0 / ( ak * fkm1 );
    fmkm1 = b * ( da * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0 );

    if ( ioe == 0 )
    {
      sum = sum + fmkm2;
    }
    else
    {
      sum = sum + fmkm1;
    }
    ak = 1.0 / ( ak * fk );
    fk = fk + 2.0;
  }

  if ( ioe == 0 )
  {
    value = alnorm ( d, true ) + sum * g3;
  }
  else
  {
    value = alnorm ( drb, true ) + 2.0 * ( sum + tfn ( drb, a ) );
  }

  return value;
}
//****************************************************************************80

//****************************************************************************80

double tfn ( double x, double fx )

//****************************************************************************80
//
//  Purpose:
//
//    TFN calculates the T-function of Owen.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by JC Young, Christoph Minder.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    MA Porter, DJ Winstanley,
//    Remark AS R30:
//    A Remark on Algorithm AS76:
//    An Integral Useful in Calculating Noncentral T and Bivariate
//    Normal Probabilities,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, page 113.
//
//    JC Young, Christoph Minder,
//    Algorithm AS 76: 
//    An Algorithm Useful in Calculating Non-Central T and 
//    Bivariate Normal Distributions,
//    Applied Statistics,
//    Volume 23, Number 3, 1974, pages 455-457.
//
//  Parameters:
//
//    Input, double X, FX, the parameters of the function.
//
//    Output, double TFN, the value of the T-function.
//
{
# define NG 5

  double fxs;
  int i;
  double r[NG] = {
    0.1477621, 
    0.1346334, 
    0.1095432, 
    0.0747257, 
    0.0333357 };
  double r1;
  double r2;
  double rt;
  double tp = 0.159155;
  double tv1 = 1.0E-35;
  double tv2 = 15.0;
  double tv3 = 15.0;
  double tv4 = 1.0E-05;
  double u[NG] = {
    0.0744372, 
    0.2166977, 
    0.3397048,
    0.4325317, 
    0.4869533 };
  double value;
  double x1;
  double x2;
  double xs;
//
//  Test for X near zero.
//
  if ( fabs ( x ) < tv1 )
  {
    value = tp * atan ( fx );
    return value;
  }
//
//  Test for large values of abs(X).
//
  if ( tv2 < fabs ( x ) )
  {
    value = 0.0;
    return value;
  }
//
//  Test for FX near zero.
//
  if ( fabs ( fx ) < tv1 )
  {
    value = 0.0;
    return value;
  }
//
//  Test whether abs ( FX ) is so large that it must be truncated.
//
  xs = - 0.5 * x * x;
  x2 = fx;
  fxs = fx * fx;
//
//  Computation of truncation point by Newton iteration.
//
  if ( tv3 <= log ( 1.0 + fxs ) - xs * fxs )
  {
    x1 = 0.5 * fx;
    fxs = 0.25 * fxs;

    for ( ; ; )
    {
      rt = fxs + 1.0;

      x2 = x1 + ( xs * fxs + tv3 - log ( rt ) ) 
      / ( 2.0 * x1 * ( 1.0 / rt - xs ) );

      fxs = x2 * x2;

      if ( fabs ( x2 - x1 ) < tv4 )
      {
        break;
      }
      x1 = x2;
    }
  }
//
//  Gaussian quadrature.
//
  rt = 0.0;
  for ( i = 0; i < NG; i++ )
  {
    r1 = 1.0 + fxs * pow ( 0.5 + u[i], 2 );
    r2 = 1.0 + fxs * pow ( 0.5 - u[i], 2 );

    rt = rt + r[i] * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 );
  }

  value = rt * x2 * tp;

  return value;
# undef NG
}


double betain ( double x, double p, double q, double beta, int &ifault )

//****************************************************************************80
//
//  Purpose:
//
//    BETAIN computes the incomplete Beta function ratio.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    KL Majumder, GP Bhattacharjee,
//    Algorithm AS 63:
//    The incomplete Beta Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 409-411.
//
//  Parameters:
//
//    Input, double X, the argument, between 0 and 1.
//
//    Input, double P, Q, the parameters, which
//    must be positive.
//
//    Input, double BETA, the logarithm of the complete
//    beta function.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    nonzero, an error occurred.
//
//    Output, double BETAIN, the value of the incomplete
//    Beta function ratio.
//
{
  double acu = 0.1E-14;
  double ai;
  double cx;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double value;
  double xx;

  value = x;
  ifault = 0;
//
//  Check the input arguments.
//
  if ( p <= 0.0 || q <= 0.0 )
  {
    ifault = 1;
    return value;
  }

  if ( x < 0.0 || 1.0 < x )
  {
    ifault = 2;
    return value;
  }
//
//  Special cases.
//
  if ( x == 0.0 || x == 1.0 )
  {
    return value;
  }
//
//  Change tail if necessary and determine S.
//
  psq = p + q;
  cx = 1.0 - x;

  if ( p < psq * x )
  {
    xx = cx;
    cx = x;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    xx = x;
    pp = p;
    qq = q;
    indx = false;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = ( int ) ( qq + cx * psq );
//
//  Use the Soper reduction formula.
//
  rx = xx / cx;
  temp = qq - ai;
  if ( ns == 0 )
  {
    rx = xx;
  }

  for ( ; ; )
  {
    term = term * temp * rx / ( pp + ai );
    value = value + term;;
    temp = fabs ( term );

    if ( temp <= acu && temp <= acu * value )
    {
      value = value * exp ( pp * log ( xx )
      + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

      if ( indx )
      {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ai;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
  }

  return value;
}


//**************************************************************************
void cdftnc(int& which, double& p, double& q, double& t, double& df, double& pnonc, int& status, double& bound)
{
    //               Cumulative Distribution Function
    //                  Non-Central T distribution
    //
    //                               Function
    //
    //     Calculates any one parameter of the noncentral t distribution give
    //     values for the others.
    //
    //                               Arguments
    //
    //     WHICH --> Integer indicating which  argument
    //               values is to be calculated from the others.
    //               Legal range: 1..3
    //               iwhich = 1 : Calculate P and Q from T,DF,PNONC
    //               iwhich = 2 : Calculate T from P,Q,DF,PNONC
    //               iwhich = 3 : Calculate DF from P,Q,T
    //               iwhich = 4 : Calculate PNONC from P,Q,DF,T
    //                    INTEGER WHICH
    //
    //        P <--> The integral from -infinity to t of the noncentral t-den
    //              Input range: (0,1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: (0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //        T <--> Upper limit of integration of the noncentral t-density.
    //               Input range: ( -infinity, +infinity).
    //               Search range: [ -1E100, 1E100 ]
    //                    DOUBLE PRECISION T
    //
    //        DF <--> Degrees of freedom of the noncentral t-distribution.
    //                Input range: (0 , +infinity).
    //                Search range: [1e-100, 1E10]
    //                    DOUBLE PRECISION DF
    //
    //     PNONC <--> Noncentrality parameter of the noncentral t-distributio
    //                Input range: [-1e6, 1E6].
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //                                Method
    //
    //     Upper tail    of  the  cumulative  noncentral t is calculated usin
    //     formulae  from page 532  of Johnson, Kotz,  Balakrishnan, Coninuou
    //     Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //***********************************************************************
    double tent6;
    tent6 = 1.0E6;
    double tol;
    tol = 1.0E-8;
    double atol;
    atol = 1.0E-50;
    double zero, one, inf, inff, tent66;
    zero = 1.0E-100;
    one = 1.0E0 - 1.0E-16;
    inf = 1.0E100;
    inff = -inf;
    double ccum, cum, fx;
    int ifault;
    int dff;
    static unsigned long qhi, qleft;
    double a1, a2, a3,a4;

    dff = int(df);

    if (t > inf) t = inf;
    if (t < -inf) t = -inf;
    if (df > 1.0E10) df = 1.0E10;
    if (t != t) {
        status = -4;
        return;
    }

    if (which != 4) {
        if (pnonc < -tent6) {
            status = -6; bound = -tent6;
            return;
        }
        if (pnonc > tent6) {
            status = -6; bound = tent6;
            return;
         }
    }
    
    if (which>=1 || which<=4) goto g30;
    if (which>=1) goto g10;
    bound = 1.0E0;
    goto g20;
g10:
    bound = 5.0E0;
g20:
    status = -1;
    return;
g30:
    if (which == 1) goto g70;
    if (p>0 && p<one) goto g60;
    if (p<=0.0) goto g40;
    bound=0.;
    goto g50;
g40:
    bound = one;
g50:
    status=-2;
    return;
g60:
g70:
    if (which == 3) goto g90;
    if (df > 0.) goto g80;
    bound = 0.0;
    status = -5;
    return;
g80:
g90:
    if (which == 4) goto g100;
g100:
    if (which==1)  {
        p=prncst(t,dff,pnonc,ifault);
        q=1-p;
        status = 0;
    } 
    if (which==2) {
        t = 5.0;
        a1 = 0.5; a2 = 0.5; a3 = 5.0;
        dstinv(&inff, &inf, &a1, &a2, &a3, &atol, &tol);
        status = 0;
        dinvr(&status, &t, &fx, &qleft, &qhi);
    g110:
        if (status != 1) goto g120;
        cum = prncst(t, dff, pnonc, ifault);
        ccum = 1 - cum;
        fx = cum - p;
        dinvr(&status, &t, &fx, &qleft, &qhi);
        goto g110;
    g120:
        if (status != -1) goto g150;
        if (qleft == 1) goto g130;
        status = 1;
        bound = -inf;
        goto g140;
    g130:
        status = 2;
        bound = inf;
    g140:
    g150:
        a4 = 0;
   }
    if(which==3) {
        df = 5.0E0;
        a1 = 0.5; a2 = 0.5; a3 = 5.0;
        dstinv(&zero,&tent6,&a1,&a2,&a3,&atol,&tol);
        status = 0;
        dinvr(&status,&df,&fx,&qleft,&qhi);
g160:
        if (status !=1) goto g170;
        cum=prncst (t,dff,pnonc,ifault);
        ccum=1-cum;
        fx = cum - p;
        dinvr(&status,&df,&fx,&qleft,&qhi);
        goto g160;
g170:
        if (status!=-1) goto g200;
        if (qleft==1) goto g180;
        status = 1;
        bound = zero;
        goto g190;
g180:
        status = 2;
        bound = inf;
g190:
    g200:
        a4 = 0;
    } 
    if (which==4) { 
        pnonc = 5.0E0;
        tent66 = -tent6;
        a1 = 0.5; a2 = 0.5; a3 = 5.0;
        dstinv(&tent66,&tent6,&a1,&a2,&a3,&atol,&tol);
        status = 0;
        dinvr(&status,&pnonc,&fx,&qleft,&qhi);
g210:
        if (status!=1) goto g220;
        cum=prncst(t,dff,pnonc,ifault);
        ccum=1-cum;
        fx = cum - p;
        dinvr(&status,&pnonc,&fx,&qleft,&qhi);
        goto g210;
g220:
        if (status!=-1) goto g250;
        if (qleft==1) goto g230;
        status = 1;
        bound = 0.0E0;
        goto g240;
g230:
        status = 2;
        bound = tent6;
g240:
    g250:
        a4 = 0;
    }
    return;
}



//#########################################################
 void myswap(double& ax, double& ay) {
        double s = ax; ax = ay; ay = s;
    }
//#####################################################
double binom(double a, double b, double x) {
	//Ip(b-a,a+1,0.5)
	double z;
	z = incompletebeta(b - a, a + 1, x);
	return z;
}




//#########################################################
void cum(int n,double* x,int* r,int &k,double* fcum,double* ycum) {

//Kaplan-Meier

 int i,j;
 double s;
 
  k=0;
  for (i=1;i<=n;i++) {
      s = 1.0;
      for (j=1;j<=i;j++)  if (r[j-1]==0) s=s*(n*1.0-j*1.0)/(n*1.0- j*1.0+1.0);
         if (r[i-1]==0) {
          //fcum[k]= 1.-s;
          //fcum[k]=((1.-s)*n-0.375)/(n+0.25);
          //fcum[k]=((1.-s)*n-0.5)/n;
          fcum[k]=(1.-s)*n*1.0/(n+1.0);
          if (s==1 || s==0) fcum[k] = ((1-s)*n-0.375)/(n+0.25);
          ycum[k]=x[i-1];k++;
      }
  }
}
//###########################################################
double invnormaldistribution(double y0) {
	double result, x, y, z, y2, x0, x1, p0, q0, p1, q1, p2, q2, zz;
	int code;

	if (y0 <= 0) {
		result = -maxrealnumber;
		return result;
	}
	if (y0 >= 1) {
		result = maxrealnumber;
		return result;
	}
	code = 1;
	y = y0;
	if (y > 1.0 - expm2) {
		y = 1.0 - y;
		code = 0;
	}
	if (y > expm2) {
		y = y - 0.5;
		y2 = y * y;
		p0 = -59.9633501014107895267;
		p0 = 98.0010754185999661536 + y2 * p0;
		p0 = -56.6762857469070293439 + y2 * p0;
		p0 = 13.9312609387279679503 + y2 * p0;
		p0 = -1.23916583867381258016 + y2 * p0;
		q0 = 1;
		q0 = 1.95448858338141759834 + y2 * q0;
		q0 = 4.67627912898881538453 + y2 * q0;
		q0 = 86.3602421390890590575 + y2 * q0;
		q0 = -225.462687854119370527 + y2 * q0;
		q0 = 200.260212380060660359 + y2 * q0;
		q0 = -82.0372256168333339912 + y2 * q0;
		q0 = 15.9056225126211695515 + y2 * q0;
		q0 = -1.18331621121330003142 + y2 * q0;
		x = y + y * y2 * p0 / q0;
		x = x * s2pi;
		result = x;
		return result;
	}

	zz = log(y);
	x = sqrt(-(2.0 * zz));
	x0 = x - log(x) / x;
	z = 1.0 / x;

	if (x < 8.0) {
		p1 = 4.05544892305962419923;
		p1 = 31.5251094599893866154 + z * p1;
		p1 = 57.1628192246421288162 + z * p1;
		p1 = 44.0805073893200834700 + z * p1;
		p1 = 14.6849561928858024014 + z * p1;
		p1 = 2.18663306850790267539 + z * p1;
		p1 = -(1.40256079171354495875 * 0.1) + z * p1;
		p1 = -(3.50424626827848203418 * 0.01) + z * p1;
		p1 = -(8.57456785154685413611 * 0.0001) + z * p1;
		q1 = 1;
		q1 = 15.7799883256466749731 + z * q1;
		q1 = 45.3907635128879210584 + z * q1;
		q1 = 41.3172038254672030440 + z * q1;
		q1 = 15.0425385692907503408 + z * q1;
		q1 = 2.50464946208309415979 + z * q1;
		q1 = -(1.42182922854787788574 * 0.1) + z * q1;
		q1 = -(3.80806407691578277194 * 0.01) + z * q1;
		q1 = -(9.33259480895457427372 * 0.0001) + z * q1;
		x1 = z * p1 / q1;
	}
	else {
		p2 = 3.23774891776946035970;
		p2 = 6.91522889068984211695 + z * p2;
		p2 = 3.93881025292474443415 + z * p2;
		p2 = 1.33303460815807542389 + z * p2;
		p2 = 2.01485389549179081538 * 0.1 + z * p2;
		p2 = 1.23716634817820021358 * 0.01 + z * p2;
		p2 = 3.01581553508235416007 * 0.0001 + z * p2;
		p2 = 2.65806974686737550832 * 0.000001 + z * p2;
		p2 = 6.23974539184983293730 * 0.000000001 + z * p2;
		q2 = 1;
		q2 = 6.02427039364742014255 + z * q2;
		q2 = 3.67983563856160859403 + z * q2;
		q2 = 1.37702099489081330271 + z * q2;
		q2 = 2.16236993594496635890 * 0.1 + z * q2;
		q2 = 1.34204006088543189037 * 0.01 + z * q2;
		q2 = 3.28014464682127739104 * 0.0001 + z * q2;
		q2 = 2.89247864745380683936 * 0.000001 + z * q2;
		q2 = 6.79019408009981274425 * 0.000000001 + z * q2;
		x1 = z * p2 / q2;
	}
	x = x0 - x1;
	if (code != 0) x = -x;
	result = x;
	return result;
}
//######################################################
double normaldistribution(double x) {
	double p;
	p = 0.5 * (erf(x / 1.41421356237309504880) + 1.0);
	return p;
}
//*****************Приближенные доверительные граница для квантиля**************************
void lmtaprn(double beta,int n,double t1,double t2,double t12,double d,double* xp){

double zb,f1x,f2x,f4x,e11,e2x,e3x,e44;

zb=invnormaldistribution(beta);
f1x=t2/(n-1.0);
f2x=2*t12/sqrt(n);
f4x=1.0-f1x/2;
e3x=f4x*f4x-zb*zb*f1x;
e11=f4x*d+zb*zb*f2x/2;
e2x=d*d-zb*zb*t1;
e44=sqrt(abs(e11*e11-e2x*e3x));
xp[0]=(e11-e44)/e3x;
xp[1]=(e11+e44)/e3x;
}
//************************************************************************
double invnontap(double beta,int f,double d) {
	double zb, z, f4x;
 zb=invnormaldistribution(beta);
 f4x=1.0-1.0/(4.0*f);
 z=(f4x*d+zb*sqrt(f4x*f4x-zb*zb/(2.*f)+d*d/(2.*f)))/(f4x*f4x-zb*zb/(2.*f));
 return z;
}
//*****************************************************************************
double incompletebetaps(double a, double b, double x, double maxgam) {
			   double ai, u, v, t1, t, s, z;
			   int n;
			   ai = 1 / a;
			   u = (1 - b) * x;
			   v = u / (a + 1);
			   t1 = v;
			   t = u;
			   n = 2;
			   s = 0;
			   z = machineepsilon * ai;
			   while (abs(v) > z) {
				   u = (n - b) * x / n;
				   t = t * u; v = t / (a + n);
				   v = v / pow(10, 25);
				   s = s + v;
				   n++;
				   v = v * pow(10, 25);
			   }
			   s = s * pow(10, 25);
			   s = s + t1;
			   s = s + ai;
			   u = a * log(x);
			   if (((a + b) < maxgam) && (abs(u) < log(maxrealnumber))) {
				   t = tgamma(a + b) / (tgamma(a) * tgamma(b));
				   s = s * t * pow(x, a);
			   }
			   else {
				   t = lgamma(a + b) - lgamma(a) - lgamma(b) + u + log(s);
				   if (t < log(minrealnumber)) {
					   s = 0;
				   }
				   else {
					   s = exp(t);
				   }
			   }
			   return s;
		   }
//*****************************************************************************
double incompletebetafe(double a, double b, double x, double big, double biginv) {
			   double k1, k2, k3, k4, k5, k6, k7, k8;
			   k1 = a; k2 = a + b; k3 = a; k4 = a + 1; k5 = 1; k6 = b - 1; k7 = k4; k8 = a + 2;
			   double pkm2, qkm2, pkm1, qkm1, ans, r, thresh, xk, pk, qk, t, result;
			   int n;
			   pkm2 = 0; qkm2 = 1; pkm1 = 1; qkm1 = 1; ans = 1; r = 1; n = 0;
			   thresh = 3 * machineepsilon;
			   do {
				   xk = -(x * k1 * k2 / (k3 * k4));
				   pk = pkm1 + pkm2 * xk;
				   qk = qkm1 + qkm2 * xk;
				   pkm2 = pkm1;
				   pkm1 = pk;
				   qkm2 = qkm1;
				   qkm1 = qk;
				   xk = x * k5 * k6 / (k7 * k8);
				   pk = pkm1 + pkm2 * xk;
				   qk = qkm1 + qkm2 * xk;
				   pkm2 = pkm1;
				   pkm1 = pk;
				   qkm2 = qkm1;
				   qkm1 = qk;
				   if (qk != 0) {
					   r = pk / qk;
				   }
				   if (r != 0) {
					   t = abs((ans - r) / r);
					   ans = r;
				   }
				   else {
					   t = 1;
				   }
				   if (t < thresh) {
					   break;
				   }
				   k1 = k1 + 1;
				   k2 = k2 + 1;
				   k3 = k3 + 2;
				   k4 = k4 + 2;
				   k5 = k5 + 1;
				   k6 = k6 - 1;
				   k7 = k7 + 2;
				   k8 = k8 + 2;
				   if ((abs(qk) + abs(pk)) > big) {
					   pkm2 *= biginv;
					   pkm1 *= biginv;
					   qkm2 *= biginv;
					   qkm1 *= biginv;
				   }
				   if ((abs(qk) < biginv) || (abs(pk) < biginv)) {
					   pkm2 *= big;
					   pkm1 *= big;
					   qkm2 *= big;
					   qkm1 *= big;
				   }
				   n++;
			   } while (n < 300);
			   result = ans;
			   return result;
		   }
//*****************************************************************************
double incompletebetafe2(double a, double b, double x, double big, double biginv) {
			   double k1, k2, k3, k4, k5, k6, k7, k8, z;
			   double pkm2, qkm2, pkm1, qkm1, ans, r, thresh, xk, pk, qk, t, result;
			   int n;
			   k1 = a; k2 = b - 1; k3 = a; k4 = a + 1; k5 = 1; k6 = a + b; k7 = a + 1; k8 = a + 2;
			   pkm2 = 0; qkm2 = 1; pkm1 = 1; qkm1 = 1;
			   z = x / (1 - x);
			   ans = 1; r = 1; n = 0;
			   thresh = 3 * machineepsilon;
			   do {
				   xk = -(z * k1 * k2 / (k3 * k4));
				   pk = pkm1 + pkm2 * xk;
				   qk = qkm1 + qkm2 * xk;
				   pkm2 = pkm1;
				   pkm1 = pk;
				   qkm2 = qkm1;
				   qkm1 = qk;
				   xk = z * k5 * k6 / (k7 * k8);
				   pk = pkm1 + pkm2 * xk;
				   qk = qkm1 + qkm2 * xk;
				   pkm2 = pkm1;
				   pkm1 = pk;
				   qkm2 = qkm1;
				   qkm1 = qk;
				   if (qk != 0) {
					   r = pk / qk;
				   }
				   if (r != 0) {
					   t = abs((ans - r) / r);
					   ans = r;
				   }
				   else {
					   t = 1;
				   }
				   if (t < thresh) {
					   break;
				   }
				   k1++;
				   k2--;
				   k3 = k3 + 2;
				   k4 = k4 + 2;
				   k5++;
				   k6++;
				   k7 = k7 + 2;
				   k8 = k8 + 2;
				   if ((abs(qk) + abs(pk)) > big) {
					   pkm2 *= biginv;
					   pkm1 *= biginv;
					   qkm2 *= biginv;
					   qkm1 *= biginv;
				   }
				   if ((abs(qk) < biginv) || (abs(pk) < biginv)) {
					   pkm2 *= big;
					   pkm1 *= big;
					   qkm2 *= big;
					   qkm1 *= big;
				   }
				   n++;
			   } while (n < 300);
			   result = ans;
			   return result;
		   }
//*****************************************************************************
double incompletebeta(double a, double b, double x) {

			   double big, biginv, maxgam, minlog, maxlog, result, w, xc, t, y;
			   int flag;

			   big = 4.5035996273705 * pow(10, 15);
			   biginv = 2.22044604925031 * pow(10, -16);
			   maxgam = 171.624376956303;
			   minlog = log(minrealnumber);
			   maxlog = log(maxrealnumber);
			   if (x == 0) {
				   result = 0;
				   return result;
			   }
			   if (x == 1) {
				   result = 1;
				   return result;
			   }
			   flag = 0;
			   if ((b * x <= 1) && (x <= 0.95)) {
				   result = incompletebetaps(a, b, x, maxgam);
				   return result;
			   }
			   w = 1 - x;
			   if (x > (a / (a + b))) {
				   flag = 1;
				   t = a;
				   a = b;
				   b = t;
				   xc = x;
				   x = w;
			   }
			   else {
				   xc = w;
			   }
			   if ((flag == 1) && ((b * x <= 1)) && (x <= 0.95)) {
				   t = incompletebetaps(a, b, x, maxgam);
				   if (t <= machineepsilon) {
					   result = 1 - machineepsilon;
				   }
				   else {
					   result = 1 - t;
				   }
				   return result;
			   }
			   y = x * (a + b - 2) - (a - 1);
			   if (y < 0) {
				   w = incompletebetafe(a, b, x, big, biginv);
			   }
			   else {
				   w = incompletebetafe2(a, b, x, big, biginv) / xc;
			   }
			   y = a * log(x);
			   t = b * log(xc);
			   if (((a + b) < maxgam) && (abs(y) < maxlog) && (abs(t) < maxlog)) {
				   t = pow(xc, b);
				   t *= pow(x, a);
				   t /= a;
				   t *= w;
				   t *= (tgamma(a + b) / (tgamma(a) * tgamma(b)));
				   if (flag == 1) {
					   if (t <= machineepsilon) {
						   result = 1 - machineepsilon;
					   }
					   else {
						   result = 1 - t;
					   }
				   }
				   else {
					   result = t;
				   }
				   return result;
			   }
			   y += t + lgamma(a + b) - lgamma(a) - lgamma(b);
			   y += log(w / a);
			   if (y < minlog) {
				   t = 0;
			   }
			   else {
				   t = exp(y);
			   }
			   if (flag == 1) {
				   if (t <= machineepsilon) {
					   t = 1 - machineepsilon;
				   }
				   else {
					   t = 1 - t;
				   }
			   }
			   result = t;
			   return result;
		   }
//*****************************************************************************
double fdistribution(double a, double b, double x) {
			   //w=(a*x)/(b+a*x);
			   double result;
			   result = incompletebeta(0.5 * a, 0.5 * b, a * x / (b + a * x));
			   return result;
		   }
//*****************************************************************************
double studenttdistr(int k, double t) {
			   double result, x, z, xsqk, p, tz, f;
			   int rk, j;
			   if (t == 0) {
				   result = 0.5;
				   return result;
			   }
			   if (t < 0) {
				   x = -t;
			   }
			   else
			   {
				   x = t;
			   }
			   rk = k;
			   z = 1 + x * x / rk;
			   if ((k % 2) != 0) {
				   xsqk = x / sqrt(rk);
				   p = atan(xsqk);
				   if (k > 1) {
					   f = 1;
					   tz = 1;
					   j = 3;
					   while ((j <= (k - 2)) && ((tz / f) > machineepsilon)) {
						   tz *= (j - 1) / (z * j);
						   f += tz;
						   j += 2;
					   }
					   p += f * xsqk / z;
				   }
				   p *= 2 / pi;
			   }
			   else {
				   f = 1; tz = 1; j = 2;
				   while (j <= (k - 2) && ((tz / f) > machineepsilon)) {
					   tz *= (j - 1) / (z * j);
					   f += tz;
					   j += 2;
				   }
				   p = f * x / sqrt(z * rk);
			   }
			   if (t < 0) {
				   result = 0.5 - 0.5 * p;
				   return result;
			   }
			   result = 0.5 + 0.5 * p;
			   return result;
		   }
//*****************************************************************************
double invincompletebeta(double a, double b, double y) {
			   double result, aaa, bbb, y0, d, yyy, x, x0, x1, lgm, yp, di, dithresh, yl, yh, xt, dir;
			   int  newt, newtcycle, breaknewtcycle, breakihalvecycle;
			   int i, rflg, nflg, mainlooppos, ihalve, ihalvecycle;

			   result = 0; aaa = 0;  bbb = 0; y0 = 0; d = 0; yyy = 0; x = 0;
			   x0 = 0; x1 = 0; lgm = 0; yp = 0; di = 0; dithresh = 0; yl = 0; yh = 0;
			   xt = 0; i = 0; rflg = 0; dir = 0; nflg = 0;
			   mainlooppos = 0; ihalve = 0; ihalvecycle = 0;
			   newt = 0; newtcycle = 0; breaknewtcycle = 0; breakihalvecycle = 0;
			   i = 0;
			   if (y == 0)
			   {
				   result = 0;
				   return result;
			   }
			   if (y == 1.0)
			   {
				   result = 1;
				   return result;
			   }
			   dithresh = 0;
			   rflg = 0;
			   aaa = 0;
			   bbb = 0;
			   y0 = 0;
			   x = 0;
			   yyy = 0;
			   lgm = 0;
			   dir = 0;
			   di = 0;
			   x0 = 0;
			   yl = 0;
			   x1 = 1.0;
			   yh = 1.0;
			   nflg = 0;
			   mainlooppos = 0;
			   ihalve = 1;
			   ihalvecycle = 2;
			   newt = 3;
			   newtcycle = 4;
			   breaknewtcycle = 5;
			   breakihalvecycle = 6;
			   while (true) {
				   if (mainlooppos == 0)
				   {
					   if (a <= 1.0 || b <= 1.0)
					   {
						   dithresh = 0.000001;
						   rflg = 0;
						   aaa = a;
						   bbb = b;
						   y0 = y;
						   x = aaa / (aaa + bbb);
						   yyy = incompletebeta(aaa, bbb, x);
						   mainlooppos = ihalve;
						   continue;
					   }
					   else
					   {
						   dithresh = 0.0001;
					   }
					   yp = -invnormaldistribution(y);
					   if (y > 0.5)
					   {
						   rflg = 1;
						   aaa = b;
						   bbb = a;
						   y0 = 1.0 - y;
						   yp = -yp;
					   }
					   else
					   {
						   rflg = 0;
						   aaa = a;
						   bbb = b;
						   y0 = y;
					   }
					   lgm = (yp * yp - 3.0) / 6.0;
					   x = 2.0 / (1.0 / (2.0 * aaa - 1.0) + 1.0 / (2.0 * bbb - 1.0));
					   d = yp * sqrt(x + lgm) / x - (1.0 / (2.0 * bbb - 1.0) - 1.0 / (2.0 * aaa - 1.0)) * (lgm + 5.0 / 6.0 - 2.0 / (3.0 * x));
					   d = 2.0 * d;
					   if (d < log(minrealnumber))
					   {
						   x = 0;
						   break;
					   }
					   x = aaa / (aaa + bbb * exp(d));
					   yyy = incompletebeta(aaa, bbb, x);
					   yp = (yyy - y0) / y0;
					   if (abs(yp) < 0.2)
					   {
						   mainlooppos = newt;
						   continue;
					   }
					   mainlooppos = ihalve;
					   continue;
				   }
				   if (mainlooppos == ihalve)
				   {
					   dir = 0;
					   di = 0.5;
					   i = 0;
					   mainlooppos = ihalvecycle;
					   continue;
				   }
				   if (mainlooppos == ihalvecycle)
				   {
					   if (i <= 99)
					   {
						   if (i != 0)
						   {
							   x = x0 + di * (x1 - x0);
							   if (x == 1.0)
							   {
								   x = 1.0 - machineepsilon;
							   }
							   if (x == 0)
							   {
								   di = 0.5;
								   x = x0 + di * (x1 - x0);
								   if (x == 0)
								   {
									   break;
								   }
							   }
							   yyy = incompletebeta(aaa, bbb, x);
							   yp = (x1 - x0) / (x1 + x0);
							   if (abs(yp) < dithresh)
							   {
								   mainlooppos = newt;
								   continue;
							   }
							   yp = (yyy - y0) / y0;
							   if (abs(yp) < dithresh)
							   {
								   mainlooppos = newt;
								   continue;
							   }
						   }
						   if (yyy < y0)
						   {
							   x0 = x;
							   yl = yyy;
							   if (dir < 0)
							   {
								   dir = 0;
								   di = 0.5;
							   }
							   else
							   {
								   if (dir > 3)
								   {
									   di = 1.0 - (1.0 - di) * (1.0 - di);
								   }
								   else
								   {
									   if (dir > 1)
									   {
										   di = 0.5 * di + 0.5;
									   }
									   else
									   {
										   di = (y0 - yyy) / (yh - yl);
									   }
								   }
							   }
							   dir = dir + 1;
							   if (x0 > 0.75)
							   {
								   if (rflg == 1)
								   {
									   rflg = 0;
									   aaa = a;
									   bbb = b;
									   y0 = y;
								   }
								   else
								   {
									   rflg = 1;
									   aaa = b;
									   bbb = a;
									   y0 = 1.0 - y;
								   }
								   x = 1.0 - x;
								   yyy = incompletebeta(aaa, bbb, x);
								   x0 = 0;
								   yl = 0;
								   x1 = 1.0;
								   yh = 1.0;
								   mainlooppos = ihalve;
								   continue;
							   }
						   }
						   else
						   {
							   x1 = x;
							   if (rflg == 1 && x1 < machineepsilon)
							   {
								   x = 0;
								   break;
							   }
							   yh = yyy;
							   if (dir > 0)
							   {
								   dir = 0;
								   di = 0.5;
							   }
							   else
							   {
								   if (dir < -3)
								   {
									   di = di * di;
								   }
								   else
								   {
									   if (dir < -1)
									   {
										   di = 0.5 * di;
									   }
									   else
									   {
										   di = (yyy - y0) / (yh - yl);
									   }
								   }
							   }
							   dir = dir - 1;
						   }
						   i = i + 1;
						   mainlooppos = ihalvecycle;
						   continue;
					   }
					   else
					   {
						   mainlooppos = breakihalvecycle;
						   continue;
					   }
				   }
				   if (mainlooppos == breakihalvecycle)
				   {
					   if (x0 >= 1.0)
					   {
						   x = 1.0 - machineepsilon;
						   break;
					   }
					   if (x <= 0)
					   {
						   x = 0;
						   break;
					   }
					   mainlooppos = newt;
					   continue;
				   }
				   if (mainlooppos == newt)
				   {
					   if (nflg != 0)
					   {
						   break;
					   }
					   nflg = 1;
					   lgm = lgamma(aaa + bbb) - lgamma(aaa) - lgamma(bbb);
					   i = 0;
					   mainlooppos = newtcycle;
					   continue;
				   }
				   if (mainlooppos == newtcycle)
				   {
					   if (i <= 7)
					   {
						   if (i != 0)
						   {
							   yyy = incompletebeta(aaa, bbb, x);
						   }
						   if (yyy < yl)
						   {
							   x = x0;
							   yyy = yl;
						   }
						   else
						   {
							   if (yyy > yh)
							   {
								   x = x1;
								   yyy = yh;
							   }
							   else
							   {
								   if (yyy < y0)
								   {
									   x0 = x;
									   yl = yyy;
								   }
								   else
								   {
									   x1 = x;
									   yh = yyy;
								   }
							   }
						   }
						   if (x == 1.0 || x == 0)
						   {
							   mainlooppos = breaknewtcycle;
							   continue;
						   }
						   d = (aaa - 1.0) * log(x) + (bbb - 1.0) * log(1.0 - x) + lgm;
						   if (d < log(minrealnumber))
						   {
							   break;
						   }
						   if (d > log(maxrealnumber))
						   {
							   mainlooppos = breaknewtcycle;
							   continue;
						   }
						   d = exp(d);
						   d = (yyy - y0) / d;
						   xt = x - d;
						   if (xt <= x0)
						   {
							   yyy = (x - x0) / (x1 - x0);
							   xt = x0 + 0.5 * yyy * (x - x0);
							   if (xt <= 0)
							   {
								   mainlooppos = breaknewtcycle;
								   continue;
							   }
						   }
						   if (xt >= x1)
						   {
							   yyy = (x1 - x) / (x1 - x0);
							   xt = x1 - 0.5 * yyy * (x1 - x);
							   if (xt >= 1)
							   {
								   mainlooppos = breaknewtcycle;
								   continue;
							   }
						   }
						   x = xt;
						   if (abs(d / x) < 128.0 * machineepsilon)
						   {
							   break;
						   }
						   i = i + 1;
						   mainlooppos = newtcycle;
						   continue;
					   }
					   else
					   {
						   mainlooppos = breaknewtcycle;
						   continue;
					   }
				   }
				   if (mainlooppos == breaknewtcycle)
				   {
					   dithresh = 256.0 * machineepsilon;
					   mainlooppos = ihalve;
					   continue;
				   }
			   }
			   if (rflg != 0)
			   {
				   if (x <= machineepsilon)
				   {
					   x = 1.0 - machineepsilon;
				   }
				   else
				   {
					   x = 1.0 - x;
				   }
			   }
			   result = x;
			   return result;
		   }
 //*************************************************************************
		   double invstudenttdistr(int k, double p) {
			   double t, z;
			   if (p == 0.5) {
				   t = 0;
				   return t;
			   }
			   z = 1 - 2 * p;
			   z = invincompletebeta(0.5, 0.5 * k, abs(z));
			   t = sqrt(k * z / (1 - z));
			   if (p < 0.5) t = -t;
			   return t;
		   }
//*************************************************************************
		   double invfdistribution(int a, int  b, double y) {
			   double result, w;
			   result = 0; w = 0;
			   w = incompletebeta(0.5 * b, 0.5 * a, 0.5);
			   if ((w > y) || (y < 0.001)) {
				   w = invincompletebeta(0.5 * b, 0.5 * a, y);
				   result = (b - b * w) / (a * w);
			   }
			   else {
				   w = invincompletebeta(0.5 * a, 0.5 * b, 1.0 - y);
				   result = b * w / (a * (1.0 - w));
			   }
			   return result;
		   }
		   //*****************************************************************************
		   int minint(int m1, int m2) {
			   if (m1 < m2) {
				   return m1;
			   }
			   else {
				   return m2;
			   }
		   }
		   //*****************************************************************************
		   int maxint(int m1, int m2) {
			   if (m1 > m2) {
				   return m1;
			   }
			   else {
				   return m2;
			   }
		   }
//##############################################################
		   double incompletegammac(double a, double x) {
			   double result, ax, c, ans, y, z, r,zx;
			   double pkm2, qkm2, pkm1, qkm1, yc, pk, qk, t;
			   if (x <= 0 || a <= 0) {
				   result = 1;
				   return result;
			   }
			   if (x < 1 || x < a)
			   {
				  zx= incompletegamma(a, x);
				  result = 1 - zx;
				   return result;
			   }
			   ax = a * log(x) - x - lgamma(a);
			   if (ax < -709.78271289338399)
			   {
				   result = 0;
				   return result;
			   }
			   ax = exp(ax);
			   y = 1 - a;
			   z = x + y + 1;
			   c = 0;
			   pkm2 = 1;
			   qkm2 = x;
			   pkm1 = x + 1;
			   qkm1 = z * x;
			   ans = pkm1 / qkm1;
			   do
			   {
				   c = c + 1;
				   y = y + 1;
				   z = z + 2;
				   yc = y * c;
				   pk = pkm1 * z - pkm2 * yc;
				   qk = qkm1 * z - qkm2 * yc;
				   if (qk != 0)
				   {
					   r = pk / qk;
					   t = abs((ans - r) / r);
					   ans = r;
				   }
				   else
				   {
					   t = 1;
				   }
				   pkm2 = pkm1;
				   pkm1 = pk;
				   qkm2 = qkm1;
				   qkm1 = qk;
				   if (abs(pk) > igammabignumber)
				   {
					   pkm2 = pkm2 * igammabignumberinv;
					   pkm1 = pkm1 * igammabignumberinv;
					   qkm2 = qkm2 * igammabignumberinv;
					   qkm1 = qkm1 * igammabignumberinv;
				   }
			   } while (t > igammaepsilon);
			   result = ans * ax;
			   return result;
		   }
//#########################################################
double incompletegamma(double a,double x) {
			   double result, ax, c, ans,r;
			   if (x <= 0 || a <= 0) {
				   result = 0;
				   return result;
			   }
			   if (x > 1 && x > a) {
				   result = 1 - incompletegammac(a, x);
				   return result;
			   }
			   ax = a * log(x) - x - lgamma(a);
			   if (ax < -709.782712893384) {
				   result = 0;
				   return result;
			   }
			   ax = exp(ax);
			   r = a; c = 1; ans = 1;
			   do {
				   r++; c *= x / r; ans += c;
			   } while ((c / ans) > igammaepsilon);
			   result = ans * ax / a;
			   return result;
		   }
//###########################################################
double invincompletegammac(double a, double y0) {
	double x0, yl, x1, dithresh, yh, d, y, x, lgm,result,dir;
	int i;
	x0 = igammabignumber; yl = 0;x1 = 0; yh = 1.0;
	dithresh = 5 * igammaepsilon;d = 1.0/(9.0 * a);
	 y = 1.0 - d - invnormaldistribution(y0) *sqrt(d);
	 x = a * y * y * y; lgm = lgamma(a);
	 i = 0;
	 while (i < 10)
			   {
				   if (x > x0 || x < x)
				   {
					   d = 0.0625;
					   break;
				   }
				   y = incompletegammac(a, x);
				   if (y<yl || y>yh)
				   {
					   d = 0.0625;
					   break;
				   }
				   if (y < y0)
				   {
					   x0 = x;
					   yl = y;
				   }
				   else
				   {
					   x1 = x;
					   yh = y;
				   }
				   d = (a - 1.0) * log(x) - x - lgm;
				   if (d < -709.78271289338399)
				   {
					   d = 0.0625;
					   break;
				   }
				   d = -exp(d);
				   d = (y - y0) / d;
				   if (abs(d / x) < igammaepsilon)
				   {
					   result = x;
					   return result;
				   }
				   x = x - d;
				   i = i + 1;
			   }
			   if (x0 == igammabignumber) {
				   if (x <= 0) x = 1;
				   while (x0 == igammabignumber) {
					   x = (1 + d) * x;
					   y = incompletegammac(a, x);
					   if (y < y0) {
						   x0 = x;
						   yl = y;
						   break;
					   }
					   d = d + d;
				   }
			   }
			   d = 0.5;
			   dir = 0;
			   i = 0;
			   while (i < 400) {
				   x = x1 + d * (x0 - x1);
				   y = incompletegammac(a, x);
				   lgm = (x0 - x1) / (x1 + x0);
				   if (abs(lgm) < dithresh) break;
				   lgm = (y - y0) / y0;
				   if (abs(lgm) < dithresh) break;
				   if (x <= 0)  break;
				   if (y >= y0) {
					   x1 = x;
					   yh = y;
					   if (dir < 0) {
						   dir = 0; d = 0.5;
					   }
					   else
					   {
						   if (dir > 1) {
							   d = 0.5 * d + 0.5;
						   }
						   else {
							   d = (y0 - yl) / (yh - yl);
						   }
					   }
					   dir = dir + 1;
				   }
				   else
				   {
					   x0 = x;
					   yl = y;
					   if (dir > 0)
					   {
						   dir = 0;
						   d = 0.5;
					   }
					   else
					   {
						   if (dir < -1)
						   {
							   d = 0.5 * d;
						   }
						   else
						   {
							   d = (y0 - yl) / (yh - yl);
						   }
					   }
					   dir = dir - 1;
				   }
				   i = i + 1;
			   }
			   result = x;
			   return result;
		   }
//*****************************************************************************
double invchisquaredistribution(double v,double y) {
		double result;
		result=2*invincompletegammac(0.5*v, y);
		return result;
		   }
//*************************************************************
double chisquaredistribution(double v, double x) {
	double result;
		result = incompletegamma(v / 2.0, x / 2.0);
		return result;
	}
//####################################################
void fatiq(int k, double* w, double* cx, double* y, double& pv, double& an, double& alpha, double& q, double& da, double& db) {
	int i, j, kn;
	double s1, s1x, s2, s3,bx, xx, pvx, xcp, a, b;

	double step;

	step = 0.1;
	kn = int(cx[k - 1] / step);
	if (pv != 0) kn = 1;
	q = 100000.;
	s2 = 0; s3 = 0;
	for (i = 0; i < k; i++) {
		s2 += w[i];
		s3 += w[i] * log(y[i]);
	}
	a = s3 / s2;
	for (j = 0; j < kn; j++) {
		pvx = j * (cx[k - 1] - step) / (kn - 1);
		if (pv != 0) pvx = pv;
		xcp = 0.;
		for (i = 0; i < k; i++)  xcp += w[i] * log(cx[i] - pvx);
		xcp = xcp / s2;
		s1 = 0.; s3 = 0.;
		for (i = 0; i < k; i++) {
			s1 += w[i] * pow(log(cx[i] - pvx) - xcp, 2);
			s3 += w[i] * (log(cx[i] - pvx) - xcp) * log(y[i]);
		}
		if (alpha != 0) b = -1. / alpha;
		if (alpha == 0) b = s3 / s1;
		s3 = 0;
		for (i = 0; i < k; i++) {
			s3 += w[i] * pow((log(y[i]) - a - b * (log(cx[i] - pvx) - xcp)), 2);
		}
		if (s3 < q) {
			q = s3; pv = pvx; bx = b; xx = xcp; s1x = s1;
		}
	}
	alpha = -1. / bx;
	an = exp(xx - a / bx);
	da = q / s2;
	db = q / s1x;
}
//###########################################################

/*
double mydistr(int dindex, int dtype, double df1, double df2, double p) {
    double z;
    if (dindex == 0) {
        //Normal distribution
        normal df(0, 1);
        if (dtype == 0) z = cdf(df, p);
        if (dtype == 1) z = quantile(df, p);
    }
    else if (dindex == 1) { //chisquare distribution
        chi_squared df(df1);
        if (dtype == 0) z = cdf(df, p);
        if (dtype == 1) z = quantile(df, p);
    }
    else if (dindex == 2) { //F distribution
        fisher_f df(df1, df2);
        if (dtype == 0) z = cdf(df, p);
        if (dtype == 1) z = quantile(df, p);
    }
    else if (dindex == 3) { // t distribution
        students_t df(df1);
        if (dtype == 0) z = cdf(df, p);
        if (dtype == 1) z = quantile(df, p);
    }
    return z;
}
double mydistrnc(int dindex, int dtype, double df1, double df2, double p,double dnc) {
    double z;
    if (dindex == 1) { // chisquare distribution
        non_central_chi_squared df(df1, dnc);
        if (dtype == 0) z = cdf(df, p);
        if (dtype == 1) z = quantile(df, p);
    }
    if (dindex == 2) { // F distribution
        non_central_f df(df1,df2, dnc);
        if (dtype == 0) z = cdf(df, p);
        if (dtype == 1) z = quantile(df, p);
    }
    if (dindex == 3) { // t distribution
        non_central_t df(df1, dnc);
        if (dtype == 0) z = cdf(df, p);
        if (dtype == 1) z = quantile(df, p);
    }
    return z;
}
double invnormaldistribution_1(double p) {
	double t2, t, x, d;
	if (p <= 0) { p = 0; return -5.; }
	if (p >= 1) { p = 1; return 5.; }
	d = p;
	if (d > 0.5) d = 1. - d;
	t2 = log(1. / (d * d)); t = sqrt(t2);
	x = t - (2.515517 + 0.802853 * t + 0.010328 * t2) / (1.0 + 1.432788 * t + 0.189269 * t2 + 0.001308 * t * t2);
	if (p <= 0.5) x = -x;
	d = 0.3989423 * exp(-x * x / 2);
	return x;
}
*/

//***********************************************************

void qsort_1(double *a, int lo, int hi) {
  int h, l;
  double p,t;

  if (lo < hi) {
    l = lo;
    h = hi;
    p = a[hi];

    do {
      while ((l < h) && (a[l] <= p)) 
          l = l+1;
      while ((h > l) && (a[h] >= p))
          h = h-1;
      if (l < h) {
          t = a[l];
          a[l] = a[h];
          a[h] = t;
      }
    } while (l < h);

    a[hi] = a[l];
    a[l] = p;

    qsort_1( a, lo, l-1 );
    qsort_1( a, l+1, hi );
  }
}


int Partition(double* ar, int lo, int hi) {
    double pivot = ar[lo];
    int i = lo - 1;
    int j = hi + 1;

    for (;;) {
        do
            i++;
        while (ar[i] < pivot);

        do
            j--;
        while (ar[j] > pivot);

        if (i >= j)
            return j;

        myswap(ar[i], ar[j]);
    }

    for (int j = lo; j < hi; ++j)
        if (ar[j] <= pivot)
            myswap(ar[++i], ar[j]);

    if (ar[hi] < ar[i + 1])
        myswap(ar[i + 1], ar[hi]);

    return i + 1;
}

void QuickSort(double* ar, int lo, int hi) {
    if (lo < hi) {
        int prt = Partition(ar, lo, hi);
        QuickSort(ar, lo, prt);
        QuickSort(ar, prt + 1, hi);
    }
}

void QuickSort1(double* ar, int hi) {
    if (hi > 0) {
        double pivot = ar[0];
        int i = -1;
        int j = hi + 1;

        for (;;) {
            do
                i++;
            while (ar[i] < pivot);

            do
                j--;
            while (ar[j] > pivot);

            if (i >= j) {
                i = j;
                goto LabelDown;
            }

           myswap(ar[i], ar[j]);
        }

        for (int j = 0; j < hi; ++j)
            if (ar[j] <= pivot)
                myswap(ar[++i], ar[j]);

        if (ar[hi] < ar[i + 1])
            myswap(ar[i + 1], ar[hi]);
        i++;

    LabelDown:
        QuickSort1(ar, i);
        QuickSort1(ar + i + 1, hi - i - 1);
    }
}
//******************************************************************
int compare(const void* a, const void* b) {
    const double* x = (double*)a;
    const double* y = (double*)b;
    if (*x > *y) return 1;
    else if (*x < *y) return -1;
    return 0;
}
//##############################################################

void ordern(int n, double pr, double ps, double& er, double& vrs) {
    double p, pr1, ps1, xr, xr1, xr2, xr3, xr4, xr5, xr6, dr, qr, qs;
    double xs1, xs2, xs3, xs4, xs5, xs6, ds, xs;
    double z1, z2, z3, z4, z5, z6, z7;

    p = 1;
    xr = invnormaldistribution(pr);
    xs = invnormaldistribution(ps);
    qr = 1. - pr;
    qs = 1. - ps;
    pr1 = pr * p; ps1 = ps * p;
    xr = invnormaldistribution(pr1);
    xs = invnormaldistribution(ps1);
    dr = spi * exp(-xr * xr / 2.);
    ds = spi * exp(-xs * xs / 2.);
    xr1 = p / dr; xr2 = xr * xr1 * xr1;
    xr3 = (2. * xr * xr + 1.) * pow((p / dr), 3);
    xr4 = (6. * xr * xr * xr + 7. * xr) * pow((p / dr), 4);
    xr5 = (24. * pow(xr, 4) + 46. * xr * xr + 7.) * pow((p / dr), 5);
    xr6 = (120. * pow(xr, 5) + 326. * xr * xr * xr + 127. * xr) * pow((p / dr), 6);
    xs1 = p / ds; xs2 = xs * xs1 * xs1;
    xs3 = (2. * xs * xs + 1.) * pow((p / ds), 3);
    xs4 = (6. * xs * xs * xs + 7. * xs) * pow((p / ds), 4);
    xs5 = (24. * pow(xs, 4) + 46. * xs * xs + 7.) * pow((p / ds), 5);
    xs6 = (120. * pow(xs, 5) + 326. * xs * xs * xs + 127. * xs) * pow(p / ds, 6);

    er = xr + pr * qr * xr2 / (2. * (n + 2.)) + pr * qr * ((qr - pr) * xr3 / 3. + pr * qr * xr4 / 8.) / pow((n + 2.), 2) + pr * qr * (-(qr - pr) * xr3 / 3. + (pow((qr - pr), 2) - pr * qr) * xr4 / 4. + qr * pr * (qr - pr) * xr5 / 6. + pow((qr * pr), 2) * xr6 / 48.) / pow((n + 2.), 3);

    z1 = (qr - pr) * xr2 * xs1 + (qs - ps) * xr1 * xs2 + pr * qr * xr3 * xs1 / 2. + ps * qs * xr1 * xs3 / 2. + pr * qs * xr2 * xs2 / 2.;
    z1 = z1 * pr * qs / pow((n + 2.), 2);
    z2 = -(qr - pr) * xr2 * xs1 - (qs - ps) * xr1 * xs2 + (pow((qr - pr), 2) - pr * qr) * xr3 * xs1;
    z3 = (pow((qs - ps), 2) - ps * qs) * xr1 * xs3 + (1.5 * (qr - pr) * (qs - ps) + 0.5 * ps * qr - 2. * pr * qs) * xr2 * xs2;
    z4 = (5. / 6.) * pr * qr * (qr - pr) * xr4 * xs1 + (5. / 6.) * ps * qs * (qs - ps) * xr1 * xs4 + (pr * qs * (qr - pr) + .5 * pr * qr * (qs - ps)) * xr3 * xs2;
    z5 = (pr * qs * (qs - ps) + 0.5 * ps * qs * (qr - pr)) * xr2 * xs3 + (1. / 8.) * pow((pr * qr), 2) * xr5 * xs1 + (1. / 8.) * pow((ps * qs), 2) * xr1 * xs5;
    z6 = 0.25 * pr * pr * qr * qs * xr4 * xs2 + 0.25 * pr * ps * qs * qs * xr2 * xs4 + (2. * (pr * pr * qs * qs) + 3. * pr * qr * ps * qs) * xr3 * xs3 / 12.;
    z7 = z2 + z3 + z4 + z5 + z6;
    vrs = z1 + pr * qs * z7 / pow((n + 2.), 3) + pr * qs * xr1 * xs1 / (n + 2.);
}
//#############################################################
void orderw(int n, double pr, double ps, double &er, double &vrs) {
    double  qr, qs,xr, xr1, xr2, xr3, xr4, xr5, xr55, xr6;
    double xs1, xs2, xs3, xs4, xs5;
    double z1, z2, z3, z4, z5, z6, z7, a1, b1, c1, d1;

    xr = log(log(1. / (1. - pr)));
    qr = 1. - pr;
    //xs = log(log(1. / (1. - ps)));
    qs = 1. - ps;
    xr1 = 1. / (log(1. / (1. - pr)) * (1. - pr));
    xr2 = xr1 * (1. / (1. - pr) - xr1);
    xr3 = xr2 * xr2 / xr1 + xr1 * (1. / pow((1. - pr), 2) - xr2);
    xr4 = (3. * xr1 * xr2 * xr3 - 2. * pow(xr2, 3)) / pow(xr1, 2) + xr1 * (2. / pow((1. - pr), 3) - xr3);
    xr55 = (-12. * xr1 * xr2 * xr2 * xr3 + 3. * xr1 * xr1 * xr3 * xr3 + 4. * xr1 * xr1 * xr2 * xr4 + 6. * pow(xr2, 4));
    xr5 = xr55 / pow(xr1, 3) + xr1 * (6. / pow(1. - pr, 4) - xr4);
    a1 = -12. * pow(xr2, 3) * xr3 - 12. * xr1 * (2. * xr2 * xr3 * xr3 + xr2 * xr2 * xr4);
    b1 = 6. * xr1 * xr2 * xr3 * xr3 + 6. * xr1 * xr1 * xr3 * xr4;
    c1 = 8. * xr1 * xr2 * xr2 * xr4 + 4. * xr1 * xr1 * (xr3 * xr4 + xr2 * xr5);
    d1 = 24. * pow(xr2, 3) * xr3;
    xr6 = (pow(xr1, 3) * (a1 + b1 + c1 + d1) - 3. * xr1 * xr1 * xr2 * xr55) / pow(xr1, 6) + xr2 * (6. / pow((1. - pr), 4) - xr4) + xr1 * (24. / pow((1. - pr), 5) - xr5);
    xs1 = 1. / (log(1. / (1. - ps)) * (1. - ps));
    xs2 = xs1 * (1. / (1. - ps) - xs1);
    xs3 = xs2 * xs2 / xs1 + xs1 * (1. / pow((1. - ps), 2) - xs2);
    xs4 = (3. * xs1 * xs2 * xs3 - 2. * pow(xs2, 3)) / (xs1 * xs1) + xs1 * (2. / pow((1. - ps), 3) - xs3);
    xs5 = (-12. * xs1 * xs2 * xs2 * xs3 + 3. * xs1 * xs1 * xs3 * xs3 + 4. * xs1 * xs1 * xs2 * xs4 + 6. * pow(xs2, 4)) / pow(xs1, 3) + xs1 * (6. / pow((1. - ps), 4) - xs4);

    z1 = pr * qr * xr2 / (2. * (n + 2.));
    z2 = pr * qr * ((qr - pr) * xr3 / 3. + pr * qr * xr4 / 8.) / pow((n + 2.), 2);
    z3 = -(qr - pr) * xr3 / 3. + (pow((qr - pr), 2) - pr * qr) * xr4 / 4.;
    z4 = qr * pr * (qr - pr) * xr5 / 6. + qr * qr * pr * pr * xr6 / 48.;
    er = xr + z1 + z2 + pr * qr * (z3 + z4) / pow((n + 2.), 3);

    z1 = (qr - pr) * xr2 * xs1 + (qs - ps) * xr1 * xs2 + pr * qr * xr3 * xs1 / 2. + ps * qs * xr1 * xs3 / 2. + pr * qs * xr2 * xs2 / 2.;
    z1 = z1 * pr * qs / pow((n + 2.), 2);
    z2 = -(qr - pr) * xr2 * xs1 - (qs - ps) * xr1 * xs2 + (pow((qr - pr), 2) - pr * qr) * xr3 * xs1;
    z3 = (pow((qs - ps), 2) - ps * qs) * xr1 * xs3 + (1.5 * (qr - pr) * (qs - ps) + 0.5 * ps * qr - 2. * pr * qs) * xr2 * xs2;
    z4 = (5. / 6.) * pr * qr * (qr - pr) * xr4 * xs1 + (5. / 6.) * ps * qs * (qs - ps) * xr1 * xs4 + (pr * qs * (qr - pr) + .5 * pr * qr * (qs - ps)) * xr3 * xs2;
    z5 = (pr * qs * (qs - ps) + .5 * ps * qs * (qr - pr)) * xr2 * xs3 + (1. / 8.) * pr * pr * qr * qr * xr5 * xs1 + (1. / 8.) * ps * ps * qs * qs * xr1 * xs5;
    z6 = 0.25 * pr * pr * qr * qs * xr4 * xs2 + 0.25 * pr * ps * qs * qs * xr2 * xs4 + (2. * pr * pr * qs * qs + 3. * pr * qr * ps * qs) * xr3 * xs3 / 12.;
    z7 = z2 + z3 + z4 + z5 + z6;
    vrs = z1 + pr * qs * z7 / pow((n + 2.), 3) + pr * qs * xr1 * xs1 / (n + 2.);
}
//#################################################################
 void randomsample(int k, int* r) {
		   int i,m,n;
		   if (k == 0) return;
		   n = pow(2, k);
		   srand(time(0));
		   for (i = 0; i <n; i++) r[i] = int(n * rand() / RAND_MAX) + 1;

		   m = -1;
		   while (m <n) {
		   m1: m++;
			   r[m] = int(n * rand() / RAND_MAX) + 1;
			   for (i = 0; i < n; i++) {
				   if (m != i) {
					   if (r[m] == r[i]) {
						   m--; goto m1;
					   }
				   }
			   }
		   }
   }
//#####################Матрица факторов x ##################################
void planmatrix(int k, int** x) {
		   int m,i, j, n,p,km;
		   int y[20];

		   n = pow(2, k);
		   for (i = 0; i < k+1; i++) y[i] =0;
//*********************************************************
		   i = 0; p = 0; km = 0;
		   while (p <= k) {
			   for (m = 1; m <= k; m++) x[km][m - 1] = y[m];  
			   i++; p = 1; j = i;
			   while (j % 2 == 0) {
				   j = j / 2; p++;
			   }
			   if (p <= k) y[p] =1 -y[p];
			   km++;
		   }
//********************************************
 }
//############MLE Normal Minimized Function############################
double fun1(double* xsimpl) {
    double s1,s2,s3,s4,z,psi,p,d,c1,c2;
    int i,kx;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;
    if (xsimpl[1]<=0) return 10000;
    for (i=0;i<sm.n;i++) {
            z = (sm.x[i] - xsimpl[0]) / xsimpl[1];
            d = spi * exp(-z * z / 2.);
            p = normaldistribution(z);
            psi = d / (1. - p);
            s1 +=(1.-sm.r[i])*(sm.x[i] - xsimpl[0]);
            s2 += (1.-sm.r[i])*pow(sm.x[i]-xsimpl[0],2);
            s3 += sm.r[i]*psi;
            s4 += sm.r[i]*psi*z;
            kx+=1-sm.r[i];
    }
    c1=s1+xsimpl[1]*s3;
    c2=s2+pow(xsimpl[1], 2)*(s4-kx);
    z=c1*c1+c2*c2;
    return z;
}
//#################MLE Weibull Minimized Function#######################
double fun2(double* xsimpl) {
 double s1,s2,s3,z,b,c;
 int i,k,n;
 n = sm.n;
   s1=0;s2=0;s3=0;k=0;
   b=xsimpl[0];
 for(i=0;i<n;i++) {
   k+=(1-sm.r[i]);
   s1+=pow(sm.x[i],b);
 }
   c=s1/k;

for(i=0;i<n;i++) {
  z=(pow(sm.x[i],b))/c;
  s3+=z*log(z);
  s2+=(1-sm.r[i])*log(z);
}
 c=s3-s2-k;
 return c*c;
}

//##########################################################

void CovMatrixMleW(int n,double*x,int* r,double c, double b, double **&v) {

    int i, k;
    double s1, s2, z, cpw, ckow;
    cpw = log(c); ckow = 1 / b; s1 = 0; s2 = 0; k = 0;
    for (i = 0; i < n; i++) {
        z = (log(x[i]) - cpw) / ckow;
        s1 += (1 -r[i]) * z;
        s2 += z * z * exp(z);
        k += (1 - r[i]);
    }
    v[0][0] =double(k)/double(n); v[0][1] = (k + s1) / n; v[1][0] = (k + s1) / n; v[1][1] = (k + s2) / n;
    v=InverseMatrix(v,2);
 
}
//############################################
void CovMatrixMleN(int n,double*x,int* r,double a, double s, double **&v) {

    double z, p, d, s1, s2, s3, psi;
    int j, k;
    s1 = 0; s2 = 0; s3 = 0; k = 0;

    for (j = 0; j < n; j++) {
        z = (x[j] - a) / s;
        p = normaldistribution(z);
        d = s2pi * exp(-x[j] * x[j] / 2);
        psi = d / (1 - p);
        s1 += r[j] * psi * (psi - z);
        s2 += r[j] * psi * z * (z * (psi - z) - 1);
        s3 += r[j] * psi * (z * (psi - z) - 1);
        k += (1 - r[j]);
    }

    v[0][0] = (k + s1) / n; v[0][1] = s3 / n;
    v[1][0] = s3/n; v[1][1] = (2 * k + s2) / n;
    v=InverseMatrix(v,2);
   
}
//#######################Nalder-Mead################################################
void simpl(double*xsimpl, double step, double eps, int &lim,int &ier,int nx,double &q,double(*fn)(double*)) {
    double x1[500], sum[25], xnx, step1, step2, difer, dif;
    double  alfa, beta, gama, suml, sums, sumh, sum2;
    int istep, k1, k2, k3, k4, vn, i, j, k, l, l1, l2, l3, kount;
    int ik, in, index;

    k = 0; q = 0;
    double dop[25];

    for (i=0;i<500;i++) x1[i]=0;
    for (i=0;i<25; i++) {dop[i]=0;sum[i]=0;}
        
    alfa = 1.0; beta = 0.45; gama = 2.8; istep = 0; ier = 0; k1 = nx + 1;
    k2 = nx + 2; k3 = nx + 3; k4 = nx + 4;
    vn = nx; xnx = 1. / (1.0*vn);
    step1=step/ (1.0*vn * sqrt(2.)) * (sqrt(vn + 1.) + vn - 1.);
    step2 = step / (1.0*vn * sqrt(2.)) * (sqrt(vn + 1.) - 1.);
    for (i = 2; i <= k1; i++) {
        l = (i - 1) * nx; l1 = l + i - 1;
        for (j = 1; j <= nx; j++) {
            l2 = l + j; x1[l2] = xsimpl[j-1] + step2;
        }
        x1[l1] = xsimpl[i-2]+step1;
    }
    for (j=1;j <= nx;j++) x1[j]=xsimpl[j-1];
    for (i=1;i<=k1;i++) {
        l=(i-1)*nx+1;
        for (ik=1;ik<=nx;ik++) {
            dop[ik-1]=x1[l];l++;
        }
        sum[i] = (fn(dop));
    }
a28:
    sumh = sum[1]; index = 1;
    for (i = 2;i<=k1;i++) {
        if (sum[i] >sumh) {
            sumh=sum[i];index=i;
        }
    }
    suml = sum[1]; kount = 1;
    for (i = 2; i <= k1; i++) {
        if (suml > sum[i]) {
            suml = sum[i]; kount = i;
        }
    }
    istep++; difer = 0.;
    for (i = 1; i <= k1; i++) {
        dif = sum[i] - sum[kount]; difer = difer + dif * dif;
    }
    difer=xnx*sqrt(difer);
    if (suml<=eps) goto a32;
    if (suml>eps) goto a31;
a32:
    if (difer<= eps) goto a30;
    if (difer>eps) goto a26;
a31:
    if ((difer / fabs(suml)) <= eps) goto a30;
    if ((difer / fabs(suml)) > eps) goto a26;
a30:
    difer = 0.; l1 = (kount - 1) * nx;
    for (i = 1; i <= k1; i++) {
        if (kount != i) {
            l = (i - 1) * nx;
            for (j = 1; j <= nx; j++) {
                l2 = l + j; l3 = l1 + j;
                if (x1[l3] <= eps) goto a35;
                if (x1[l3] > eps) goto a34;
            a34:
                dif = (x1[l2] - x1[l3]) / x1[l3];
                goto a33;
            a35:
                dif = x1[l2] - x1[l3];
            a33:
                difer = difer + dif * dif;
            }
        }
    }
    difer = xnx * sqrt(difer);
    if (difer <= eps) goto a23;
a26:
    if (istep >= lim) goto a38;
    if (k == 1) goto a17;
    for (j = 1; j <= nx; j++) {
        sum2 = 0.;
        for (i = 1; i <= k1; i++) {
            l = (i - 1) * nx + j; sum2 = sum2 + x1[l];
        }
        l2 = (index - 1) * nx + j; l1 = k1 * nx + j; x1[l1] = (sum2 - x1[l2]) * xnx;
        l3 = l1 + nx; x1[l3] = (1. + alfa) * x1[l1] - alfa * x1[l2];
    }
    in = k2 * nx + 1;
    for (ik = 1; ik <= nx; ik++) {
        dop[ik-1] = x1[in]; in++;
    }
    sum[k3] =(fn(dop));
    if (sum[k3] < suml) goto a11;
    sums = suml;
    for (i = 1; i <= k1; i++) {
        if (index == i) goto a12;
        if (sum[i] <= sums) goto a12;
        sums = sum[i];
    a12:
        continue;
    }
    if (sum[k3] > sums) goto a13;
    goto a14;
a11:
    for (j = 1; j <= nx; j++) {
        l = k1 * nx + j; l1 = l + nx; l2 = l1 + nx; x1[l2] = (1. - gama) * x1[l] + gama * x1[l1];
    }
    in = k3 * nx + 1;
    for (ik = 1; ik <= nx; ik++) {
        dop[ik-1] = x1[in]; in++;
    }
    sum[k4] =(fn(dop));
    if (sum[k4] < suml) goto a16;
    goto a14;
a13:
    if (sum[k3] > sumh) goto a17;
    k = 1; goto a14;
a17:
    for (j = 1; j <= nx; j++) {
        l = k3 * nx + j; l1 = (index - 1) * nx + j; l2 = l - nx - nx;
        x1[l] = beta * x1[l1] + (1. - beta) * x1[l2];
    }
    k = 0; in = k3 * nx + 1;
    for (ik = 1; ik <= nx; ik++) {
        dop[ik-1] = x1[in]; in++;
    }
    sum[k4] =(fn(dop));
    if (sumh > sum[k4]) goto a16;
    for (j = 1; j <= nx; j++) {
        l1 = (kount - 1) * nx + j;
        for (i = 1; i <= k1; i++) {
            if (i == kount) goto a20;
            l = (i - 1) * nx + j;
        a20:
            x1[l] = 0.5 * (x1[l] + x1[l1]);
        }
    }
    for (i = 1; i <= k1; i++) {
        if (i == kount) goto a29;
        l = (i - 1) * nx + 1;
        for (ik = 1; ik <= nx; ik++) {
            dop[ik-1] = x1[l]; l++;
        }
        sum[i]=(fn(dop));
    a29:
        continue;
    }
    goto a28;
a16:
    for (j = 1; j <= nx; j++) {
        l = (index - 1) * nx + j; l1 = k3 * nx + j; x1[l] = x1[l1];
        sum[index] = sum[k4];
    }
    goto a28;
a14:
    for (j = 1; j <= nx; j++) {
        l = (index - 1) * nx + j; l1 = k2 * nx + j; x1[l] = x1[l1];
        sum[index] = sum[k3];
    }
    goto a28;
a38:
    ier=1;
a23:
    for (j = 1; j <= nx; j++) {
        l = (kount - 1) * nx + j; xsimpl[j-1] = x1[l];
    }
    sum[1]=suml;lim=istep;q = sum[1];
}
//**************************Операции с матрицами****************************
double** TransMatrix(int m, int n, double** a) {
	int i, j;
	double** b;
	b = new double* [n];
	for (i = 0; i < n; i++) b[i] = new double[m];
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) b[i][j] = a[j][i];
	}
	return b;
}
//Умножение матриц
double** MultiplyMatrix(int rowsa, int colsa, int rowsb, int colsb, double** a, double** b) {
	int i, j, k;
	double t;
	double** c;
	c = new double* [rowsa];
	for (i = 0; i < rowsa; i++) c[i] = new double[colsb];

	if (colsa != rowsb) return 0;
	for (k = 0; k < colsb; k++) {
		for (i = 0; i < rowsa; i++) {
			t = 0;
			for (j = 0; j < rowsb; j++) t += a[i][j] * b[j][k];
			c[i][k] = t;
		}
	}
	return c;
}
void clearMemory(double** a, int n) { //Функция освобождения памяти, выделенной под двумерный динамический массив
	for (int i = 0; i < n; i++) {
		delete[] a[i];
	}
	delete[] a;
}
//***********************************************************************
double **InverseMatrix(double **a,int n) {
        double temp;
        int i,j,k;

        double **e;
	e = new double* [n];
	for (i = 0; i < n; i++) e[i] = new double[n];
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                e[i][j] = 0;
                 if (i == j)
                    e[i][j] = 1;
            }
 
        for (k = 0; k < n; k++) {
            temp = a[k][k];
            for (j = 0; j < n; j++) {
                a[k][j] /= temp;
                e[k][j] /= temp;
            }
 
            for (i = k + 1; i < n; i++) {
                temp = a[i][k];
                for (j = 0; j < n; j++) {
                    a[i][j] -= a[k][j] * temp;
                    e[i][j] -= e[k][j] * temp;
                }
            }
        }
       for (k = n - 1; k > 0; k--) {
            for (i = k - 1; i >= 0; i--)  {
                temp = a[i][k];
                 for (j = 0; j < n; j++)  {
                    a[i][j] -= a[k][j] * temp;
                    e[i][j] -= e[k][j] * temp;
                }
            }
        }
    return e; 
    }

void sort_1(double *x,long n) {
    double xx;
    int i, j;
    for (i = 1; i <= n; i++) {
        for (j = i + 1; j <= n; j++) {
            if (x[i] > x[j]) {
                xx = x[i]; x[i] = x[j]; x[j] = xx;
            }
        }
    }
}

//#########################################################################
void lmtexactbeta(double &beta,int n,double prob,double t) {

    double* p = new double;
    double* q = new double;
    int* which = new int;
    double* mean = new double;
    double* sd = new double;
    double* bound = new double;
    double* z = new double;
    double* phonc = new double;
    int* status = new int;
       
    double phonc1,q1,df,bound1; 
    int which1,status1;
    

    *mean = 0;
    *sd = 1;    
    *which =2;
    *p = prob;
    *q = 1 - *p;
    df=n-1.;
    which1 = 1;
    
    cdfnor(which, p, q, z, mean, sd, status, bound);
    phonc1 = *z * sqrt(n);
    cdftnc(which1,beta,q1,t,df,phonc1,status1,bound1);
    
    delete which,status,p,q,phonc,mean,sd,bound,z;

}

//#########################################################################


void lmtexact(double beta,int n,double prob,double &t) {

    double* p = new double;
    double* q = new double;
    int* which = new int;
    double* mean = new double;
    double* sd = new double;
    double* bound = new double;
    double* z = new double;
    double* phonc = new double;
    int* status = new int;

    double phonc1,q1,df,bound1; 
    int which1,status1;
    

    *mean = 0;
    *sd = 1;    
    *which = 2;
    *p = prob;
    *q = 1 - *p;
    df=n-1.;
    which1 = 2;
 
    cdfnor(which, p, q, z, mean, sd, status, bound);
    phonc1 = *z * sqrt(n);
    cdftnc(which1, beta, q1, t, df, phonc1, status1, bound1);
    
    delete which,status,p,q,phonc,mean,sd,bound,z;

}





//##################################################################
int testr(int a, int b, int c) {
    if (c < 0) return 0;
    if (a == 0 || b == 0) {
        if (c == 0) return 1;
        return 0;
    }
    return -1;
}
//*****************************************************
int signet(int a, int b) {
    if (b < 0) return 0;
    if (a == 0) {
        if (b == 0) return 1;
        return 0;
    }
    return -1;
}
//*****Точное распределение критерия Уилкоксона по рекурентной формуле*********
void wilc_Rec(int m, int n,double *p,double *wrange) {
    int mn, i, j, k;
    double s,h1,h2;
    mn = m * n+1;

    double***w;
    w = new double**[m];
    for (i = 0; i <= m; ++i) {
        w[i] = new double*[n];
        for (j = 0; j <= n; ++j)  w[i][j]=new double[mn];
    }

    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            for (k = 0; k <= mn; k++)  w[i][j][k] = 0;
        }
    }
    
    s = 0;
    for (k=0;k<mn;k++) {
        for (i=1;i<=m;i++) {
            for (j=1;j<= n;j++) {
                h1=testr(i,j-1,k-i);
                if (h1== -1) h1=(w[i][j-1][k-i]);
                h2=testr(i-1,j,k);
                if (h2== -1) h2=(w[i-1][j][k]);
                w[i][j][k]= (h1* j + h2*i)/(1.00*(i+j));
            }
        }
        s+=w[m][n][k];
        p[k]=s;
        wrange[k] = (k + m * (m + 1) / 2);
    }
}
//**Точное распределение критерия знаковых рангов Уилкоксона по рекурентной формуле(в книгах ошибка, смотри test_signe)***
void signe_Rec(int n,double *p) {
    double z;
    int i,j,k,s,h1,h2,mn;
    mn=n*(n+1)/2+1;
    int**w;
    w = new int* [n+1];
    for (i = 0;i <=n; ++i) {
        w[i] = new int[mn+1];
    }
    for (i = 0; i <=n; i++) {
        for (j = 0; j<=mn; j++) w[i][j] = 0;
    }
    s=0;
    z = pow(2, n);
    for (k=0;k<mn;k++) {
        for (i=1;i<= n;i++) {
            h1=signet(i-1,k);
            if (h1== -1) h1=w[i-1][k];
            h2=signet(i-1,k-i);
            if (h2==-1) h2=w[i-1][k-i];
            w[i][k]=h1+h2;
        }
        s+= w[n][k];
        p[k]=s/z;
    }
}

//****************Сортировка от 0*********************************************
void sort_0(double *x, long n) {
    double xx;
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            if (x[i] > x[j]) {
                xx = x[i]; x[i] = x[j]; x[j] = xx;
            }
        }
    }
}
//#######################################################################
//************cns2=n!/n1!*n2!...nkx!;mm[kx]=n1+n2+...+nkx******************
double cns2(int kx, int n, int m[]) {
    int i, j;
    double s1, s2, s3,w;
    s1 = 0;
    for (i = 1; i <= n; i++)  s1 += log(double(i));
        s3 = 0;
    for (j = 0; j < kx; j++) {
        s2 = 0;
        for (i = 1; i <= m[j]; i++)  s2 += log(double(i));
            s3+=s2;
    }
    w =exp(s1-s3);
    return w;
}
//***********cnm=n!/m!*(n-m)!********************************
double cnm(int n, int m) {
    double s1, s2;
    int i;
    s1 = 0; s2 = 0;
    for (i = m + 1; i <= n; i++) s1 += log(i);
    for (i = 1; i <= n - m; i++) s2 += log(i);
    return exp(s1 - s2);
}
//####################################################
void rndvalue(int n,double* x) {
int i;
for (i=1;i<=n;i++) x[i]=rand()/double(RAND_MAX);
}
 //#########################################################
void xmetka(double* x, int k, int* m, int* metka) {
    int i, j, n, mt;
    double xx;
    n = 0;
    for (i = 0; i < k; i++) {
        for (j = 0; j < m[i]; j++) metka[j + n] = i;
        n += m[i];
    }
    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (x[i] > x[j]) {
                xx = x[i]; x[i] = x[j]; x[j] = xx;
                mt = metka[i]; metka[i] = metka[j]; metka[j] = mt;
            }
        }
    }
}
//################################################################
void wilcoxon(double* x, int k, int n1,int n2, double alpha, double &rstat, double &zw, double &zwcrit) {
    int i, j, n, mtest, nmin;
    double ew, dw, z, t;
    int m[2];
    m[0]=n1;m[1]=n2;
    n=n1+n2;
    int* metka = new int[n];
    xmetka(x, k, m, metka);
    nmin = m[0]; mtest = 0;
    for (i = 0; i < k; i++) {
        if (m[i] < nmin) {
            nmin = m[i]; mtest = i;
        }
    }
    for (i = 0; i < k; i++) {
        rstat = 0.0;
        for (j = 0; j < n; j++)  if (metka[j] == mtest) rstat += (j + 1.0);
    }
    //rstat - Rank sum 
    ew = nmin * (n + 1.) / 2.; //E{w}
    dw = m[0] * m[1] * (n + 1.) / 12.; //D{w}
    zw = abs((rstat - ew + 0.5) / sqrt(dw)); //z(w)
    zw = zw * (1. + sqrt((n - 2) / (n - 1. - zw * zw))) / 2.; //zcor(w)
    z = invnormaldistribution(alpha / 2.);
    t = invstudenttdistr(n - 2, alpha / 2.);
    zwcrit = (-z - t) / 2.;  //w(alpha)
}
//###############################################################
void leman(double* x, int k, int n1, int n2,double alpha, double &rstat, double &zr, double &zrcrit) {
    int i,k1,k2, n;
    double er, dr,r1,r2;
    int m[2];
    m[0]=n1;m[1]=n2;
    n = m[0] + m[1];
    int* metka = new int[n];
    
    xmetka(x,k,m,metka);
    r1=0;r2=0;k1=0;k2 = 0;
    for (i=1;i<= n;i++) {
        if (metka[i-1]== 0) {
            r1+=pow(i-k1,2); k1++;
        }
        if (metka[i-1]==1) {
            r2+= pow(i-k2,2); k2++;
        }
    }
    rstat=(r1 * 1.0 / m[1] + r2 * 1.0 / m[0] + 1. / 6.) / (1.0 * m[0] * m[1]) - 2. / 3.;
    er = (1. + 1. / n) / 6.;
    dr = (1. + 1. / n) * (1. + 1. / n - 0.75 * (1. / m[0] + 1. / m[1])) / 45.;
    zr = (m[0] * m[1] * rstat / n - er) / sqrt(45. * dr) + 1. / 6.;
    zrcrit = invnormaldistribution(1. - alpha);
}
//################################################################
void series(double* x, int k, int n1,int n2, double alpha, int* col, int& ksr, int& kx, double* p, int* w) {

    int i, n, mx, nx, kk;
    double wx, s, pc;
    int m[2];
    m[0]=n1;m[1]=n2;
    n = m[0] + m[1];
    int* metka = new int[n + 1];
    int* aa = new int[n + 1];

    xmetka(x, k, m, metka);
    for (i = 1; i <= n; i++) aa[i] = metka[i - 1];
    ksr = seriesstatistic(n, aa, col);

    mx = m[0];
    if (m[1] < m[0]) mx = m[1];
    nx = n - mx; kx = 2 * mx;
    p[kx - 1] = 1; w[kx - 1] = kx + 1;
    wx = cnm(n, mx);
    s = 0;
    for (i = 1; i <= kx - 1; i++) {
        w[i - 1] = i + 1;
        if (w[i - 1] % 2 == 0) {
            kk = w[i - 1] / 2;
            pc = 2. * cnm(mx - 1, kk - 1) * cnm(nx - 1, kk - 1);
        }
        else {
            kk = (w[i - 1] - 1) / 2;
            pc = cnm(mx - 1, kk) * cnm(nx - 1, kk - 1) + cnm(mx - 1, kk - 1) * cnm(nx - 1, kk);
        }
        s += pc; p[i - 1] = s / wx;
    }
}
//###############################################################
void krusk(double* x, int k, int* m, double alpha, double& hstat, double& hcrit) {
    int i, j, n, mtest;
    double s, r, hstat1, falpha, chialpha;
    int metka[200];

    n = 0;
    for (i = 0; i < k; i++) {
        for (j = 0; j < m[i]; j++) metka[j + n] = i;
        n = n + m[i];
    }
    xmetka(x, k, m, metka);
    s = 0;
    for (i = 0; i < k; i++) {
        mtest = i; r = 0;
        for (j = 0; j < n; j++)  if (metka[j] == mtest) r += (j + 1);
        s += r * r / m[i];
    }
    hstat = 12. * s / (n * (n + 1.0)) - 3.0 * (n + 1.0);
    hstat1 = (0.5 * hstat * (1.0 + (n - k) / (n - 1.0 - hstat)));
    falpha = invfdistribution(k - 1, n - k, alpha); //F(k-1,n-k,alpha)
    chialpha = invchisquaredistribution(k - 1, alpha);  //Chi(k-1,alpha)
    hcrit = (0.5 * (falpha * (k - 1.0) + chialpha)); //Hcr 
}
//################################################
void ansari(double* x, int k, int n1, int n2, double alpha, double& astat, double& zw, double& zwcrit) {
    int mtest, i, mm, nn,n;
    double em, dm;
    int m[2];
    n = n1 + n2;
    int* metka = new int[n];
    m[0]=n1;m[1]=n2;

    xmetka(x, k, m, metka);
    mtest = 1; mm = m[1]; nn = m[0];
    if (m[0] < m[1]) {
        mm = m[0]; nn = m[1]; mtest = 0;
    }

    astat = 0;
    for (i = 0; i < n; i++) {
        if (metka[i] == mtest) astat += abs((i + 1.) - (n + 1.) / 2.0);
    }
    astat = mm * (n + 1.) / 2. - astat;

    if ((mm + nn) % 2 == 0) {
        em = mm * (n + 2.) / 4.;
        dm = mm * nn * (n - 2.) * (n + 2.) / (48. * (n - 1.0));
    }
    else {
        em = mm * (n + 1.0) * (n + 1.0) / (4. * n);
        dm = mm * nn * (n + 1.0) * (n * n + 3.) / (48. * n * n);
    }
    zw = abs(astat - em) / sqrt(dm);
    zwcrit = invnormaldistribution(1. - alpha / 2.);
}
//####################################################
void mood(double* x, int k, int n1, int n2, double alpha, double& astat, double& zw, double& zwcrit) {
    int mtest, i, mm, nn,n;
    double em, dm, corr;
    int m[2];
    n = n1 + n2;
    int* metka = new int[n];
    m[0]=n1;m[1]=n2;

    xmetka(x, k, m, metka);
    mtest = 1; mm = m[1]; nn = m[0];
    if (m[0] < m[1]) {
        mm = m[0]; nn = m[1]; mtest = 0;
    }

    astat = 0;
    for (i = 0; i < n; i++) {
        if (metka[i] == mtest) astat += pow(i + 1. - (n + 1.) / 2., 2);
    }
    em = mm * (n + 1.) * (n - 1.) / 12.;
    dm = mm * nn * (n + 1.) * (n + 2.) * (n - 2.) / 180.;
    corr = 0.5 * sqrt(45. / ((n + 1.) * (n + 2) * (n - 2))) / sqrt(dm * mm * nn);
    zw = abs((astat - em + 0.5) / sqrt(dm)) + corr;
    zwcrit = invnormaldistribution(1. - alpha / 2.);
}
//####################################################################
void david(double* x, int k, int n1, int n2, double alpha, double& astat, double& zw, double& zwcrit) {
    int mtest, i, mm, nn,n;
    double em, dm;
    int m[2];
    n = n1 + n2;
    int* metka = new int[n];
    double* rr = new double[n + 1];
   
    m[0]=n1;m[1]=n2;

    xmetka(x, k, m, metka);
    mtest = 1; mm = m[1]; nn = m[0];
    if (m[0] < m[1]) {
        mm = m[0]; nn = m[1]; mtest = 0;
    }

    for (i = 1; i <= n; i++) {
        rr[i] = (n / 2. - i + 1.);
        if (i > n / 2) rr[i] = (i - n / 2.);
    }
    astat = 0;
    for (i = 0; i < n; i++) if (metka[i] == mtest) astat += rr[i + 1];

    if ((mm + nn) % 2 == 0) {
        em = mm * (n + 2.) / 4.;
        dm = mm * nn * (n - 2.) * (n + 2.) / (48. * (n - 1.0));
    }
    else {
        em = mm * (n + 1.0) * (n + 1.0) / (4. * n);
        dm = mm * nn * (n + 1.0) * (n * n + 3.) / (48. * n * n);
    }

    zw = abs(astat - em) / sqrt(dm);
    zwcrit = invnormaldistribution(1. - alpha / 2.);
}
//####################################################################
void msigne(int n, double* x, double alpha, double median, double& bothtails, double& lefttail, double& righttail) {
    int i, gt, ne, msmal;
    gt = 0; ne = 0;
    for (i = 0; i < n; i++) {
        if (x[i] > median) gt += 1;
        if (x[i] != median) ne += 1;
    }

    if (ne == 0) {
        bothtails = 0; lefttail = 0; righttail = 0;
        return;
    }
    msmal = int(std::fmin(gt, ne - gt));
    bothtails = 2. * binom(msmal, ne, 0.5);
    lefttail = binom(gt, ne, 0.5);
    righttail = binom(ne - gt, ne, 0.5);

}
//#############################################################
void wsigne(int n, double* x, double* y, double alpha, double& wstat, double* wrange, double* wdistr) {
    int i, j, k;
    double* z = new double[n];
    double* t = new double[n];
    double zz;

    k = n * (n + 1) / 2 + 1;
    for (i = 0; i < n; i++) {
        z[i] = abs(x[i] - y[i]);
        t[i] = x[i] - y[i];
    }

    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (z[i] > z[j]) {
                zz = z[i]; z[i] = z[j]; z[j] = zz;
                zz = t[i]; t[i] = t[j]; t[j] = zz;
            }
        }
    }
    wstat = 0;
    for (i = 0; i < n; i++)if (t[i] > 0) wstat += i + 1;
    for(i = 0; i < k; i++) wrange[i] = i;
signe_Rec(n, wdistr);
}
//#######################################################
void sigel(double* x, int k, int n1, int n2, double alpha, double* wrange, double* wdistr, int& rmin, int& rmax, double& zstat, double& zcrit) {

    int i, j, jj, kk, kx, mmin, mmax, mx, nx,n,m[2];
    n=n1+n2;
    m[0]=n1;m[1]=n2;
    int* metka = new int[n];
    int* a = new int[n + 1];
    double sr,mr;

    xmetka(x, k, m, metka);
//*****************Transformation  Siegel**************************
    kk = 1; j = 0; jj = 0; a[0] = 1;
    for (i = 2; i <= n; i++) {
        if (j == 0 || j == 1) {
            kx = n - jj; jj++; a[kx - 1] = i;
        }
        if (j == 2 || j == 3) {
            kk++;a[kk - 1] = i;
        }
        j++;
        if (j > 3) j = 0;
    }
//**********************************************************
    rmin = 0; mmin = 1; mmax = 0; rmax = 0;
    if (m[0] <= m[1]) {
        mmin = 0; mmax = 1;
    }
    for (i = 0; i < n; i++) {
        if (metka[i] == mmin) rmin += a[i];
        if (metka[i] == mmax) rmax += a[i];
    }

    //rmax=n*(n+1)/2-rmin; rmin+rmax=n*(n+1)/2;

    mx = m[1]; nx = m[0];
    if (m[0] < m[1]) {
        mx = m[0]; nx = m[1];
    }

    wilc_Rec(mx, nx, wdistr, wrange);

    sr = sqrt(mx * nx * (n + 1.) / 12.);
    mr = abs(rmin - 0.5 * mx * (n + 1));
    zstat = (mr - 0.5) / sr;
    zcrit = invnormaldistribution(1. - alpha / 2.);

    /* Уточнение Iman & Davenport
    double z,t;
    zstat=zstat*(1. + sqrt((n - 2) / (n - 1. - zstat*zstat))) / 2.;
    z = invnormaldistribution(alpha / 2.);
    t = invstudenttdistr(n - 2, alpha / 2.);
    zcrit = (-z - t) / 2.;
    */
}


//*********************Ansari**start1********************************************
int start1(int n, int* f) {
    int i, lout;

    lout = int(1 + n / 2);
    for (i = 1; i <= lout; i++) f[i] = 2;
    if ((n % 2) == 0) f[lout] = 1;
    return(lout);
}
//**********************Ansari start2**************************************
int start2(int n, int* f) {
    int one, two, three, four, i, j, a, b, lt1, ndo, nu, lout;

    one = 1; two = 2; three = 3; four = 4;
    nu = n - n % 2;
    j = nu + 1; lout = j; lt1 = lout + 1;
    ndo = int(lt1 / 2);
    a = one; b = three;
    for (i = 1; i <= ndo; i++) {
        f[i] = a; f[j] = a; j = j - 1; a = a + b; b = four - b;
    }
    if (nu == n) return(lout);
    nu = ndo + 1;
    for (i = nu; i <= lout; i++) f[i] = f[i] + two;
    f[lt1] = two; lout = lt1;
    return(lout);
}
//***********Точное распределение критерия Ансари-Брэдли*(аналитика)*************
long int ansari_exact(int *m,int* astat, double* pw) {
    int min_val, max_val, i,nrows, a0;
    double sum, s;
    int* w;

    min_val = m[0];
    max_val = m[1];
    nrows = 1 + m[0]*m[1];
    w = new int[nrows];
    
    nrows = gscale(min_val, max_val, w);
    a0 = int((min_val + 1) / 2) * (1 + int(min_val / 2));
    
     sum = 0;
    for (i = 1; i <= nrows; i++) {
        astat[i - 1] = a0 + i - 1;
       sum += w[i];
    }
    
     s = 0;
    for (i = 0; i < nrows; i++) {
        s = s + w[i+1];
        pw[i] =  s / sum;
    }
    
    delete[] w;
    return long int(nrows);
}
//************ Критерий серий точное распределение**(аналитика)************
int series_exact(int* m, double *wrange,double *pw) {
    int nn,k, kcur, m1, m2, n, i, ksr;
    double pc, w2, s;
    double *w;

    nn = m[0] * m[1] + 1;
    w = new double[nn];

    k = 2;
    m1 = m[0]; m2 = m[1]; n = 0;
    for (i = 0; i < k; i++)  n = n + m[i];
    w2 = cnm(n, m1);
    s = 0;
    for (i = 1; i <= 5000; i++) {
        kcur = i + 1;
        if (kcur % 2 == 0) {
            pc = double(2. * cnm(m1 - 1, kcur / 2 - 1) * cnm(m2 - 1, kcur / 2 - 1));
        }
        else {
            pc = double(cnm(m1 - 1, (kcur - 1) / 2) * cnm(m2 - 1, (kcur - 1) / 2 - 1) + cnm(m1 - 1, (kcur - 1) / 2 - 1) * cnm(m2 - 1, (kcur - 1) / 2));
        }
        s += pc;
        w[i-1]=pc;
        pw[i-1]=(s / w2);
        wrange[i-1]=double(kcur);
        if (pw[i - 1] > 1) break;
    }
    ksr = i - 1;
    return ksr;
}
//********Критерий Уилкоксона (точное распределение)*(аналитика)***********************

long int wilcoxon_exact(int *m, double *wrange, double *pw) {

    int* work, minmn, maxmn, inx;
    int* w;
    int n1, kk, k, j, i, sum;
    long int mn;

    mn = m[0]*m[1]+1; 
    maxmn = fmax(m[0],m[1]);
    minmn = fmin(m[0],m[1]);
    n1 = maxmn + 1;
    work = new int[mn + 2];
    w = new int[mn + 2];
   
    for (i = 1; i <= n1; i++)  w[i] = 1;
    n1++;
    for (i = n1; i <= mn; i++) w[i] = 0;
    work[1] = 0; inx = maxmn;

    for (i = 2; i <= minmn; i++) {
        work[i] = 0; inx = inx + maxmn; n1 = inx + 2; kk = 1 + inx / 2; k = i;
        for (j = 1; j <= kk; j++) {
            k++;  n1--; sum = w[j] + work[j];
            w[j] = sum; work[k] = sum - w[n1]; w[n1] = sum;
        }
    }

    for (i = 0; i < mn; i++)  pw[i]=double(w[i + 1]);
    sum = 0;
    for (i = 0; i < mn; i++) {
        sum += pw[i];
        wrange[i]=double(i + minmn * (minmn + 1.) / 2.); //Wilcoxon
        //wrange.push_back((double(i)); //Mann-Whitney
        pw[i] = sum;
    }
    for (i = 0; i < mn; i++) pw[i] = pw[i] / sum;
    delete[] work;
    return mn;
}

int seriesstatistic(int n, int* aa, int* col) {
    int ks, i, ksr;

    ks = 1; ksr = 1;col[ksr] =1;
    for (i=1;i<n;i++) {
        if (aa[i]==aa[i+1]) {
            ks++;col[ksr] = ks;
        }
        else {
            ks = 1; ksr++;col[ksr] = ks;
        }
    }
    return ksr;
}


//########################################################################
void standart(int n,double *x,double &mean,double &s){
	 	int i;
	 	mean=0;s=0;
		for(i=0;i<n;i++){
			mean+=x[i];
			s+=x[i]*x[i];
		}
		mean/=n;
		s=(s-pow(mean,2)*n)/(n-1);
		s=pow(s,0.5);
	}
//*************************************************************************
long int TestPerm(int kk,int* m,vector<double>&h,double(*critfun)(int*, int, int, int*)) {
	
        int* a;
        int n,k, j, l,i, r,km;
        double z;
        long int num;
        num=0;

        //kk=sizeof(m)/sizeof(m[0]);
        n = 0;
        for (i = 0; i < kk; i++)  n += m[i];
        a = new int[n + 1];
        km = 0;
        for (i = 0; i < kk; i++) {
            for (j = 0; j < m[i]; j++) a[j + km] = i + 1;
            km = km + m[i];
        }


        z =critfun(a, kk, n, m);
        h.push_back(z);

        while(1 > 0) {
            j = n - 2;
        while (j >= 0 && a[j] >= a[j + 1]) j--;
        if (j < 0) return num;
        k = n - 1;
        while (a[j] >= a[k]) k--;
        swap1(a[j], a[k]);
        l = j + 1; r = n - 1;
        while (l < r)  swap1(a[l++], a[r--]);
        num++;
        z =critfun(a, kk, n, m);
        h.push_back(z);
      }
        return num;
    }
//#######################################################
double ansaristatistic(int* a,int kx,int n, int* m) {  //Ansari-Breadly
    int i, j, msmal;
    double r;
    msmal = m[0];
    r = 0;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= msmal; j++)    if (a[i] == j) r += ((n*1.0 + 1.0)/2. - abs(i*1.0 - (n*1.0 + 1.0)/2.0));
    }
    return r;
}
//*********************************************
double wilcoxonstatistic(int* a,int kx, int n,int* m) {
    int msmal, i, j;
    double r;
    msmal = m[0];
    r = 0;
    for (i=1;i<= n;i++) {
        for (j=1;j<=msmal;j++) if (a[i]==j) r += i;
    }
    return r;
}
//################################################################
void qsortRecursive(double *mas, long int size) {
    long int i = 0;
    long int j = size - 1;
    double tmp,mid;
    mid = mas[size / 2];

    do {
        while (mas[i] < mid)  i++;
        while (mas[j] > mid)  j--;
        if (i <= j) {
            tmp = mas[i];mas[i] = mas[j];mas[j] = tmp; i++; j--;
        }
    } while (i <= j);
    if (j > 0)  qsortRecursive(mas, j + 1);
    if (i < size) qsortRecursive(&mas[i], size - i);
}
//*************************************************************
double fisherstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r,vrs;
      r=0;
       for (i=0;i<n;i++) { 
         ordern(n,(i+1.)/(n+1.),(i+1.)/(n+1.),er,vrs);
         if (a[i]==1) r+=er;
        }
         return r; 
     }
//**************************************************************
double vandervardenstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r;
      r=0;
       for (i=0;i<n;i++) { 
         er=invnormaldistribution((i+1.)/(m[0]+m[1]+1.));
         if (a[i]==1) r+=er;
        }
         return r; 
     }
//**************************************************************
double caponstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r,vrs;
      r=0;
       for (i=0;i<n;i++) { 
        ordern(n,(i+1.)/(n+1.),(i+1.)/(n+1.),er,vrs);
        er=vrs+er*er;
        if (a[i]==1) r+=er;
       }
         return r; 
     }
//**************************************************************
double klotzstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r;
      r=0;
       for (i=0;i<n;i++) { 
          er=invnormaldistribution((i+1.)/(m[0]+m[1]+1.));
          if (a[i]==1) r+=er*er;
        }
        
   return r;
  }
//********************************************************
double moodstatistic(int* a, int kx,int n,int* m) {
    int i, msmal;
    double r;
    msmal = 2;
    if (m[0] < m[1]) msmal = 1;
    r= 0;
    for (i = 0; i <n; i++)  if (a[i]==msmal) r+=pow((i+1.- (n + 1.) / 2.),2);
    return r;
}
//********************************************************
double davidstatistic(int* a, int kx,int n,int* m) {
    int i, msmal;
    double r;
    msmal = 2;
    if (m[0] < m[1]) msmal = 1;

    r = 0;
    for (i=0;i<n;i++) {
         if (a[i]==msmal) {
          if(i<=n/2) r+=(n/2.-i);
          if(i>n/2) r+=(i-n/2.-1);
        }
    }
    return r;
}
//******************************************************************
double lemanstatistic(int* a, int kx, int n, int* m) {
    int i, k1,k2;
    double zleman,r1,r2;
    r1 = 0.; r2 = 0.; k1 = 0; k2 = 0;
    for (i = 0; i < n; i++) {
        if (a[i] == 1) {
            r1=r1+ pow(i - k1,2); k1++;
        }
        if (a[i] == 2) {
            r2 += pow(i - k2,2); k2++;
        }
    }
    zleman = (r1 * m[0] + r2 * m[1] + m[0] * m[1] / 6.) / (m[0] * m[0] * m[1] * m[1]) - 2. / 3.;
    return zleman;
}
//**********************************************************************
double kruskalstatistic(int* aa, int kx, int n, int* m) {
    int i,j;
    double s,r;
    s=0.;
     for (j=0;j<kx;j++) {
        r=0.;
        for (i=0;i<n;i++) if (aa[i]==j+1) r+=i+1;
        s=s+r*r/m[j];
    }
    return 12.*s/(1.0*n*(n+1.))-3.*(1.0*n+1.0);
}
//***************Ansari frqadd***********************************************
int* frqadd(int* f1, int* f2, int l1in, int l1out, int l2, int nstart) {

    int i1, i2, nxt, *fadd;
    
     fadd = new int[3];
    i2 = 1;
    for (i1 = nstart; i1 <= l1in; i1++) {
        f1[i1] = f1[i1] + 2 * f2[i2];
        i2++;
    }
    nxt = l1in + 1;
    l1out = l2 + nstart - 1;
    for (i1 = nxt; i1 <= l1out; i1++) {
        f1[i1] = 2 * f2[i2];
        i2++;
    }
    nstart++;
    fadd[1] = l1out; fadd[2] = nstart;
    return(fadd);
}
//***************Ansari imply*************************************************
int* imply(int* f1, int* f2, int l1in, int l1out, int l2, int noff) {
    int sum, diff, i2, i1, j2, j1, j2min, ndo;
	int *fimply;

     fimply = new int[4];
     
    i2 = 1 - noff; j1 = l1out; j2 = l1out - noff; l2 = j2;
    j2min = int((j2 + 1) / 2);
    ndo = int((l1out + 1) / 2);

    for (i1 = 1; i1 <= ndo; i1++) {
        if (i2 > 0) {
            sum = f1[i1] + f2[i2];
            f1[i1] = sum;
        }
        else {
            sum = f1[i1];
        }
        i2 = i2 + 1;
        if (j2 >= j2min) {
            if (j1 <= l1in) {
                diff = sum - f1[j1];
            }
            else {
                diff = sum;
            }

            f2[i1] = diff; f2[j2] = diff; j2 = j2 - 1;
        }
        f1[j1] = sum; j1 = j1 - 1;
    }
    fimply[1] = l1out; fimply[2] = l2; fimply[3] = noff;
    return(fimply);
}
//**************Ansari gscale*******************************************
int gscale(int test, int other, int* pw) {

    int  i, m, lres, mm1, nm1, nm2, ier, mnow, ks, j, ndo;
    int n, ln1, ln2, nc, l1out, l2out, n2b1, n2b2, ln3, kk, ai, z;
    int* a2, * a3, * fadd, * fimply;
    bool symm;

    ln1 = 0; ln2 = 0; l1out = 0; l2out = 0; ln1 = 0; ln2 = 0; ln3 = 0; n2b1 = 0; n2b2 = 0; kk = 0; nc = 0;
    ndo = 0; ier = 0; ks = 0; j = 0; i = 0;

    m = int(fmin(test, other));
    if (m < 0) return(0);
    n = int(fmax(test, other));
    lres = 1 + int((m * n) / 2);
    fadd = new int[lres];
    fimply = new int[lres];
    a2 = new int[lres];
    a3 = new int[lres];

    symm = false;
    z = (m + n) % 2;
    if (z == 0) symm = true;
    mm1 = m - 1;
    //*****************************************************         
    if (m <= 2) {
        if (mm1 < 0) {
            pw[1] = 1; return(lres);
        }

        if (mm1 == 0) ln1 = start1(n, pw);
        if (mm1 > 0)  ln1 = start2(n, pw);
        if (symm || (other > test)) return(lres);
        j = lres;
        ndo = int(lres / 2);
        for (i = 1; i <= ndo; i++) {
            ai = pw[i];
            pw[i] = pw[j];
            pw[j] = ai;
            j = j - 1;
        }
        return(lres);
    }
    //***********************************************************
    nm1 = n - 1; nm2 = n - 2; mnow = 3; nc = 3; ier = 0;
    //**************************************************************
    while (true) {
        if (ier == 0) {
            if ((n % 2) != 1) {
                n2b1 = 3;
                n2b2 = 2;
                ln1 = start2(n, pw);
                ln3 = start2(nm2, a3);
                ln2 = start1(nm1, a2);
                //***********************************************************************          
                fadd = frqadd(a2, a3, ln2, l2out, ln3, n2b2);
                l2out = fadd[1]; n2b2 = fadd[2];
                ln2 = ln2 + nm1;
                fimply = imply(a2, a3, l2out, ln2, j, nc);
                ln2 = fimply[1]; j = fimply[2]; nc = fimply[3];
                //***************************************************************************
                nc = nc + 1;
                if (mnow == m) break;
                mnow = mnow + 1;
            }
            else {
                n2b1 = 2;
                n2b2 = 3;
                ln1 = start1(n, pw);
                ln2 = start2(nm1, a2);
            }
        }
        //*******************************************************************
        fadd = frqadd(pw, a2, ln1, l1out, ln2, n2b1);
        l1out = fadd[1]; n2b1 = fadd[2];
        ln1 = ln1 + n;
        fimply = imply(pw, a3, l1out, ln1, ln3, nc);
        ln1 = fimply[1]; ln3 = fimply[2]; nc = fimply[3];
        nc = nc + 1;
        if (mnow == m) break;
        mnow = mnow + 1;
        fadd = frqadd(a2, a3, ln2, l2out, ln3, n2b2);
        l2out = fadd[1]; n2b2 = fadd[2];
        ln2 = ln2 + nm1;
        fimply = imply(a2, a3, l2out, ln2, j, nc);
        ln2 = fimply[1]; j = fimply[2]; nc = fimply[3];
        //***************************************************************************
        nc = nc + 1;
        if (mnow == m) break;
        mnow = mnow + 1;
        ier = 1;
    }
    //********************************************************************
    if (symm) return(lres);
    ks = int((m + 3) / 2);
    j = 1;
    for (i = ks; i <= lres; i++) {
        if (i > ln1) {
            pw[i] = a2[j];
        }
        else {
            pw[i] = pw[i] + a2[j];
        }
        j = j + 1;
    }
    if (other < test) return(lres);
    j = lres;
    ndo = int(lres / 2);
    for (i = 1; i <= ndo; i++) {
        ai = pw[i];
        pw[i] = pw[j];
        pw[j] = ai;
        j = j - 1;
    }
    return(lres);
}
//************************************************************
void swap1(int& ax, int& ay) {
        int s = ax; ax = ay; ay = s;
    }
void swap(double& ax, double& ay) {
            double s = ax; ax = ay; ay = s;
        }
//********************************************************
long int count_perm(int k, int* m) {
        int i, n;
        long int knum;
        double z, w;
        w = 1.; n = 0;
        for (i = 0; i < k; i++) {
            z = w * tgamma((m[i]) + 1.); w = z; n = n + m[i];
        }
        w = tgamma(n + 1.) / w;
        knum = long(w);
        return knum;
    }
//############################################################


/*
 *  Purpose:
 *
 *    NELMIN minimizes a function using the Nelder-Mead algorithm.
 *
 *  Discussion:
 *
 *    This routine seeks the minimum value of a user-specified function.
 *
 *    Simplex function minimisation procedure due to Nelder+Mead(1965),
 *    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
 *    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
 *    25, 97) and Hill(1978, 27, 380-2).
 *
 *    The function to be minimized must be defined by a function(al) of
 *    the form
 *
 *      real fn ( const std::array<real,n>& x )
 *
 *    where "real" can be any floating-point type, e.g. double.
 *
 *    This routine does not include a termination test using the
 *    fitting of a quadratic surface.
 *
 *  Licensing:
 *
 *    This code is distributed under the GNU LGPL license.
 *
 *  Author:
 *
 *    Original FORTRAN77 version by R ONeill.
 *    C version by John Burkardt (last modified 28 October 2010).
 *    Port to modern C++ by Piotr Różański. Numerical behaviour unchanged; changes restricted to secondary issues:
 *    a) output variables moved to returned struct type
 *    b) function is now passed as std::function which allows using objects as well as traditional functions
 *    c) floating-point type and number of variables are now template arguments
 *    d) std::array is now used instead of pointers to raw arrays
 *    e) overall comments and code formatting
 *
 *  Reference:
 *
 *    John Nelder, Roger Mead,
 *    A simplex method for function minimization,
 *    Computer Journal,
 *    Volume 7, 1965, pages 308-313.
 *
 *    R ONeill,
 *    Algorithm AS 47:
 *    Function Minimization Using a Simplex Procedure,
 *    Applied Statistics,
 *    Volume 20, Number 3, 1971, pages 338-345.
 *
 *    This C++ version is available at
 *    https://github.com/develancer/nelder-mead/
 */
//#ifndef PTR_NELDER_MEAD_H
//#define PTR_NELDER_MEAD_H

#include <array>
#include <climits>
#include <functional>

/**
 * Plain data object with output information from the run of nelder_mead routine.
 *
 * @tparam real floating-point type to be used, e.g. double
 * @tparam n the number of variables
 */
template<typename real, int n>
struct nelder_mead_result {
    std::array<real,n> xmin;
    real ynewlo;
    int icount;
    int numres;
    int ifault;
};

/**
 * This routine seeks the minimum value of a user-specified function.
 *
 * @tparam real floating-point type to be used, e.g. double
 * @tparam n the number of variables
 * @param fn the function to be minimized
 * @param start a starting point for the iteration
 * @param reqmin the terminating limit for the variance of function values
 * @param step determines the size and shape of the initial simplex;
 * the relative magnitudes of its elements should reflect the units of the variables
 * @param konvge the convergence check is carried out every konvge iterations
 * @param kcount the maximum number of function evaluations
 * @return structure with output information
 */
template<typename real, int n>
nelder_mead_result<real,n> nelder_mead(
        const std::function<real(const std::array<real,n> &)> &fn,
        std::array<real,n> start,
        real reqmin,
        const std::array<real,n> &step,
        int konvge = 1,
        int kcount = INT_MAX
) {
    const real ccoeff = 0.5;
    real del;
    const real ecoeff = 2.0;
    const real eps = 0.001;
    int ihi;
    int ilo;
    int jcount;
    int l;
    const real rcoeff = 1.0;
    real rq;
    real x;
    real y2star;
    real ylo;
    real ystar;
    real z;

    nelder_mead_result<real,n> result;

    // Check the input parameters.
    if (reqmin <= 0.0 || n < 1 || konvge < 1) {
        result.ifault = 1;
        return result;
    }

    std::array<real,n> p[n + 1];
    std::array<real,n> pstar, p2star, pbar;
    real y[n + 1];

    result.icount = 0;
    result.numres = 0;

    jcount = konvge;
    del = 1.0;
    rq = reqmin * n;

    // Initial or restarted loop.
    for (;;) {
        p[n] = start;
        y[n] = fn(start);
        result.icount++;

        for (int j = 0; j < n; j++) {
            x = start[j];
            start[j] += step[j] * del;
            p[j] = start;
            y[j] = fn(start);
            result.icount++;
            start[j] = x;
        }
        // The simplex construction is complete.

        // Find highest and lowest Y values.
        // YNEWLO = Y(IHI) indicates the vertex of the simplex to be replaced.
        ylo = y[0];
        ilo = 0;

        for (int i = 1; i < n + 1; i++) {
            if (y[i] < ylo) {
                ylo = y[i];
                ilo = i;
            }
        }
        // Inner loop.
        for (;;) {
            if (kcount <= result.icount) {
                break;
            }
            result.ynewlo = y[0];
            ihi = 0;

            for (int i = 1; i < n + 1; i++) {
                if (result.ynewlo < y[i]) {
                    result.ynewlo = y[i];
                    ihi = i;
                }
            }
            // Calculate PBAR, the centroid of the simplex vertices
            // excepting the vertex with Y value YNEWLO.
            for (int i = 0; i < n; i++) {
                z = 0.0;
                for (int j = 0; j < n + 1; j++) {
                    z += p[j][i];
                }
                z -= p[ihi][i];
                pbar[i] = z / n;
            }
            // Reflection through the centroid.
            for (int i = 0; i < n; i++) {
                pstar[i] = pbar[i] + rcoeff * (pbar[i] - p[ihi][i]);
            }
            ystar = fn(pstar);
            result.icount++;

            // Successful reflection, so extension.
            if (ystar < ylo) {
                for (int i = 0; i < n; i++) {
                    p2star[i] = pbar[i] + ecoeff * (pstar[i] - pbar[i]);
                }
                y2star = fn(p2star);
                result.icount++;

                // Check extension.
                if (ystar < y2star) {
                    p[ihi] = pstar;
                    y[ihi] = ystar;
                } else {
                    // Retain extension or contraction.
                    p[ihi] = p2star;
                    y[ihi] = y2star;
                }
            } else {
                // No extension.
                l = 0;
                for (int i = 0; i < n + 1; i++) {
                    if (ystar < y[i]) {
                        l++;
                    }
                }

                if (1 < l) {
                    p[ihi] = pstar;
                    y[ihi] = ystar;
                } else if (l == 0) {
                    // Contraction on the Y(IHI) side of the centroid.
                    for (int i = 0; i < n; i++) {
                        p2star[i] = pbar[i] + ccoeff * (p[ihi][i] - pbar[i]);
                    }
                    y2star = fn(p2star);
                    result.icount++;
                    // Contract the whole simplex.
                    if (y[ihi] < y2star) {
                        for (int j = 0; j < n + 1; j++) {
                            for (int i = 0; i < n; i++) {
                                p[j][i] = (p[j][i] + p[ilo][i]) * 0.5;
                            }
                            result.xmin = p[j];
                            y[j] = fn(result.xmin);
                            result.icount++;
                        }
                        ylo = y[0];
                        ilo = 0;

                        for (int i = 1; i < n + 1; i++) {
                            if (y[i] < ylo) {
                                ylo = y[i];
                                ilo = i;
                            }
                        }
                        continue;
                    } else {
                        // Retain contraction.
                        p[ihi] = p2star;
                        y[ihi] = y2star;
                    }
                } else if (l == 1) {
                    // Contraction on the reflection side of the centroid.
                    for (int i = 0; i < n; i++) {
                        p2star[i] = pbar[i] + ccoeff * (pstar[i] - pbar[i]);
                    }
                    y2star = fn(p2star);
                    result.icount++;
                    // Retain reflection?
                    if (y2star <= ystar) {
                        p[ihi] = p2star;
                        y[ihi] = y2star;
                    } else {
                        p[ihi] = pstar;
                        y[ihi] = ystar;
                    }
                }
            }
            // Check if YLO improved.
            if (y[ihi] < ylo) {
                ylo = y[ihi];
                ilo = ihi;
            }
            jcount--;

            if (0 < jcount) {
                continue;
            }
            // Check to see if minimum reached.
            if (result.icount <= kcount) {
                jcount = konvge;

                z = 0.0;
                for (int i = 0; i < n + 1; i++) {
                    z += y[i];
                }
                x = z / (n + 1);

                z = 0.0;
                for (int i = 0; i < n + 1; i++) {
                    real yx = y[i] - x;
                    z += yx * yx;
                }

                if (z <= rq) {
                    break;
                }
            }
        }
        // Factorial tests to check that YNEWLO is a local minimum.
        result.xmin = p[ilo];
        result.ynewlo = y[ilo];

        if (kcount < result.icount) {
            result.ifault = 2;
            break;
        }

        result.ifault = 0;

        for (int i = 0; i < n; i++) {
            del = step[i] * eps;
            result.xmin[i] += del;
            z = fn(result.xmin);
            result.icount++;
            if (z < result.ynewlo) {
                result.ifault = 2;
                break;
            }
            result.xmin[i] -= del + del;
            z = fn(result.xmin);
            result.icount++;
            if (z < result.ynewlo) {
                result.ifault = 2;
                break;
            }
            result.xmin[i] += del;
        }

        if (result.ifault == 0) {
            break;
        }
        // Restart the procedure.
        start = result.xmin;
        del = eps;
        result.numres++;
    }

    return result;
}

//#endif // PTR_NELDER_MEAD_H


 double Normal_function_to_minimize(const array<double, 2>& xsimpl) {
        double s1, s2, s3, s4, z, psi, p, d, c1, c2;
        int i, kx;
        s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;
        if (xsimpl[1] <= 0) return 10000;
        for (i = 0; i < sm.n; i++) {
            z = (sm.x[i] - xsimpl[0]) / xsimpl[1];
            d = spi * exp(-z * z / 2.);
            p = normaldistribution(z);
            psi = d / (1. - p);
            s1 += (1. - sm.r[i]) * (sm.x[i] - xsimpl[0]);
            s2 += (1. - sm.r[i]) * pow(sm.x[i] - xsimpl[0], 2);
            s3 += sm.r[i] * psi;
            s4 += sm.r[i] * psi * z;
            kx += 1 - sm.r[i];
        }
        c1 = s1 + xsimpl[1] * s3;
        c2 = s2 + pow(xsimpl[1],2) * (s4 - kx);
        z = c1 * c1 + c2 * c2;
        return(z);
    }

//#################MLE Weibull Minimized Censorized Function#######################
double Weibull_function_to_minimize(const array<double,1>& xsimpl) {
    
 double s1,s2,s3,z,b,c;
 int i,k,n;
 n = sm.n;
   s1=0;s2=0;s3=0;k=0;
   b=xsimpl[0];
 for(i=0;i<n;i++) {
   k+=(1-sm.r[i]);
   s1+=pow(sm.x[i],b);
 }
   c=s1/k;

for(i=0;i<n;i++) {
  z=(pow(sm.x[i],b))/c;
  s3+=z*log(z);
  s2+=(1-sm.r[i])*log(z);
}
 c=s3-s2-k;
 return c*c;
}
//##############MLE Up-down Minimize Function Normal,Weibull#################################################################
//
 double UpDownMinimize(const array<double,2>& xsimpl) {

  double z,p,d,s1,s2,fiz;
  static double pi = 3.1415926535898e0;
  int i;

  if(xsimpl[1]<=0) return(10000000.0);

  s1=0;s2=0;
  for(i=0;i<smf.k;i++) {
    z=(smf.x[i]-xsimpl[0])/xsimpl[1];
   if(smf.ts=="Normal") {
      p=normaldistribution(z);
      d=exp(-z*z/2)/sqrt(2.*pi);
   }
    if(smf.ts=="Weibull") {
     p = 1.-exp(-exp(z));
     d=exp(z-(exp(z)));
    }
  if(p<=0 || p>=1) return(10000000.0);
   fiz=(smf.p[i]-p)*smf.n[i]*d/(p*(1.-p));
   s1+=fiz;
   s2+=fiz*z;
 }
 return(s1*s1+s2*s2);
}
//####################################################################
 void CovMatrixUpDown(string ts,double stepx,double bint,double **&v) {
	double *zcov,*pcov,*qcov,*dcov,*wcov,*f1;
	double z2, s2, s1;
	int kv, m, m1, kk,i1, j1,i,j;
	static double pi = 3.1415926535898e0;
	if (ts == "Normal") m =int(4/stepx)+1;
	if (ts == "Weibull") m =int(3/stepx)+1;
	m1 = 2 * m;s2 = 0;kk = 0;kv = 2;z2 = 0;s1 = 1;
	zcov = new double[m1];
	pcov = new double[m1];
	qcov = new double[m1];
	dcov = new double[m1];
	wcov = new double[m1];
	for (i = 0; i < m; i++) {
 		zcov[i] = bint - stepx * (m - i);
		zcov[i+m] = bint + stepx * (i);
	}
	for (i = 0; i < m1; i++) {
		if (ts == "Normal") {
		pcov[i] = normaldistribution(zcov[i]);
		dcov[i] = exp(-zcov[i] * zcov[i] / 2) / sqrt(2. * pi);
		qcov[i] = 1. - pcov[i];
		}
		if (ts == "Weibull") {
		pcov[i] = 1. - exp(-exp(zcov[i]));
		dcov[i] = exp(zcov[i] - (exp(zcov[i])));
		qcov[i] = 1. - pcov[i];
		}
	}
	    j=0;
//#########################################
	for(i=0;i<m1;i++) {
	  s1=1.;
	  if(kk==1) {
		for(i1=j;i1<i;i1++) s1*=qcov[i1]/pcov[i1];
		  }
	 else {
	//###########################################
	 if(pcov[i]<qcov[i]) {
		 for(i1=i;i1<m1;i1++) {
			if(pcov[i1]>=qcov[i1]) break;
			s1*=pcov[i1]/qcov[i1];
		  }
	}
	   else {
	kk=1;j=i;
	  }
	//###########################################
	  }
	  wcov[i]=s1;s2+=wcov[i];
	}
//###########################################
	f1 = new double[kv];
	for(i=0;i<m1;i++) {
	 f1[0]=-dcov[i];
	 f1[1] = -dcov[i] * zcov[i];
	 z2=wcov[i]/(pcov[i]*pcov[i]*qcov[i]);
	   for(i1=0;i1<kv;i1++) {
		 for(j1=0;j1<kv;j1++) {
		   if(wcov[i]==0 || pcov[i]==0 || pcov[i]==1) {
			 s1=0;
			} else {
			 s1=f1[i1]*f1[j1]*z2;
			  }
			  v[i1][j1]=v[i1][j1]+s1;
		   }
		  }
	}
	for(i=0;i<kv;i++) {
		  for(j=0;j<kv;j++) v[i][j]=v[i][j]/(2.*s2);
	}
	  v=InverseMatrix(v,kv);
	  delete[] zcov, pcov, qcov,wcov,f1;
      }
//##################################################################
void MleastSquare(int n,int k,double **x,double **y,double **&db,double **&b,double *&yr) {

   int i,j;
  double s;
 
  db = InverseMatrix((MultiplyMatrix(k, n, n, k, TransMatrix(n, k, x), x)), k);  //covariance matrix factors (k x k)
  b = MultiplyMatrix(k, n, n, 1, MultiplyMatrix(k, k, k, n, db, TransMatrix(n, k, x)), y); // coef
  for (i = 0; i < n; i++) {
   s = 0;
    for (j = 0; j < k; j++)  s += b[j][0] * x[i][j];
      yr[i] = s;
   }

}
//##################################################################
void MleastSquare_weight(int n,int k,double **x,double **y,double **v,double **&db,double **&b,double *&yr) {

   int i,j;
  double s,**xt;

   xt=new double *[n];
   for (i = 0; i < n; i++)  xt[i] = new double[n];
   v=InverseMatrix(v,n); // (n x n)
   xt=MultiplyMatrix(k,n,n,n,TransMatrix(n, k, x), v); // (k x n)

  db = InverseMatrix(MultiplyMatrix(k, n, n, k, xt, x),k);  //covariance matrix factors (k x k)
  b = MultiplyMatrix(k, n, n, 1, MultiplyMatrix(k, k, k, n, db, xt), y); // coef b[k][0]
  for (i = 0; i < n; i++) {
   s = 0;
    for (j = 0; j < k; j++)  s += b[j][0] * x[i][j];
      yr[i] = s;
   }
   delete [] xt;
}


//###########################################################################

//*********************************************************************
 double sf35r(double x) {
  double t,s,s1,s2;   
  t = abs(x);
  s=0;
   if (t<6.5) {
      s1 =exp(-t*t)*((((((0.56419*t+6.802899)*t+38.71143)*t+131.1266)*t+278.5978)*t+355.969)*t+224.1828);
      s2=(((((((t+12.05784)*t+69.11384)*t+238.4503)*t+527.5538)*t+741.5214)*t+608.9322)*t+224.1828);
      s=s1/s2;
   }
      if (x>=0) return s;
      if (x<0) return (2-s);
 }
//***********************************************************************
double sf49r( double x) {
  double z;
  z=0.5*sf35r(-0.7071067*x);
  return z;
}
//*********************************************************************
double sf53r(double y,double z,double eps) {
   double sys076,sys017,expov,c,ep1,t,b,a,ta,hsqb,bexp,asq,a4,b4,a4b4,ahsqb,ab4,f,sum,g,g1,ber,ter,d1;
   double d2,d,aeps;

   sys076 = 88.72283;
   sys017 = 0.1591549;
   expov = 88.72283;
   c = 0.1591549;
   ep1 = eps;
   if (eps==0) ep1=0.000001;
   t=0;
   b=abs(y);
   a=abs(z);
   if (a==0) return t;
   ta =atan(a);
   if (a*b>4) {
     t = sf49r(b);
     t =c*(ta+atan(1./a))-0.5*(t-0.5);
     if (z<0) t=-t;
     return t;
   }
   hsqb=0.5*b*b;
  if (hsqb>expov) return t;
   bexp=exp(-hsqb);
   asq=a*a;
   a4=asq*asq;
   b4=hsqb*hsqb;
   a4b4=a4*b4;
   ahsqb=a*hsqb;
   ab4=a*b4*0.5;
    f=1;
    sum=0;
    g=3;

  while(1<0) {
  g1=g;
  ber=0;
  ter=ab4;

  while (1<0) {
 ber=ber+ter;
 if (ter<=ber*ep1) break;
  ter=ter*hsqb/g1;
  g1=g1+1;
 }

  d1=(ber+ahsqb)/f;
  d2=ber*asq/(f+2);
  d=d1-d2;
  sum=sum+d;
  t=ta-sum*bexp;
  aeps=ep1*t;
  ahsqb=ahsqb*a4b4/((g-1)*g);
  ab4=ab4*a4b4/((g+1)*g);
   f=f+4;
   g=g+2;
   if (d2*bexp<aeps) break; 
 }
   t=t*c;
 
  if (z<0) t=-t;
  return t;
} 
  
//*************************************************************************  
double sf54r(double x,double d,double idf) {
  double sys059,c,a,b,a2,a3,a1,sys089,sys029,tval,df,i1,sb,da,dsb,dasb,p1,f2,f1,sum,idfm2,az,fz,fkm1,p,z;
   int l;


sys059 = maxrealnumber;
c = 0.1591549;
a2 = 0.1591549;
a3 = 2.506628;
a1 = 0.3989423;
sys089 = 0.00001;
sys029 = 0.7071068;

if (idf<=0) return sys059;
df=idf;
tval=x;
i1=idf-2*round(idf/2);   
a=tval/sqrt(df);
b=df/(df+tval*tval);
sb =sqrt(b);
da=d*a;
dsb=d*sb;
dasb=a*dsb;
p1 = sf49r(dasb);
f2=a*sb*exp(-0.5*dsb*dsb)*p1*a1;
f1=b*(da*f2+a*a2*exp(-0.5*d*d));
sum = 0;
if (idf!=1.0) {
 if (i1>0) {
  sum = f1;
 } 
else {
  sum = f2;
}
if (idf>=4.) {
idfm2 = idf - 2;
az = 1.0;
fz = 2.0;
 for (l=2;l<=int(idfm2);l+=2) {
     fkm1 = fz - 1; 
     f2 = b * (da * az * f1 + f2) * fkm1 / fz;
      az = 1/(az * fkm1);
      f1 = b * (da * az * f2 + f1) * fz / (fz + 1);
      if (i1<=0) {
       sum = sum + f2;
      }
      else {
        sum = sum + f1;
      }
      az = 1/(az*fz);
      fz = fz + 2;
}
}
}
       if (i1>0) {
         p1=0.5 * sf35r(sys029*dsb);
         p=sf53r(dsb,a,sys089);
         z=p1+2*(p+sum);
        }
        else {
         p1=0.5*sf35r(sys029*d);
         z=p1+sum*a3;
       }
       
       if (z<0) return 0;
       return z;
} 

