/* C port of files in slatec_src.tgz, with time stamp from 1999 Nov 18:
     d9lgmc.f
     dcsevl.f
     dgamlm.f
     dgamma.f
     dlngam.f
     dlnrel.f
     initds.f
     d1mach.f

   http://cm.bell-labs.com/netlib/slatec/
 */

#include <math.h>
#include <stdio.h>
#include "lib_c.h"


#ifdef _WIN32
#define round(dbl) dbl >= 0.0 ? (int)(dbl + 0.5) : ((dbl - (double)(int)dbl) <= -0.5 ? (int)dbl : (int)(dbl - 0.5))
#endif

static const char* slatec_error_ptr = NULL;
static char slatec_error_buf[1024] = {'\0'};

static void
xermsg(
  const char* librar,
  const char* subrou,
  const char* messg,
  int nerr,
  int level)
{
  if (level < 2) return;
  if (slatec_error_ptr == NULL) {
    snprintf(slatec_error_buf, sizeof(slatec_error_buf), "%s: %s: %s (nerr=%d, level=%d)",
      librar, subrou, messg, nerr, level);
    slatec_error_ptr = slatec_error_buf;
  }
}

const char*
slatec_error(void) { return slatec_error_ptr; }

void
slatec_clear_error(void)
{
  slatec_error_buf[0] = '\0';
  slatec_error_ptr = NULL;
}

static double
min_double(double x, double y)
{
  if (x < y) return x;
  return y;
}

static double
max_double(double x, double y)
{
  if (x > y) return x;
  return y;
}

static int
initds(const double* os, int nos, double eta)
{
  int i, ii;
  double err;
  if (nos < 1) {
    xermsg("slatec", "initds",
      "number of coefficients is less than 1", 2, 1);
    return 0; /* to avoid use of uninitialized i */
  }
  err = 0.;
  for(ii=1;ii<=nos;ii++) {
    i = nos + 1 - ii;
    err = err + fabs(os[i-1]);
    if (err > eta) break;
  }
  if (i == nos) {
    xermsg("slatec", "initds",
      "chebyshev series too short for specified accuracy", 1, 1);
  }
  return i;
}

static double
d1mach(int i)
{
  double machine_constants[] = {
    2.2250738585072014e-308,
    1.7976931348623157e+308,
    1.1102230246251565e-16,
    2.2204460492503131e-16,
    0.301029995663981195
  };
  if (i < 1 || i > 5) {
    xermsg("slatec", "d1mach",
      "i out of bounds", 1, 2);
    return 1.0;
  }
  return machine_constants[i-1];
}

static double
dcsevl(double x, const double* cs, int n)
{
  static int first = 1;
  static double onepl;
  int i, ni;
  double b0, b1, b2, twox;
  if (first) {
    onepl = 1.0 + d1mach(4);
    first = 0;
  }
  if (n < 1) {
    xermsg("slatec", "dcsevl",
      "number of terms <= 0", 2, 2);
    return 1.0;
  }
  if (n > 1000) {
    xermsg("slatec", "dcsevl",
      "number of terms > 1000", 3, 2);
    return 1.0;
  }
  if (fabs(x) > onepl) {
    xermsg("slatec", "dcsevl",
      "x outside the interval (-1,+1)", 1, 1);
  }
  b1 = 0.0;
  b2 = 0.0;
  b0 = 0.0;
  twox = 2.0*x;
  for(i=1;i<=n;i++) {
    b2 = b1;
    b1 = b0;
    ni = n + 1 - i;
    b0 = twox*b1 - b2 + cs[ni-1];
  }
  return 0.5*(b0-b2);
}

static double
d9lgmc(double x)
{
  static int first=1, nalgm;
  static double xbig, xmax;
  static const double algmcs [] = {
     .1666389480451863247205729650822e+0,
    -.1384948176067563840732986059135e-4,
     .9810825646924729426157171547487e-8,
    -.1809129475572494194263306266719e-10,
     .6221098041892605227126015543416e-13,
    -.3399615005417721944303330599666e-15,
     .2683181998482698748957538846666e-17,
    -.2868042435334643284144622399999e-19,
     .3962837061046434803679306666666e-21,
    -.6831888753985766870111999999999e-23,
     .1429227355942498147573333333333e-24,
    -.3547598158101070547199999999999e-26,
     .1025680058010470912000000000000e-27,
    -.3401102254316748799999999999999e-29,
     .1276642195630062933333333333333e-30
  };
  if (first) {
    nalgm = initds(algmcs, 15, d1mach(3));
    xbig = 1.0/sqrt(d1mach(3));
    xmax = exp(min_double(log(d1mach(2)/12.), -log(12.*d1mach(1))));
    first = 0;
  }
  if (x < 10.) {
    xermsg("slatec", "d9lgmc",
      "x must be ge 10", 1, 2);
    return 1.0;
  }
  if (x < xmax) {
    if (x < xbig) {
      return dcsevl(2.0*(100./(x*x))-1., algmcs, nalgm) / x;
    }
    return 1./(12.*x);
  }
  xermsg("slatec", "d9lgmc",
    "x so big d9lgmc underflows", 2, 1);
  return 0.;
}

static void
dgamlm(double* xmin, double* xmax)
{
  int i;
  double alnbig, alnsml, xln, xold;
  alnsml = log(d1mach(1));
  (*xmin) = -alnsml;
  for(i=1;i<=10;i++) {
    xold = (*xmin);
    xln = log((*xmin));
    (*xmin) = (*xmin) - (*xmin)*(((*xmin)+0.5)*xln - (*xmin) - 0.2258 + alnsml)
      / ((*xmin)*xln+0.5);
    if (fabs((*xmin)-xold)<0.005) goto lbl_20;
  }
  xermsg("slatec", "dgamlm",
    "unable to find xmin", 1, 2);
  return;
 lbl_20:
  (*xmin) = -(*xmin) + 0.01;
  alnbig = log(d1mach(2));
  (*xmax) = alnbig;
  for(i=1;i<=10;i++) {
    xold = (*xmax);
    xln = log((*xmax));
    (*xmax) = (*xmax) - (*xmax)*(((*xmax)-0.5)*xln - (*xmax) + 0.9189 - alnbig)
      / ((*xmax)*xln-0.5);
    if (fabs((*xmax)-xold)<0.005) goto lbl_40;
  }
  xermsg("slatec", "dgamlm",
    "unable to find xmax", 2, 2);
  return;
 lbl_40:
  (*xmax) = (*xmax) - 0.01;
  (*xmin) = max_double((*xmin), -(*xmax)+1.);
}

double
slatec_dgamma(double x)
{
  static int first=1, ngam;
  static double xmin, xmax, dxrel;
  static const double pi = 3.14159265358979323846264338327950;
  static const double sq2pil = 0.91893853320467274178032973640562;
  static const double gamcs[] = {
     .8571195590989331421920062399942e-2,
     .4415381324841006757191315771652e-2,
     .5685043681599363378632664588789e-1,
    -.4219835396418560501012500186624e-2,
     .1326808181212460220584006796352e-2,
    -.1893024529798880432523947023886e-3,
     .3606925327441245256578082217225e-4,
    -.6056761904460864218485548290365e-5,
     .1055829546302283344731823509093e-5,
    -.1811967365542384048291855891166e-6,
     .3117724964715322277790254593169e-7,
    -.5354219639019687140874081024347e-8,
     .9193275519859588946887786825940e-9,
    -.1577941280288339761767423273953e-9,
     .2707980622934954543266540433089e-10,
    -.4646818653825730144081661058933e-11,
     .7973350192007419656460767175359e-12,
    -.1368078209830916025799499172309e-12,
     .2347319486563800657233471771688e-13,
    -.4027432614949066932766570534699e-14,
     .6910051747372100912138336975257e-15,
    -.1185584500221992907052387126192e-15,
     .2034148542496373955201026051932e-16,
    -.3490054341717405849274012949108e-17,
     .5987993856485305567135051066026e-18,
    -.1027378057872228074490069778431e-18,
     .1762702816060529824942759660748e-19,
    -.3024320653735306260958772112042e-20,
     .5188914660218397839717833550506e-21,
    -.8902770842456576692449251601066e-22,
     .1527474068493342602274596891306e-22,
    -.2620731256187362900257328332799e-23,
     .4496464047830538670331046570666e-24,
    -.7714712731336877911703901525333e-25,
     .1323635453126044036486572714666e-25,
    -.2270999412942928816702313813333e-26,
     .3896418998003991449320816639999e-27,
    -.6685198115125953327792127999999e-28,
     .1146998663140024384347613866666e-28,
    -.1967938586345134677295103999999e-29,
     .3376448816585338090334890666666e-30,
    -.5793070335782135784625493333333e-31
  };
  int i, n;
  double result, sinpiy, y;
  if (first) {
    ngam = initds(gamcs, 42, 0.1*d1mach(3));
    dgamlm(&xmin, &xmax);
    dxrel = sqrt(d1mach(4));
    first = 0;
  }
  y = fabs(x);
  if (y > 10.) goto lbl_50;
  n = (long) x;
  if (x < 0.) n = n - 1;
  y = x - n;
  n = n - 1;
  result = 0.9375 + dcsevl(2.*y-1., gamcs, ngam);
  if (n == 0) {
    return result;
  }
  if (n > 0) goto lbl_30;
  n = -n;
  if (x == 0.) {
    xermsg("slatec", "dgamma",
      "x is 0", 4, 2);
    return 1.0;
  }
  if (x < 0.0 && x+n-2 == 0.) {
    xermsg("slatec", "dgamma",
      "x is a negative integer", 4, 2);
    return 1.0;
  }
  if (x < (-0.5) && fabs((x-(long)(x-0.5))/x) < dxrel) {
    xermsg("slatec", "dgamma",
      "answer lt half precision because x too near negative integer", 1, 1);
  }
  for(i=1;i<=n;i++) {
    result = result/(x+i-1);
  }
  return result;
 lbl_30:
  for(i=1;i<=n;i++) {
    result = (y+i) * result;
  }
  return result;
 lbl_50:
  if (x > xmax) {
    xermsg("slatec", "dgamma",
      "x so big gamma overflows", 3, 2);
    return 1.0;
  }
  result = 0.;
  if (x < xmin) {
    xermsg("slatec", "dgamma",
      "x so small gamma underflows", 2, 1);
  }
  if (x < xmin) return result;
  result = exp((y-0.5)*log(y) - y + sq2pil + d9lgmc(y));
  if (x > 0.) return result;
  if (fabs((x-(long)(x-0.5))/x) < dxrel) {
    xermsg("slatec", "dgamma",
      "answer lt half precision, x too near negative int", 1, 1);
  }
  sinpiy = sin(pi*y);
  if (sinpiy == 0.) {
    xermsg("slatec", "dgamma",
      "x is a negative integer", 4, 2);
    return 1.0;
  }
  result = -pi/(y*sinpiy*result);
  return result;
}

double
slatec_dlngam(double x)
{
  static int first=1;
  static double xmax, dxrel;
  static const double sq2pil = 0.91893853320467274178032973640562;
  static const double sqpi2l = .225791352644727432363097614947441;
  static const double pi = 3.14159265358979323846264338327950;
  double result, sinpiy, temp, y;
  result = 0.; // silence compiler warnings
  if (first) {
    temp = 1./log(d1mach(2));
    xmax = temp*d1mach(2);
    dxrel = sqrt(d1mach(4));
    first = 0;
  }
  y = fabs(x);
  if (y > 10.) goto lbl_20;
  result = log(fabs(slatec_dgamma(x)));
  return result;
 lbl_20:
  if (y > xmax) {
    xermsg("slatec", "dlngam",
      "abs(x) so big dlngam overflows", 2, 2);
    return 1.0;
  }
  if (x>0.) result = sq2pil + (x-0.5)*log(x) - x + d9lgmc(y);
  if (x>0.) return result;
  sinpiy = fabs(sin(pi*fmod(y,2.0)));
  if (sinpiy == 0.) {
    xermsg("slatec", "dlngam",
      "x is a negative integer", 3, 2);
    return 1.0;
  }
  if (fabs((x-round(x))/x) < dxrel) {
    xermsg("slatec", "dlngam",
      "answer lt half precision because x too near negative int", 1, 1);
  }
  result = sqpi2l + (x-0.5)*log(y) - x - log(sinpiy) - d9lgmc(y);
  return result;
}

/*! \brief Natural logarithm of (1.0+x).
 */
/*! This routine should be used when X is small and accurate to
    calculate the logarithm accurately (in the relative error sense)
    in the neighborhood of 1.0.
 */
double
slatec_dlnrel(double x)
{
  static const double alnrcs[] = {
    +.10378693562743769800686267719098e+1,
    -.13364301504908918098766041553133e+0,
    +.19408249135520563357926199374750e-1,
    -.30107551127535777690376537776592e-2,
    +.48694614797154850090456366509137e-3,
    -.81054881893175356066809943008622e-4,
    +.13778847799559524782938251496059e-4,
    -.23802210894358970251369992914935e-5,
    +.41640416213865183476391859901989e-6,
    -.73595828378075994984266837031998e-7,
    +.13117611876241674949152294345011e-7,
    -.23546709317742425136696092330175e-8,
    +.42522773276034997775638052962567e-9,
    -.77190894134840796826108107493300e-10,
    +.14075746481359069909215356472191e-10,
    -.25769072058024680627537078627584e-11,
    +.47342406666294421849154395005938e-12,
    -.87249012674742641745301263292675e-13,
    +.16124614902740551465739833119115e-13,
    -.29875652015665773006710792416815e-14,
    +.55480701209082887983041321697279e-15,
    -.10324619158271569595141333961932e-15,
    +.19250239203049851177878503244868e-16,
    -.35955073465265150011189707844266e-17,
    +.67264542537876857892194574226773e-18,
    -.12602624168735219252082425637546e-18,
    +.23644884408606210044916158955519e-19,
    -.44419377050807936898878389179733e-20,
    +.83546594464034259016241293994666e-21,
    -.15731559416479562574899253521066e-21,
    +.29653128740247422686154369706666e-22,
    -.55949583481815947292156013226666e-23,
    +.10566354268835681048187284138666e-23,
    -.19972483680670204548314999466666e-24,
    +.37782977818839361421049855999999e-25,
    -.71531586889081740345038165333333e-26,
    +.13552488463674213646502024533333e-26,
    -.25694673048487567430079829333333e-27,
    +.48747756066216949076459519999999e-28,
    -.92542112530849715321132373333333e-29,
    +.17578597841760239233269760000000e-29,
    -.33410026677731010351377066666666e-30,
    +.63533936180236187354180266666666e-31
  };
  static int first = 1;
  static int nlnrel;
  static double xmin;
  if (first) {
    first = 0;
    nlnrel = initds(alnrcs, 43, 0.1*d1mach(3));
    xmin = -1 + sqrt(d1mach(4));
  }
  if (x <= -1) {
    xermsg("slatec", "dlnrel",
      "x is le -1" , 2, 2);
    return 1.0;
  }
  if (x < xmin) {
    xermsg("slatec", "dlnrel",
      "answer lt half precision because x too near -1", 1, 1);
  }
  if (fabs(x) <= 0.375) {
    return x*(1 - x*dcsevl(x/.375, alnrcs, nlnrel));
  }
  return log(1+x);
}

/*! \brief Binomial coefficient for integer arguments n and m.
 */
/*! The result is (n!)/((m!)(n-m)!).
 */
double
slatec_dbinom(unsigned n, unsigned m)
{
  unsigned i, k;
  double result, xn, xk, xnk, corr;
  static const double sq2pil = 0.91893853320467274178032973640562;
  static int first  = 1;
  static double bilnmx;
  static double fintmx;
  if (first) {
    first = 0;
    bilnmx = log(d1mach(2)) - 0.0001;
    fintmx = 0.9/d1mach(3);
  }
  if (n < m) {
    xermsg("slatec", "dbinom",
      "n lt m", 2, 2);
    return 1.0;
  }
  k = n-m;
  if (k > m) k = m;
  if (k <= 20) {
    if (k*log((double)(n > 1 ? n : 1)) <= bilnmx) {
      result = 1;
      if (k == 0) return result;
      for(i=1;i<=k;i++) {
        xn = n - i + 1;
        xk = i;
        result = result * (xn/xk);
      }
      if (result < fintmx) result = floor(result+0.5);
      return result;
    }
  }
  if (k < 9){ /* approx is not valid and answer is close to the overflow lim */
    xermsg("slatec", "dbinom",
      "result overflows because n and/or m too big", 3, 2);
    return 1.0;
  }
  xn = (double) n + 1;
  xk = (double) k + 1;
  xnk = (double) n - (double) k + 1;
  corr = d9lgmc(xn) - d9lgmc(xk) - d9lgmc(xnk);
  result = xk*log(xnk/xk) - xn*slatec_dlnrel(-(xk-1)/xn)
    -0.5*log(xn*xnk/xk) + 1 - sq2pil + corr;
  if (result > bilnmx) {
    xermsg("slatec", "dbinom",
      "result overflows because n and/or m too big", 3, 2);
    return 1.0;
  }
  result = exp(result);
  if (result < fintmx) result = floor(result+0.5);
  return result;
}
