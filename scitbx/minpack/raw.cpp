#include <scitbx/minpack/raw.h>

namespace scitbx { namespace minpack { namespace raw {

  double
  dpmpar(int i)
  {
    SCITBX_ASSERT(1 <= i && i <= 3);
    if (i == 1) return 2.22044604926e-16;
    if (i == 2) {return 2.22507385852e-308;}
    else        {return 1.79769313485e+308;}
  }

  /* given an n-vector x, this function calculates the
     euclidean norm of x.

     n is a positive integer input variable.

     x is an input array of length n.

     the euclidean norm is computed by accumulating the sum of
     squares in three different sums. the sums of squares for the
     small and large components are scaled so that no overflows
     occur. non-destructive underflows are permitted. underflows
     and overflows do not occur in the computation of the unscaled
     sum of squares for the intermediate components.
     the definitions of small, intermediate and large components
     depend on two constants, rdwarf and rgiant. the main
     restrictions on these constants are that rdwarf**2 not
     underflow and rgiant**2 not overflow. the constants
     given here are suitable for every known computer.

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more
   */
  double
  enorm(int n, ref1<double> const& x)
  {
    SCITBX_ASSERT(x.n >= n);
    static const double one = 1;
    static const double zero = 0;
    static const double rdwarf = 3.834e-20;
    static const double rgiant = 1.304e19;
    double s1 = zero;
    double s2 = zero;
    double s3 = zero;
    double x1max = zero;
    double x3max = zero;
    double floatn = n;
    SCITBX_ASSERT(floatn != 0);
    double agiant = rgiant/floatn;
    for(int i=1;i<=n;i++) {
      double xabs = std::fabs(x(i));
      if (xabs > rdwarf && xabs < agiant) goto lbl_70;
        if (xabs <= rdwarf) goto lbl_30;
          //
          // sum for large components.
          //
          if (xabs <= x1max) goto lbl_10;
            s1 = one + pow2(s1*(x1max/xabs));
            x1max = xabs;
            goto lbl_20;
          lbl_10:
            s1 = s1 + pow2((xabs/x1max));
          lbl_20:
          goto lbl_60;
        lbl_30:
          //
          // sum for small components.
          //
          if (xabs <= x3max) goto lbl_40;
            s3 = one + pow2(s3*(x3max/xabs));
            x3max = xabs;
            goto lbl_50;
          lbl_40:
            if (xabs != zero) s3 = s3 + pow2((xabs/x3max));
          lbl_50:
        lbl_60:
        goto lbl_80;
      lbl_70:
        //
        // sum for intermediate components.
        //
        s2 = s2 + xabs*xabs;
      lbl_80:;
    }
    //
    // calculation of norm.
    //
    double result;
    if (s1 == zero) goto lbl_100;
      result = x1max*std::sqrt(s1+(s2/x1max)/x1max);
      goto lbl_130;
    lbl_100:
      if (s2 == zero) goto lbl_110;
        if (s2 >= x3max) {
          result = std::sqrt(s2*(one+(x3max/s2)*(x3max*s3)));
        }
        else {
          result = std::sqrt(x3max*((s2/x3max)+(x3max*s3)));
        }
        goto lbl_120;
      lbl_110:
        result = x3max*std::sqrt(s3);
      lbl_120:
    lbl_130:
    return result;
  }

  /* this subroutine uses householder transformations with column
     pivoting (optional) to compute a qr factorization of the
     m by n matrix a. that is, qrfac determines an orthogonal
     matrix q, a permutation matrix p, and an upper trapezoidal
     matrix r with diagonal elements of nonincreasing magnitude,
     such that a*p = q*r. the householder transformation for
     column k, k = 1,2,...,min(m,n), is of the form

                           t
           i - (1/u(k))*u*u

     where u has zeros in the first k-1 positions. the form of
     this transformation and the method of pivoting first
     appeared in the corresponding linpack subroutine.

       m is a positive integer input variable set to the number
         of rows of a.

       n is a positive integer input variable set to the number
         of columns of a.

       a is an m by n array. on input a contains the matrix for
         which the qr factorization is to be computed. on output
         the strict upper trapezoidal part of a contains the strict
         upper trapezoidal part of r, and the lower trapezoidal
         part of a contains a factored form of q (the non-trivial
         elements of the u vectors described above).

       lda is a positive integer input variable not less than m
         which specifies the leading dimension of the array a.

       pivot is a logical input variable. if pivot is set true,
         then column pivoting is enforced. if pivot is set false,
         then no column pivoting is done.

       ipvt is an integer output array of length lipvt. ipvt
         defines the permutation matrix p such that a*p = q*r.
         column j of p is column ipvt(j) of the identity matrix.
         if pivot is false, ipvt is not referenced.

       lipvt is a positive integer input variable. if pivot is false,
         then lipvt may be as small as 1. if pivot is true, then
         lipvt must be at least n.

       rdiag is an output array of length n which contains the
         diagonal elements of r.

       acnorm is an output array of length n which contains the
         norms of the corresponding columns of the input matrix a.
         if this information is not needed, then acnorm can coincide
         with rdiag.

       wa is a work array of length n. if pivot is false, then wa
         can coincide with rdiag.

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more
   */
  void
  qrfac(
    const int m,
    const int n,
    ref2<double> const& a,
    const int lda,
    const bool pivot,
    ref1<int> const& ipvt,
    const int /*lipvt*/,
    ref1<double> const& rdiag,
    ref1<double> const& acnorm,
    ref1<double> const& wa)
  {
    SCITBX_ASSERT(a.ni == lda);
    SCITBX_ASSERT(a.nj == n);
    static const double one = 1;
    static const double p05 = 5.0e-2;
    static const double zero = 0;
    //
    // epsmch is the machine precision.
    //
    double epsmch = dpmpar(1);
    //
    // compute the initial column norms and initialize several arrays.
    //
    for(int j=1;j<=n;j++) {
      acnorm(j) = enorm(m, ref1<double>(&a(1,j), m));
      rdiag(j) = acnorm(j);
      wa(j) = rdiag(j);
      if (pivot) ipvt(j) = j;
    }
    //
    // reduce a to r with householder transformations.
    //
    int minmn = std::min(m, n);
    for(int j=1;j<=minmn;j++) {
      if (!pivot) goto lbl_40;
      {
        //
        // bring the column of largest norm into the pivot position.
        //
        int kmax = j;
        for(int k=j;k<=n;k++) {
          if (rdiag(k) > rdiag(kmax)) kmax = k;
        }
        if (kmax == j) goto lbl_40;
        {
          for(int i=1;i<=m;i++) {
            double temp = a(i,j);
            a(i,j) = a(i,kmax);
            a(i,kmax) = temp;
          }
          rdiag(kmax) = rdiag(j);
          wa(kmax) = wa(j);
          int k = ipvt(j);
          ipvt(j) = ipvt(kmax);
          ipvt(kmax) = k;
        }
      }
      lbl_40:
      //
      // compute the householder transformation to reduce the
      // j-th column of a to a multiple of the j-th unit vector.
      //
      double ajnorm = enorm(m-j+1,ref1<double>(&a(j,j), m-j+1));
      if (ajnorm == zero) goto lbl_100;
      {
        if (a(j,j) < zero) ajnorm = -ajnorm;
        for(int i=j;i<=m;i++) {
          a(i,j) = a(i,j)/ajnorm;
        }
        a(j,j) = a(j,j) + one;
        //
        // apply the transformation to the remaining columns
        // and update the norms.
        //
        int jp1 = j + 1;
        if (n < jp1) goto lbl_100;
        for(int k=jp1;k<=n;k++) {
          double sum = zero;
          for(int i=j;i<=m;i++) {
            sum = sum + a(i,j)*a(i,k);
          }
          double temp = sum/a(j,j);
          for(int i=j;i<=m;i++) {
            a(i,k) = a(i,k) - temp*a(i,j);
          }
          if (!pivot || rdiag(k) == zero) goto lbl_80;
          temp = a(j,k)/rdiag(k);
          rdiag(k) = rdiag(k)*std::sqrt(std::max(zero,one-pow2(temp)));
          if (p05*pow2(rdiag(k)/wa(k)) > epsmch) goto lbl_80;
          rdiag(k) = enorm(m-j,ref1<double>(&a(jp1,k),m-j));
          wa(k) = rdiag(k);
          lbl_80:;
        }
      }
      lbl_100:
      rdiag(j) = -ajnorm;
    }
  }

  /* given an m by n matrix a, an n by n diagonal matrix d,
     and an m-vector b, the problem is to determine an x which
     solves the system

           a*x = b ,     d*x = 0 ,

     in the least squares sense.

     this subroutine completes the solution of the problem
     if it is provided with the necessary information from the
     qr factorization, with column pivoting, of a. that is, if
     a*p = q*r, where p is a permutation matrix, q has orthogonal
     columns, and r is an upper triangular matrix with diagonal
     elements of nonincreasing magnitude, then qrsolv expects
     the full upper triangle of r, the permutation matrix p,
     and the first n components of (q transpose)*b. the system
     a*x = b, d*x = 0, is then equivalent to

                  t       t
           r*z = q *b ,  p *d*p*z = 0 ,

     where x = p*z. if this system does not have full rank,
     then a least squares solution is obtained. on output qrsolv
     also provides an upper triangular matrix s such that

            t   t               t
           p *(a *a + d*d)*p = s *s .

     s is computed within qrsolv and may be of separate interest.

       n is a positive integer input variable set to the order of r.

       r is an n by n array. on input the full upper triangle
         must contain the full upper triangle of the matrix r.
         on output the full upper triangle is unaltered, and the
         strict lower triangle contains the strict upper triangle
         (transposed) of the upper triangular matrix s.

       ldr is a positive integer input variable not less than n
         which specifies the leading dimension of the array r.

       ipvt is an integer input array of length n which defines the
         permutation matrix p such that a*p = q*r. column j of p
         is column ipvt(j) of the identity matrix.

       diag is an input array of length n which must contain the
         diagonal elements of the matrix d.

       qtb is an input array of length n which must contain the first
         n elements of the vector (q transpose)*b.

       x is an output array of length n which contains the least
         squares solution of the system a*x = b, d*x = 0.

       sdiag is an output array of length n which contains the
         diagonal elements of the upper triangular matrix s.

       wa is a work array of length n.

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more
   */
  void
  qrsolv(
    const int n,
    ref2<double> const& r,
    const int ldr,
    ref1<int> const& ipvt,
    ref1<double> const& diag,
    ref1<double> const& qtb,
    ref1<double> const& x,
    ref1<double> const& sdiag,
    ref1<double> const& wa)
  {
    SCITBX_ASSERT(ipvt.n == n);
    SCITBX_ASSERT(r.ni == ldr);
    SCITBX_ASSERT(r.nj == n);
    SCITBX_ASSERT(diag.n == n);
    SCITBX_ASSERT(qtb.n == n);
    SCITBX_ASSERT(x.n == n);
    SCITBX_ASSERT(sdiag.n == n);
    SCITBX_ASSERT(wa.n >= n);
    static const double p5 = 5.0e-1;
    static const double p25 = 2.5e-1;
    static const double zero = 0;
    //
    // copy r and (q transpose)*b to preserve input and initialize s.
    // in particular, save the diagonal elements of r in x.
    //
    for(int j=1;j<=n;j++) {
      for(int i=j;i<=n;i++) {
        r(i,j) = r(j,i);
      }
      x(j) = r(j,j);
      wa(j) = qtb(j);
    }
    //
    // eliminate the diagonal matrix d using a givens rotation.
    //
    for(int j=1;j<=n;j++) {
      //
      //  prepare the row of d to be eliminated, locating the
      //  diagonal element using p from the qr factorization.
      //
      int l = ipvt(j);
      if (diag(l) == zero) goto lbl_90;
      {
        for(int k=j;k<=n;k++) {
          sdiag(k) = zero;
        }
        sdiag(j) = diag(l);
        //
        // the transformations to eliminate the row of d
        // modify only a single element of (q transpose)*b
        // beyond the first n, which is initially zero.
        //
        double qtbpj = zero;
        for(int k=j;k<=n;k++) {
          //
          // determine a givens rotation which eliminates the
          // appropriate element in the current row of d.
          //
          if (sdiag(k) == zero) goto lbl_70;
          {
            double cos, sin;
            if (std::fabs(r(k,k)) >= std::fabs(sdiag(k))) goto lbl_40;
            {
              double cotan = r(k,k)/sdiag(k);
              sin = p5/std::sqrt(p25+p25*pow2(cotan));
              cos = sin*cotan;
              goto lbl_50;
            }
            lbl_40:
            {
              double tan = sdiag(k)/r(k,k);
              cos = p5/std::sqrt(p25+p25*pow2(tan));
              sin = cos*tan;
            }
            lbl_50:
            //
            // compute the modified diagonal element of r and
            // the modified element of ((q transpose)*b,0).
            //
            r(k,k) = cos*r(k,k) + sin*sdiag(k);
            double temp = cos*wa(k) + sin*qtbpj;
            qtbpj = -sin*wa(k) + cos*qtbpj;
            wa(k) = temp;
            //
            // accumulate the tranformation in the row of s.
            //
            int kp1 = k + 1;
            if (n < kp1) goto lbl_70;
            for(int i=kp1;i<=n;i++) {
              temp = cos*r(i,k) + sin*sdiag(i);
              sdiag(i) = -sin*r(i,k) + cos*sdiag(i);
              r(i,k) = temp;
            }
          }
          lbl_70:;
        }
      }
      lbl_90:
      //
      // store the diagonal element of s and restore
      // the corresponding diagonal element of r.
      //
      sdiag(j) = r(j,j);
      r(j,j) = x(j);
    }
    //
    // solve the triangular system for z. if the system is
    // singular, then obtain a least squares solution.
    //
    int nsing = n;
    for(int j=1;j<=n;j++) {
      if (sdiag(j) == zero && nsing == n) nsing = j - 1;
      if (nsing < n) wa(j) = zero;
    }
    if (nsing < 1) goto lbl_150;
    for(int k=1;k<=nsing;k++) {
      int j = nsing - k + 1;
      double sum = zero;
      int jp1 = j + 1;
      if (nsing < jp1) goto lbl_130;
      for(int i=jp1;i<=nsing;i++) {
        sum = sum + r(i,j)*wa(i);
      }
      lbl_130:
      wa(j) = (wa(j) - sum)/sdiag(j);
    }
    lbl_150:
    //
    // permute the components of z back to components of x.
    //
    for(int j=1;j<=n;j++) {
      int l = ipvt(j);
      x(l) = wa(j);
    }
  }

  /* given an m by n matrix a, an n by n nonsingular diagonal
     matrix d, an m-vector b, and a positive number delta,
     the problem is to determine a value for the parameter
     par such that if x solves the system

           a*x = b ,     sqrt(par)*d*x = 0 ,

     in the least squares sense, and dxnorm is the euclidean
     norm of d*x, then either par is zero and

           (dxnorm-delta) .le. 0.1*delta ,

     or par is positive and

           abs(dxnorm-delta) .le. 0.1*delta .

     this subroutine completes the solution of the problem
     if it is provided with the necessary information from the
     qr factorization, with column pivoting, of a. that is, if
     a*p = q*r, where p is a permutation matrix, q has orthogonal
     columns, and r is an upper triangular matrix with diagonal
     elements of nonincreasing magnitude, then lmpar expects
     the full upper triangle of r, the permutation matrix p,
     and the first n components of (q transpose)*b. on output
     lmpar also provides an upper triangular matrix s such that

            t   t                   t
           p *(a *a + par*d*d)*p = s *s .

     s is employed within lmpar and may be of separate interest.

     only a few iterations are generally needed for convergence
     of the algorithm. if, however, the limit of 10 iterations
     is reached, then the output par will contain the best
     value obtained so far.

       n is a positive integer input variable set to the order of r.

       r is an n by n array. on input the full upper triangle
         must contain the full upper triangle of the matrix r.
         on output the full upper triangle is unaltered, and the
         strict lower triangle contains the strict upper triangle
         (transposed) of the upper triangular matrix s.

       ldr is a positive integer input variable not less than n
         which specifies the leading dimension of the array r.

       ipvt is an integer input array of length n which defines the
         permutation matrix p such that a*p = q*r. column j of p
         is column ipvt(j) of the identity matrix.

       diag is an input array of length n which must contain the
         diagonal elements of the matrix d.

       qtb is an input array of length n which must contain the first
         n elements of the vector (q transpose)*b.

       delta is a positive input variable which specifies an upper
         bound on the euclidean norm of d*x.

       par is a nonnegative variable. on input par contains an
         initial estimate of the levenberg-marquardt parameter.
         on output par contains the final estimate.

       x is an output array of length n which contains the least
         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
         for the output par.

       sdiag is an output array of length n which contains the
         diagonal elements of the upper triangular matrix s.

       wa1 and wa2 are work arrays of length n.

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more
   */
  void
  lmpar(
    const int n,
    ref2<double> const& r,
    const int ldr,
    ref1<int> const& ipvt,
    ref1<double> const& diag,
    ref1<double> const& qtb,
    const double delta,
    double& par,
    ref1<double> const& x,
    ref1<double> const& sdiag,
    ref1<double> const& wa1,
    ref1<double> const& wa2)
  {
    SCITBX_ASSERT(r.ni == ldr);
    SCITBX_ASSERT(r.nj == n);
    SCITBX_ASSERT(ipvt.n == n);
    SCITBX_ASSERT(diag.n == n);
    SCITBX_ASSERT(qtb.n == n);
    SCITBX_ASSERT(delta != 0);
    SCITBX_ASSERT(x.n == n);
    SCITBX_ASSERT(sdiag.n == n);
    SCITBX_ASSERT(wa1.n == n);
    SCITBX_ASSERT(wa2.n >= n);
    static const double p1 = 1.0e-1;
    static const double p001 = 1.0e-3;
    static const double zero = 0;
    //
    // dwarf is the smallest positive magnitude.
    //
    static double dwarf = dpmpar(2);
    //
    // compute and store in x the gauss-newton direction. if the
    // jacobian is rank-deficient, obtain a least squares solution.
    //
    int nsing = n;
    for(int j=1;j<=n;j++) {
      wa1(j) = qtb(j);
      if (r(j,j) == zero && nsing == n) nsing = j - 1;
      if (nsing < n) wa1(j) = zero;
    }
    if (nsing < 1) goto lbl_50;
    for(int k=1;k<=nsing;k++) {
      int j = nsing - k + 1;
      wa1(j) = wa1(j)/r(j,j);
      double temp = wa1(j);
      int jm1 = j - 1;
      if (jm1 < 1) goto lbl_30;
      for(int i=1;i<=jm1;i++) {
        wa1(i) = wa1(i) - r(i,j)*temp;
      }
      lbl_30:;
    }
    lbl_50:
    for(int j=1;j<=n;j++) {
      int l = ipvt(j);
      x(l) = wa1(j);
    }
    //
    // initialize the iteration counter.
    // evaluate the function at the origin, and test
    // for acceptance of the gauss-newton direction.
    //
    int iter = 0;
    for(int j=1;j<=n;j++) {
      wa2(j) = diag(j)*x(j);
    }
    double dxnorm = enorm(n,wa2);
    double fp = dxnorm - delta;
    if (fp <= p1*delta) goto lbl_220;
    {
      //
      // if the jacobian is not rank deficient, the newton
      // step provides a lower bound, parl, for the zero of
      // the function. otherwise set this bound to zero.
      //
      double parl = zero;
      if (nsing < n) goto lbl_120;
      {
        for(int j=1;j<=n;j++) {
          int l = ipvt(j);
          wa1(j) = diag(l)*(wa2(l)/dxnorm);
        }
        for(int j=1;j<=n;j++) {
          double sum = zero;
          int jm1 = j - 1;
          if (jm1 < 1) goto lbl_100;
          for(int i=1;i<=jm1;i++) {
            sum = sum + r(i,j)*wa1(i);
          }
          lbl_100:
          wa1(j) = (wa1(j) - sum)/r(j,j);
        }
        double temp = enorm(n,wa1);
        parl = ((fp/delta)/temp)/temp;
      }
      lbl_120:
      //
      // calculate an upper bound, paru, for the zero of the function.
      //
      for(int j=1;j<=n;j++) {
        double sum = zero;
        for(int i=1;i<=j;i++) {
          sum = sum + r(i,j)*qtb(i);
        }
        int l = ipvt(j);
        wa1(j) = sum/diag(l);
      }
      double gnorm = enorm(n,wa1);
      double paru = gnorm/delta;
      if (paru == zero) paru = dwarf/std::min(delta,p1);
      //
      // if the input par lies outside of the interval (parl,paru),
      // set par to the closer endpoint.
      //
      par = std::max(par,parl);
      par = std::min(par,paru);
      if (par == zero) par = gnorm/dxnorm;
      //
      // beginning of an iteration.
      //
      lbl_150:
      {
        iter = iter + 1;
        //
        // evaluate the function at the current value of par.
        //
        if (par == zero) par = std::max(dwarf,p001*paru);
        double temp = std::sqrt(par);
        for(int j=1;j<=n;j++) {
          wa1(j) = temp*diag(j);
        }
        qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2);
        for(int j=1;j<=n;j++) {
          wa2(j) = diag(j)*x(j);
        }
        dxnorm = enorm(n,wa2);
        temp = fp;
        fp = dxnorm - delta;
        //
        // if the function is small enough, accept the current value
        // of par. also test for the exceptional cases where parl
        // is zero or the number of iterations has reached 10.
        //
        if (   std::fabs(fp) <= p1*delta
            || (parl == zero && fp <= temp && temp < zero)
            || iter == 10) goto lbl_220;
        //
        // compute the newton correction.
        //
        for(int j=1;j<=n;j++) {
          int l = ipvt(j);
          wa1(j) = diag(l)*(wa2(l)/dxnorm);
        }
        for(int j=1;j<=n;j++) {
          SCITBX_ASSERT(sdiag(j) != 0);
          wa1(j) = wa1(j)/sdiag(j);
          temp = wa1(j);
          int jp1 = j + 1;
          if (n < jp1) goto lbl_200;
          for(int i=jp1;i<=n;i++) {
            wa1(i) = wa1(i) - r(i,j)*temp;
          }
          lbl_200:;
        }
        temp = enorm(n,wa1);
        SCITBX_ASSERT(temp != 0);
        double parc = ((fp/delta)/temp)/temp;
        //
        // depending on the sign of the function, update parl or paru.
        //
        if (fp > zero) parl = std::max(parl,par);
        if (fp < zero) paru = std::min(paru,par);
        //
        // compute an improved estimate for par.
        //
        par = std::max(parl,par+parc);
        //
        // end of an iteration.
        //
        goto lbl_150;
      }
    }
    lbl_220:
    //
    // termination.
    //
    if (iter == 0) par = zero;
  }

}}} // namespace scitbx::minpack::raw
