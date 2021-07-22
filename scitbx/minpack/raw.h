#ifndef SCITBX_MINPACK_LIB_H
#define SCITBX_MINPACK_LIB_H

#include <scitbx/error.h>
#include <algorithm>
#include <cmath>

namespace scitbx {

//! C++ port of files in MINPACK
/*! Original MINPACK obtained 2005-08-28.
    http://www.netlib.org/minpack/
 */
namespace minpack {

//! Raw MINPACK interfaces with minimal adjustments.
namespace raw {

  template <typename NumType>
  NumType
  pow2(NumType const& x) { return x * x; }

  template <typename ElementType>
  struct ref1
  {
    ElementType* begin;
    int n;

    ref1(ElementType* begin_, int n_)
    :
      begin(begin_),
      n(n_)
    {}

    ElementType&
    operator()(int i) const { return begin[(i-1)]; }
  };

  template <typename ElementType>
  struct ref2
  {
    ElementType* begin;
    int ni, nj;

    ref2(ElementType* begin_, int ni_, int nj_)
    :
      begin(begin_),
      ni(ni_),
      nj(nj_)
    {}

    ElementType&
    operator()(int i, int j) const { return begin[(i-1) + (j-1)*ni]; }
  };

  double
  dpmpar(int i);

  double
  enorm(int n, ref1<double> const& x);

  void
  qrfac(
    const int m,
    const int n,
    ref2<double> const& a,
    const int lda,
    const bool pivot,
    ref1<int> const& ipvt,
    const int lipvt,
    ref1<double> const& rdiag,
    ref1<double> const& acnorm,
    ref1<double> const& wa);

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
    ref1<double> const& wa);

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
    ref1<double> const& wa2);

  //! C++ port of lmder.
  /*! the purpose of lmder is to minimize the sum of the squares of
      m nonlinear functions in n variables by a modification of
      the levenberg-marquardt algorithm. the user must provide a
      subroutine which calculates the functions and the jacobian.

        fcn is the name of the user-supplied subroutine which
          calculates the functions and the jacobian. fcn must
          be declared in an external statement in the user
          calling program, and should be written as follows.

          subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
          integer m,n,ldfjac,iflag
          double precision x(n),fvec(m),fjac(ldfjac,n)
          ----------
          if iflag = 1 calculate the functions at x and
          return this vector in fvec. do not alter fjac.
          if iflag = 2 calculate the jacobian at x and
          return this matrix in fjac. do not alter fvec.
          ----------
          return
          end

          the value of iflag should not be changed by fcn unless
          the user wants to terminate execution of lmder.
          in this case set iflag to a negative integer.

        m is a positive integer input variable set to the number
          of functions.

        n is a positive integer input variable set to the number
          of variables. n must not exceed m.

        x is an array of length n. on input x must contain
          an initial estimate of the solution vector. on output x
          contains the final estimate of the solution vector.

        fvec is an output array of length m which contains
          the functions evaluated at the output x.

        fjac is an output m by n array. the upper n by n submatrix
          of fjac contains an upper triangular matrix r with
          diagonal elements of nonincreasing magnitude such that

                 t     t           t
                p *(jac *jac)*p = r *r,

          where p is a permutation matrix and jac is the final
          calculated jacobian. column j of p is column ipvt(j)
          (see below) of the identity matrix. the lower trapezoidal
          part of fjac contains information generated during
          the computation of r.

        ldfjac is a positive integer input variable not less than m
          which specifies the leading dimension of the array fjac.

        ftol is a nonnegative input variable. termination
          occurs when both the actual and predicted relative
          reductions in the sum of squares are at most ftol.
          therefore, ftol measures the relative error desired
          in the sum of squares.

        xtol is a nonnegative input variable. termination
          occurs when the relative error between two consecutive
          iterates is at most xtol. therefore, xtol measures the
          relative error desired in the approximate solution.

        gtol is a nonnegative input variable. termination
          occurs when the cosine of the angle between fvec and
          any column of the jacobian is at most gtol in absolute
          value. therefore, gtol measures the orthogonality
          desired between the function vector and the columns
          of the jacobian.

        maxfev is a positive integer input variable. termination
          occurs when the number of calls to fcn with iflag = 1
          has reached maxfev.

        diag is an array of length n. if mode = 1 (see
          below), diag is internally set. if mode = 2, diag
          must contain positive entries that serve as
          multiplicative scale factors for the variables.

        mode is an integer input variable. if mode = 1, the
          variables will be scaled internally. if mode = 2,
          the scaling is specified by the input diag. other
          values of mode are equivalent to mode = 1.

        factor is a positive input variable used in determining the
          initial step bound. this bound is set to the product of
          factor and the euclidean norm of diag*x if nonzero, or else
          to factor itself. in most cases factor should lie in the
          interval (.1,100.). 100. is a generally recommended value.

        info is an integer output variable. if the user has
          terminated execution, info is set to the (negative)
          value of iflag. see description of fcn. otherwise,
          info is set as follows.

          info = 0  improper input parameters.

          info = 1  both actual and predicted relative reductions
                    in the sum of squares are at most ftol.

          info = 2  relative error between two consecutive iterates
                    is at most xtol.

          info = 3  conditions for info = 1 and info = 2 both hold.

          info = 4  the cosine of the angle between fvec and any
                    column of the jacobian is at most gtol in
                    absolute value.

          info = 5  number of calls to fcn with iflag = 1 has
                    reached maxfev.

          info = 6  ftol is too small. no further reduction in
                    the sum of squares is possible.

          info = 7  xtol is too small. no further improvement in
                    the approximate solution x is possible.

          info = 8  gtol is too small. fvec is orthogonal to the
                    columns of the jacobian to machine precision.

        nfev is an integer output variable set to the number of
          calls to fcn with iflag = 1.

        njev is an integer output variable set to the number of
          calls to fcn with iflag = 2.

        ipvt is an integer output array of length n. ipvt
          defines a permutation matrix p such that jac*p = q*r,
          where jac is the final calculated jacobian, q is
          orthogonal (not stored), and r is upper triangular
          with diagonal elements of nonincreasing magnitude.
          column j of p is column ipvt(j) of the identity matrix.

        qtf is an output array of length n which contains
          the first n elements of the vector (q transpose)*fvec.

        wa1, wa2, and wa3 are work arrays of length n.

        wa4 is a work array of length m.

      argonne national laboratory. minpack project. march 1980.
      burton s. garbow, kenneth e. hillstrom, jorge j. more
   */
  struct lmder
  {
    int iflag;
    bool call_back_after_iteration;
    double fnorm;
    double par;
    int iter;
    double xnorm;
    double delta;
    double gnorm;
    double pnorm;
    double fnorm1;
    double actred;
    double prered;
    double dirder;
    double ratio;

    lmder(bool call_back_after_iteration_=false)
    :
      iflag(0),
      call_back_after_iteration(call_back_after_iteration_)
    {}

    void
    run(
      const int m,
      const int n,
      ref1<double> const& x,
      ref1<double> const& fvec,
      ref2<double> const& fjac,
      const int ldfjac,
      const double ftol,
      const double xtol,
      const double gtol,
      const int maxfev,
      ref1<double> const& diag,
      const int mode,
      const double factor,
      int& info,
      int& nfev,
      int& njev,
      ref1<int> const& ipvt,
      ref1<double> const& qtf,
      ref1<double> const& wa1,
      ref1<double> const& wa2,
      ref1<double> const& wa3,
      ref1<double> const& wa4)
    {
      SCITBX_ASSERT(x.n == n);
      SCITBX_ASSERT(fvec.n == m);
      SCITBX_ASSERT(ipvt.n == n);
      SCITBX_ASSERT(fjac.ni == ldfjac);
      SCITBX_ASSERT(fjac.nj == n);
      SCITBX_ASSERT(diag.n == n);
      SCITBX_ASSERT(qtf.n == n);
      SCITBX_ASSERT(wa1.n == n);
      SCITBX_ASSERT(wa2.n == n);
      SCITBX_ASSERT(wa3.n == n);
      SCITBX_ASSERT(wa4.n == m);
      static const double zero = 0.0e0;
      static const double one = 1.0e0;
      static const double p1 = 1.0e-1;
      static const double p5 = 5.0e-1;
      static const double p25 = 2.5e-1;
      static const double p75 = 7.5e-1;
      static const double p0001 = 1.0e-4;
      //
      SCITBX_ASSERT(0 <= iflag && iflag <= 4);
      if (iflag == 1) {
        goto continue_after_evaluate_function_at_starting_point;
      }
      if (iflag == 2) {
        goto continue_after_calculate_jacobian;
      }
      if (iflag == 3) {
        goto continue_after_evaluate_function_at_x_plus_p;
      }
      if (iflag == 4) {
        goto continue_after_call_back_after_iteration;
      }
      //
      // epsmch is the machine precision.
      //
      static const double epsmch = dpmpar(1);
      //
      info = 0;
      nfev = 0;
      njev = 0;
      //
      // check the input parameters for errors.
      //
      if (   n <= 0 || m < n || ldfjac < m
          || ftol < zero || xtol < zero || gtol < zero
          || maxfev < 0 || factor <= zero) goto lbl_300;
      if (mode != 2) goto lbl_20;
      for(int j=1;j<=n;j++) {
        if (diag(j) <= zero) goto lbl_300;
      }
      lbl_20:
      //
      // evaluate the function at the starting point
      // and calculate its norm.
      //
      iflag = 1;
      return;
     continue_after_evaluate_function_at_starting_point:
      nfev = 1;
      fnorm = enorm(m,fvec);
      //
      // initialize levenberg-marquardt parameter and iteration counter.
      //
      par = zero;
      iter = 1;
      //
      // beginning of the outer loop.
      //
      lbl_30:
        //
        // calculate the jacobian matrix.
        //
        iflag = 2;
        return;
       continue_after_calculate_jacobian:
        njev = njev + 1;
        //
        // compute the qr factorization of the jacobian.
        //
        qrfac(m,n,fjac,ldfjac,true,ipvt,n,wa1,wa2,wa3);
        //
        // on the first iteration and if mode is 1, scale according
        // to the norms of the columns of the initial jacobian.
        //
        if (iter != 1) goto lbl_80;
        if (mode == 2) goto lbl_60;
        for(int j=1;j<=n;j++) {
          diag(j) = wa2(j);
          if (wa2(j) == zero) diag(j) = one;
        }
        lbl_60:
        //
        // on the first iteration, calculate the norm of the scaled x
        // and initialize the step bound delta.
        //
        for(int j=1;j<=n;j++) {
          wa3(j) = diag(j)*x(j);
        }
        xnorm = enorm(n,wa3);
        delta = factor*xnorm;
        if (delta == zero) delta = factor;
        lbl_80:
        //
        // form (q transpose)*fvec and store the first n components in
        // qtf.
        //
        for(int i=1;i<=m;i++) {
          wa4(i) = fvec(i);
        }
        for(int j=1;j<=n;j++) {
          if (fjac(j,j) == zero) goto lbl_120;
          {
            double sum = zero;
            for(int i=j;i<=m;i++) {
              sum = sum + fjac(i,j)*wa4(i);
            }
            double temp = -sum/fjac(j,j);
            for(int i=j;i<=m;i++) {
              wa4(i) = wa4(i) + fjac(i,j)*temp;
            }
          }
          lbl_120:
          fjac(j,j) = wa1(j);
          qtf(j) = wa4(j);
        }
        //
        // compute the norm of the scaled gradient.
        //
        gnorm = zero;
        if (fnorm == zero) goto lbl_170;
        for(int j=1;j<=n;j++) {
          int l = ipvt(j);
          if (wa2(l) == zero) goto lbl_150;
          {
            double sum = zero;
            for(int i=1;i<=j;i++) {
              sum = sum + fjac(i,j)*(qtf(i)/fnorm);
            }
            gnorm = std::max(gnorm,std::fabs(sum/wa2(l)));
          }
          lbl_150:;
        }
        lbl_170:
        //
        // test for convergence of the gradient norm.
        //
        if (gnorm <= gtol) info = 4;
        if (info != 0) goto lbl_300;
        //
        // rescale if necessary.
        //
        if (mode == 2) goto lbl_190;
        for(int j=1;j<=n;j++) {
          diag(j) = std::max(diag(j),wa2(j));
        }
        lbl_190:
        //
        // beginning of the inner loop.
        //
        lbl_200:
          //
          // determine the levenberg-marquardt parameter.
          //
          lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,wa3,wa4);
          //
          // store the direction p and x + p. calculate the norm of p.
          //
          for(int j=1;j<=n;j++) {
            wa1(j) = -wa1(j);
            wa2(j) = x(j) + wa1(j);
            wa3(j) = diag(j)*wa1(j);
          }
          pnorm = enorm(n,wa3);
          //
          // on the first iteration, adjust the initial step bound.
          //
          if (iter == 1) delta = std::min(delta,pnorm);
          //
          // evaluate the function at x + p and calculate its norm.
          //
          iflag = 3;
          return;
         continue_after_evaluate_function_at_x_plus_p:
          nfev = nfev + 1;
          fnorm1 = enorm(m,wa4);
          //
          // compute the scaled actual reduction.
          //
          actred = -one;
          if (p1*fnorm1 < fnorm) actred = one - pow2(fnorm1/fnorm);
          //
          // compute the scaled predicted reduction and
          // the scaled directional derivative.
          //
          for(int j=1;j<=n;j++) {
            wa3(j) = zero;
            int l = ipvt(j);
            double temp = wa1(l);
            for(int i=1;i<=j;i++) {
              wa3(i) = wa3(i) + fjac(i,j)*temp;
            }
          }
          {
            double temp1 = enorm(n,wa3)/fnorm;
            double temp2 = (std::sqrt(par)*pnorm)/fnorm;
            prered = pow2(temp1) + pow2(temp2)/p5;
            dirder = -(pow2(temp1) + pow2(temp2));
          }
          //
          // compute the ratio of the actual to the predicted
          // reduction.
          //
          ratio = zero;
          if (prered != zero) ratio = actred/prered;
          //
          // update the step bound.
          //
          if (ratio > p25) goto lbl_240;
          // explicitly dedenting this bracket to remove a misleading-indentation compiler warning:
          {
            double temp;
            if (actred >= zero) {
              temp = p5;
            }
            else {
              temp = p5*dirder/(dirder + p5*actred);
            }
            if (p1*fnorm1 >= fnorm || temp < p1) temp = p1;
            delta = temp*std::min(delta,pnorm/p1);
            par = par/temp;
            goto lbl_260;
          }
          lbl_240: {
            if (par != zero && ratio < p75) goto lbl_250;
            delta = pnorm/p5;
            par = p5*par;
            lbl_250:;
          }
          lbl_260:
          //
          // test for successful iteration.
          //
          if (ratio < p0001) goto lbl_290;
          //
          // successful iteration. update x, fvec, and their norms.
          //
          for(int j=1;j<=n;j++) {
            x(j) = wa2(j);
            wa2(j) = diag(j)*x(j);
          }
          for(int i=1;i<=m;i++) {
            fvec(i) = wa4(i);
          }
          xnorm = enorm(n,wa2);
          fnorm = fnorm1;
          iter = iter + 1;
          if (!call_back_after_iteration) goto lbl_290;
          iflag = 4;
          return;
         continue_after_call_back_after_iteration:
          lbl_290:
          //
          // tests for convergence.
          //
          if (   std::fabs(actred) <= ftol
              && prered <= ftol
              && p5*ratio <= one) info = 1;
          if (delta <= xtol*xnorm) info = 2;
          if (   std::fabs(actred) <= ftol
              && prered <= ftol
              && p5*ratio <= one
              && info == 2) info = 3;
          if (info != 0) goto lbl_300;
          //
          // tests for termination and stringent tolerances.
          //
          if (maxfev > 0 && nfev >= maxfev) info = 5;
          if (   std::fabs(actred) <= epsmch
              && prered <= epsmch
              && p5*ratio <= one) info = 6;
          if (delta <= epsmch*xnorm) info = 7;
          if (gnorm <= epsmch) info = 8;
          if (info != 0) goto lbl_300;
          //
          // end of the inner loop. repeat if iteration unsuccessful.
          //
          if (ratio < p0001) goto lbl_200;
          //
          // end of the outer loop.
          //
        goto lbl_30;
      lbl_300:
      //
      // termination, either normal or user imposed.
      //
      iflag = 0;
    }
  };

}}} // namespace scitbx::minpack::raw

#endif // SCITBX_MINPACK_LIB_H
