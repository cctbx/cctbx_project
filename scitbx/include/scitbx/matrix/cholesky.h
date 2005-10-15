#ifndef SCITBX_MATRIX_CHOLESKY_H
#define SCITBX_MATRIX_CHOLESKY_H

#include <scitbx/matrix/move.h>
#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace matrix { namespace cholesky {

  template <typename FloatType>
  af::shared<FloatType>
  decomposition(
    af::const_ref<FloatType> const& a,
    FloatType const& relative_eps=1.e-15)
  {
    unsigned n = symmetric_n_from_packed_size(a.size());
    FloatType eps;
    if (relative_eps <= 0 || a.size() == 0) {
      eps = 0;
    }
    else {
      eps = relative_eps * af::max_absolute(a);
    }
    af::shared<FloatType> result(a.size(), af::init_functor_null<FloatType>());
    FloatType *c = result.begin();
    std::size_t i_a = 0;
    std::size_t k0 = 0;
    for(unsigned k=0;k<n;k++,k0+=k) {
      std::size_t kj = k0;
      FloatType sum = 0;
      for(unsigned j=0;j<k;j++) {
        FloatType const& c_kj = c[kj++];
        sum += c_kj * c_kj;
      }
      FloatType d = a[i_a++] - sum;
      if (d <= eps) return af::shared<FloatType>();
      d = std::sqrt(d);
      c[kj++] = d;
      std::size_t i0 = kj;
      for(unsigned i=k+1;i<n;i++,i0+=i) {
        std::size_t ij = i0;
        kj = k0;
        sum = 0;
        for(unsigned j=0;j<k;j++) {
          sum += c[ij++] * c[kj++];
        }
        c[ij] = (a[i_a++] - sum) / d;
      }
    }
    return result;
  }

  //! Computes P^{T}AP + E = LDL^{T}
  /*! P. E. Gill, W. Murray, and M. H. Wright:
      Practical Optimization.
      New York: Academic Press, 1981.

      See also:

        J. Nocedal and S. Wright:
        Numerical Optimization.
        Springer, New York, 1999, pp. 145-150.

        http://lib.stat.cmu.edu/S/glm
        file: dmcdc.r
   */
  struct gill_murray_wright
  {
    af::shared<double> e;
    af::shared<std::size_t> pivots;

    gill_murray_wright() {}

    gill_murray_wright(af::ref<double> const& u)
    {
      unsigned n = symmetric_n_from_packed_size(u.size());
      e.resize(n);
      pivots.reserve(n);
      for(unsigned i=0;i<n;i++) {
        pivots.push_back(i);
      }
      double fp_eps = math::floating_point_epsilon<double>::get();
      double gamma = 0;
      double chi = 0;
      unsigned ij = 0;
      for(unsigned i=0;i<n;i++) {
        gamma = std::max(gamma, fn::absolute(u[ij++]));
        for(unsigned j=i+1;j<n;j++) {
          chi = std::max(chi, fn::absolute(u[ij++]));
        }
      }
      double delta = fp_eps * std::max(gamma + chi, static_cast<double>(1));
      double beta_sq = std::max(gamma, fp_eps);
      if (n > 1) {
        beta_sq = std::max(beta_sq, chi/std::sqrt(static_cast<double>(n*n-1)));
      }
      af::ref<double, packed_u_accessor> a(u.begin(), packed_u_accessor(n));
      for(unsigned j=0;j<n;j++) {
        unsigned jmax = j;
        for(unsigned i=j+1;i<n;i++) {
          if (fn::absolute(a(jmax,jmax)) < fn::absolute(a(i,i))) {
            jmax = i;
          }
        }
        // pivoting
        if (jmax != j) {
          packed_u_swap_rows_and_columns_in_place(u, j, jmax);
          std::swap(pivots[j], pivots[jmax]);
        }
        // compute j-th column of L^{T}
        for(unsigned i=0;i<j;i++) {
          a(i,j) /= a(i,i);
        }
        // update j-th row
        for(unsigned k=0;k<j;k++) {
          for(unsigned i=j+1;i<n;i++) {
            a(j,i) -= a(k,i) * a(k,j);
          }
        }
        // specify theta
        jmax = j+1;
        double theta_sq;
        if (jmax == n) {
          theta_sq = 0;
        }
        else {
          for(unsigned k=jmax+1;k<n;k++) {
            if (fn::absolute(a(j,jmax)) < fn::absolute(a(j,k))) {
              jmax = k;
            }
          }
          theta_sq = fn::pow2(a(j,jmax));
        }
        // compute diagonal (j,j)
        double tmp =
          std::max(delta, std::max(fn::absolute(a(j,j)), theta_sq/beta_sq));
        e[j] = tmp - a(j,j);
        a(j,j) = tmp;
        // update remaining diagonals
        for(unsigned i=j+1;i<n;i++) {
          a(i,i) -= fn::pow2(a(j,i)) / a(j,j);
        }
      }
      // scale
      ij = 0;
      for(unsigned i=0;i<n;i++) {
        double ajj = a[ij] = std::sqrt(a[ij]);
        ij++;
        for(unsigned j=i+1;j<n;j++) {
          a[ij++] *= ajj;
        }
      }
    }
  };

}}} // namespace scitbx::matrix::cholesky

#endif // SCITBX_MATRIX_CHOLESKY_H
