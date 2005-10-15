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
  struct gill_murray_wright_decomposition_in_place
  {
    af::shared<double> e;
    af::shared<std::size_t> pivots;

    gill_murray_wright_decomposition_in_place() {}

    gill_murray_wright_decomposition_in_place(
      af::ref<double> const& u)
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
      unsigned jj = 0;
      for(unsigned j=0;j<n;j++) {
        // pivoting
        {
          unsigned jmax = j;
          unsigned ii = jj;
          double uii_max = fn::absolute(u[ii]);
          ii += (n-j);
          for(unsigned i=j+1;i<n;ii+=(n-i),i++) {
            double uii = fn::absolute(u[ii]);
            if (uii_max < uii) {
                uii_max = uii;
                jmax = i;
            }
          }
          if (jmax != j) {
            packed_u_swap_rows_and_columns_in_place(u, j, jmax);
            std::swap(pivots[j], pivots[jmax]);
          }
        }
        // compute j-th column of L^{T}
        {
          unsigned ii = 0;
          unsigned ij = j;
          for(unsigned i=0;i<j;ii+=(n-i),i++,ij+=(n-i)) {
            u[ij] /= u[ii];
          }
        }
        // update j-th row
        {
          unsigned kj = j;
          for(unsigned k=0;k<j;k++,kj+=(n-k)) {
            double f = u[kj];
            unsigned ki = kj;
            unsigned ji = jj;
            for(unsigned i=j;++i<n;) {
              u[++ji] -= u[++ki] * f;
            }
          }
        }
        // specify theta
        double theta_sq;
        if (j+1 == n) {
          theta_sq = 0;
        }
        else {
          unsigned jk = jj;
          theta_sq = fn::absolute(u[++jk]);
          for(unsigned k=j+2;k<n;k++) {
            double aajk = fn::absolute(u[++jk]);
            if (theta_sq < aajk) {
                theta_sq = aajk;
            }
          }
          theta_sq *= theta_sq;
        }
        // compute diagonal (j,j)
        {
          double f =
            std::max(delta, std::max(fn::absolute(u[jj]), theta_sq/beta_sq));
          e[j] = f - u[jj];
          u[jj] = f;
          // update remaining diagonals
          unsigned ji = jj;
          jj += (n-j);
          unsigned ii = jj;
          for(unsigned i=j+1;i<n;ii+=(n-i),i++) {
            u[ii] -= fn::pow2(u[++ji]) / f;
          }
        }
      }
      // scale
      ij = 0;
      for(unsigned i=0;i<n;i++) {
        double ujj = u[ij] = std::sqrt(u[ij]);
        ij++;
        for(unsigned j=i+1;j<n;j++) {
          u[ij++] *= ujj;
        }
      }
    }
  };

}}} // namespace scitbx::matrix::cholesky

#endif // SCITBX_MATRIX_CHOLESKY_H
