#ifndef CCTBX_SYMMETRY_SEARCH_TRANSLATION_REFINEMENT_H
#define CCTBX_SYMMETRY_SEARCH_TRANSLATION_REFINEMENT_H

#include <scitbx/vec3.h>

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/math/imaginary.h>
#include <scitbx/math/weighted_covariance.h>
#include <cctbx/error.h>
#include <cctbx/math/cos_sin_table.h>

#include <cctbx/miller/f_calc_map.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/xray/hr_ht_cache.h>
#include <cmath>
#include <complex>

namespace cctbx { namespace symmetry_search {

/** Average over space group symmetries of shifted structure factors.
 This average reads
   \f[
      f_x(h) = \sum_{(R|t) \in G} f_c(hR) e^{i2\pi hRx} e^{i2\pi ht},
   \f]
 where \f$\{f_c(h)\}\f$ have an assumed P1 symmetry.
 Centring translations are actually ignored since they result in an
 overall scale only.
 */
template <typename FloatType>
struct symmetrised_shifted_structure_factors
{
  typedef FloatType real_type;
  typedef std::complex<FloatType> complex_type;
  typedef scitbx::vec3<real_type> vector_type;
  typedef af::tiny<complex_type, 3> complex_grad_type;

  /// The array of f_x(h) for the Miller indices h specified to the constructor.
  af::shared<complex_type> f_x;

  /// The array of grad f_x(h) for the Miller indices h passed to the constructor.
  af::shared<complex_grad_type> grad_f_x;

  /// Construct f_x and its gradients if specified so.
  symmetrised_shifted_structure_factors(sgtbx::space_group const &space_group,
                                        af::const_ref<miller::index<> > const
                                        &indices,
                                        miller::f_calc_map<real_type> &f_c,
                                        vector_type const &x,
                                        bool compute_gradient)
  {
    // compute shifted and then symmetrised structure factors
    using namespace scitbx::constants;
    using namespace xray::structure_factors;
    scitbx::math::imaginary_unit_t i;

    int n = indices.size();

    f_x.reserve(n);

    math::cos_sin_exact<real_type> exp_i_2pi;
    for (int k=0; k<indices.size(); ++k) {
      miller::index<> const &h = indices[k];
      hr_ht_cache<real_type> hr_ht(space_group, h);
      complex_type f_x = 0;
      complex_grad_type grad_f_x(0);
      for (int l=0; l<hr_ht.groups.size(); ++l) {
        hr_ht_group<real_type> const &g = hr_ht.groups[l];
        complex_type f  = f_c[g.hr]*exp_i_2pi(g.hr*x + g.ht);
        f_x  += f;
        if (!compute_gradient) continue;
        grad_f_x += (i*two_pi*f) * complex_grad_type(g.hr);
      }
      if (hr_ht.is_centric) f_x += std::conj(f_x)*hr_ht.f_h_inv_t;
      this->f_x.push_back(f_x);
      if (compute_gradient) {
        if (hr_ht.is_centric) grad_f_x += af::conj(grad_f_x)*hr_ht.f_h_inv_t;
        this->grad_f_x.push_back(grad_f_x);
      }
    }
  }
};


/// L.S. misfit on F^2 between observed and shifted symmetrised s.f.
/**
  It reads
  \f[
     q(x, \lambda, \mu)
     = \sum_h w(h) \left(\lambda |f_x(h)|^2 + \mu - f_o(h)^2\right)^2
     \f]
  where w(h) is typically the multiplicity of h.

  For the given fixed x, \f$q(x, \lambda, \mu)\f$ is minimised for
    varying \f$\lambda\f$ and \f$\mu\f$.
*/

template <typename FloatType>
struct ls_with_scale_and_bias
{
  typedef symmetrised_shifted_structure_factors<FloatType> sssf_t;
  typedef typename sssf_t::real_type real_type;
  typedef typename sssf_t::complex_type complex_type;
  typedef typename sssf_t::vector_type vector_type;
  typedef typename sssf_t::complex_grad_type complex_grad_type;
  typedef scitbx::vec3<real_type> real_grad_type;

  /// Least-squares overall scale and bias minimising \f$q\f$
  /** They are function of the \f$x\f$ passed to the contructor */
  real_type lambda, mu;

  /// Minimum of \f$q\f$ when \f$\lambda\f$ and \f$\mu\f$ vary.
  real_type q;

  /// Correlation between \f$\{|f_x(h)\}\f$ and \f$\{f_o(h)^2\}\f$
  real_type c;

  /// Gradient of \f$q(x, \lambda(x), \mu(x))\f$ wrt \f$x\f$.
  vector_type grad_q;

  ls_with_scale_and_bias(af::const_ref<complex_type> const &f_x,
                         af::const_ref<complex_grad_type> const &grad_f_x,
                         af::const_ref<real_type> const &f_o_sq,
                         af::const_ref<real_type> const &weight)
  : q(0), grad_q(0.)
  {
    CCTBX_ASSERT(f_x.size() == weight.size());
    CCTBX_ASSERT(f_o_sq.size() == weight.size());
    CCTBX_ASSERT(!grad_f_x.size() || grad_f_x.size() == weight.size());
    bool compute_gradient = grad_f_x.size();
    int n = weight.size();

    // array of |f_x(h)|^2 and grad_x |f_x(h)|^2
    af::shared<real_type> f_x_modulus_sq_;
    f_x_modulus_sq_.reserve(n);
    af::shared<real_grad_type> df_x_modulus_sq_;
    if (compute_gradient) df_x_modulus_sq_.reserve(n);
    for (int k=0; k<n; ++k) {
      f_x_modulus_sq_.push_back(std::norm(f_x[k]));
      if (compute_gradient) {
        real_grad_type d_norm_f_x;
        for (int l=0; l<3; ++l) {
          d_norm_f_x[l] = 2*(  f_x[k].real()*grad_f_x[k][l].real()
                             + f_x[k].imag()*grad_f_x[k][l].imag());
        }
        df_x_modulus_sq_.push_back(d_norm_f_x);
      }
    }

    // compute optimal lambda and mu
    af::const_ref<real_type> f_x_modulus_sq = f_x_modulus_sq_.const_ref();
    scitbx::math::weighted_covariance<real_type> stats(f_x_modulus_sq, // = x
                                                       f_o_sq, // = y
                                                       weight);
    lambda = stats.covariance_xy()/stats.variance_x();
    mu = stats.mean_y() - lambda*stats.mean_x();

    // compute the linearisation of q
    c = *stats.correlation();
    q = stats.variance_y()*(1 - c*c);

    if (compute_gradient) {
      af::const_ref<real_grad_type> df_x_modulus_sq = df_x_modulus_sq_.const_ref();
      for (int k=0; k<n; ++k) {
        real_type r = lambda*f_x_modulus_sq[k] + mu - f_o_sq[k];
        grad_q += weight[k]*2*r*lambda*df_x_modulus_sq[k];
      }
      grad_q /= stats.sum_weights();
    }

  }
};


}}



#endif
