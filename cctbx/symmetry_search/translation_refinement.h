#ifndef CCTBX_SYMMETRY_SEARCH_TRANSLATION_REFINEMENT_H
#define CCTBX_SYMMETRY_SEARCH_TRANSLATION_REFINEMENT_H

#include <scitbx/vec3.h>

#include <scitbx/array_family/shared.h>
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
  
template <typename FloatType>
struct goodness_of_symmetry
{
  typedef FloatType real_type;
  typedef std::complex<FloatType> complex_type;
  typedef scitbx::vec3<real_type> vector_type;
  typedef scitbx::vec3<real_type> real_grad_type;
  typedef scitbx::vec3<complex_type> complex_grad_type;
  
  goodness_of_symmetry(sgtbx::space_group const &space_group,
                       af::const_ref<miller::index<> > const &indices,
                       af::const_ref<real_type> const &f_o,
                       miller::f_calc_map<real_type> &f_c,
                       vector_type const &x)
    : q(0), dq(0.)
  {
    // compute shifted and then symmetrised structure factors
    using namespace scitbx::constants;
    using namespace xray::structure_factors;
    SCITBX_ASSERT(f_o.size() == indices.size());
    scitbx::math::imaginary_unit_t i;
    int n = indices.size();
    
    multiplicities.reserve(n);
    y_x.reserve(n);
    dy_x.reserve(n);

    math::cos_sin_exact<real_type> exp_i_2pi;
    for (int k=0; k<indices.size(); ++k) {
      miller::index<> const &h = indices[k];
      multiplicities.push_back(space_group.multiplicity(h, f_c.anomalous_flag()));
      hr_ht_cache<real_type> hr_ht(space_group, h);
      complex_type f_x = 0;
      complex_grad_type df_x(0, 0, 0);
      for (int l=0; l<hr_ht.groups.size(); ++l) {
        hr_ht_group<real_type> const &g = hr_ht.groups[l];
        complex_type f  = f_c[g.hr]*exp_i_2pi(g.hr*x + g.ht);
        complex_grad_type df(g.hr);
        df *= i * two_pi * f;
        f_x  += f;
        df_x += df;
      }
      real_type abs_f_x = std::abs(f_x);
      if (abs_f_x == 0) continue;
      real_grad_type d_abs_f_x;
      for (int k=0; k<3; ++k) {
        d_abs_f_x[k] = f_x.real()*df_x[k].imag() + f_x.imag()*df_x[k].real();
      }
      d_abs_f_x /= 2*abs_f_x;
      y_x.push_back(abs_f_x);
      dy_x.push_back(d_abs_f_x);
    }
    
    // compute optimal lambda and mu
    scitbx::math::weighted_covariance<real_type> stats(
      y_x.ref(), // = x
      f_o, // = y 
      multiplicities.ref());
    lambda = stats.covariance_xy()/stats.variance_x();
    mu = stats.mean_y() - lambda*stats.mean_x();
    
    // compute the linearisation
    real_type c = *stats.correlation();
    q = stats.variance_y()*(1 - c*c);
    for (int k=0; k<indices.size(); ++k) {
      real_type r = lambda*y_x[k] + mu - f_o[k];
      dq += r*lambda*dy_x[k];
    }
  }
  
  af::shared<real_type> multiplicities;
  af::shared<real_type> y_x;
  af::shared<real_grad_type> dy_x;
  real_type lambda, mu;
  real_type q;
  vector_type dq;
};

}}



#endif

