#ifndef CCTBX_ELTBX_GAUSSIAN_FIT_H
#define CCTBX_ELTBX_GAUSSIAN_FIT_H

#include <cctbx/eltbx/xray_scattering.h>
#include <scitbx/array_family/shared.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {

  class difference_gaussian : public gaussian
  {
    public:
      difference_gaussian() {}

      difference_gaussian(
        gaussian const& reference_gaussian,
        gaussian const& target_gaussian)
      :
        gaussian(target_gaussian),
        reference_gaussian_(reference_gaussian)
      {
        CCTBX_ASSERT(reference_gaussian.n_ab() >= target_gaussian.n_ab());
      }

      gaussian const&
      reference_gaussian() const { return reference_gaussian_; }

      difference_gaussian
      apply_shifts(
        af::const_ref<double> const& shifts,
        double b_min)
      {
        CCTBX_ASSERT(   shifts.size() == n_ab() * 2
                     || shifts.size() == n_ab() * 2 + 1);
        difference_gaussian result;
        std::size_t i=0;
        for(;i<n_ab();i++) {
          result.a_.push_back(a(i) + shifts[i]);
        }
        for(std::size_t j=0;j<n_ab();j++,i++) {
          result.b_.push_back(std::max(b_min, b(j) + shifts[i]));
        }
        if (i == shifts.size()) result.c_ = 0;
        else                    result.c_ = shifts[i];
        result.reference_gaussian_ = reference_gaussian_;
        return result;
      }

      double
      target_term_at_d_star_sq(double d_star_sq) const
      {
        return at_d_star_sq(d_star_sq)
             - reference_gaussian_.at_d_star_sq(d_star_sq);
      }

      double
      target_at_d_star_sq(double d_star_sq, double weight) const
      {
        return weight * scitbx::fn::pow2(target_term_at_d_star_sq(d_star_sq));
      }

      af::shared<double>
      target_terms_at_points(
        af::const_ref<double> const& d_star_sq)
      {
        af::shared<double> result((af::reserve(d_star_sq.size())));
        for(std::size_t i_point=0;i_point<d_star_sq.size();i_point++) {
          result.push_back(target_term_at_d_star_sq(d_star_sq[i_point]));
        }
        return result;
      }

      af::shared<double>
      sum_of_gradients_at_points(
        af::const_ref<double> const& d_star_sq,
        af::const_ref<double> const& weights,
        af::const_ref<double> const& target_terms,
        bool include_constant_term)
      {
        CCTBX_ASSERT(weights.size() == d_star_sq.size());
        CCTBX_ASSERT(target_terms.size() == d_star_sq.size());
        af::shared<double> result(n_ab() * 2 + 1, 0);
        if (!include_constant_term) result.pop_back();
        af::ref<double> g = result.ref();
        for(std::size_t i_point=0;i_point<d_star_sq.size();i_point++) {
          gaussian grg = gradients_at_d_star_sq(d_star_sq[i_point]);
          double twtt = 2 * weights[i_point] * target_terms[i_point];
          std::size_t i=0;
          for(;i<n_ab();i++) {
            g[i] += twtt * grg.a(i);
          }
          for(std::size_t j=0;j<n_ab();j++,i++) {
            g[i] += twtt * grg.b(j);
          }
          if (include_constant_term) {
            g[i] += twtt * grg.c();
          }
        }
        return result;
      }

    protected:
      gaussian reference_gaussian_;
  };

}}} // cctbx::eltbx::xray_scattering

#endif // CCTBX_ELTBX_GAUSSIAN_FIT_H
