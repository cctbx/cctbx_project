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
      target_at_d_star_sq(int power, double d_star_sq, double weight) const
      {
        double ptt = scitbx::fn::pow2(target_term_at_d_star_sq(d_star_sq));
        if (power == 2) return weight * ptt;
        if (power == 4) return weight * ptt * ptt;
        throw error("power must be 2 or 4.");
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
        int power,
        af::const_ref<double> const& d_star_sq,
        af::const_ref<double> const& weights,
        af::const_ref<double> const& target_terms,
        bool include_constant_term)
      {
        CCTBX_ASSERT(power == 2 || power == 4);
        CCTBX_ASSERT(weights.size() == d_star_sq.size());
        CCTBX_ASSERT(target_terms.size() == d_star_sq.size());
        af::shared<double> result(n_ab() * 2 + 1, 0);
        if (!include_constant_term) result.pop_back();
        af::ref<double> g = result.ref();
        for(std::size_t i_point=0;i_point<d_star_sq.size();i_point++) {
          gaussian grg = gradients_at_d_star_sq(d_star_sq[i_point]);
          double tt = target_terms[i_point];
          double nwptt = 2 * weights[i_point] * tt;
          if (power == 4) nwptt *= 2 * tt * tt;
          std::size_t i=0;
          for(;i<n_ab();i++) {
            g[i] += nwptt * grg.a(i);
          }
          for(std::size_t j=0;j<n_ab();j++,i++) {
            g[i] += nwptt * grg.b(j);
          }
          if (include_constant_term) {
            g[i] += nwptt * grg.c();
          }
        }
        return result;
      }

    protected:
      gaussian reference_gaussian_;
  };

  class gaussian_fit : public gaussian
  {
    public:
      gaussian_fit() {}

      gaussian_fit(
        af::shared<double> const& stols,
        af::shared<double> const& target_values,
        af::shared<double> const& sigmas,
        gaussian const& start_gaussian)
      :
        gaussian(start_gaussian),
        size_init_(stols.size()),
        stols_(stols),
        target_values_(target_values),
        sigmas_(sigmas)
      {
        CCTBX_ASSERT(target_values.size() == stols.size());
        CCTBX_ASSERT(sigmas.size() == stols.size() || sigmas.size() == 0);
        if (sigmas.size() == 0) sigmas_.resize(stols.size(), 1);
      }

      gaussian_fit(
        af::shared<double> const& stols,
        gaussian const& reference_gaussian,
        af::shared<double> const& sigmas,
        gaussian const& start_gaussian)
      :
        gaussian(start_gaussian),
        size_init_(stols.size()),
        stols_(stols),
        sigmas_(sigmas)
      {
        CCTBX_ASSERT(sigmas.size() == stols.size() || sigmas.size() == 0);
        target_values_.reserve(stols.size());
        if (sigmas.size() == 0) sigmas_.resize(stols.size(), 1);
        af::const_ref<double> s = stols_.const_ref();
        for(std::size_t i=0;i<s.size();i++) {
          target_values_.push_back(reference_gaussian.at_stol(s[i]));
        }
      }

      af::shared<double>
      stols() const { return stols_; }

      af::shared<double>
      target_values() const { return target_values_; }

      af::shared<double>
      sigmas() const { return sigmas_; }

      // Not available in Python.
      void
      size_assert_intrinsic() const
      {
        CCTBX_ASSERT(stols_.size() == size_init_);
        CCTBX_ASSERT(target_values_.size() == size_init_);
        CCTBX_ASSERT(sigmas_.size() == size_init_);
      }

      af::shared<double>
      fitted_values() const
      {
        size_assert_intrinsic();
        af::const_ref<double> stols = stols_.const_ref();
        af::shared<double> result((af::reserve(stols.size())));
        for(std::size_t i=0;i<stols.size();i++) {
          result.push_back(at_stol(stols[i]));
        }
        return result;
      }

      af::shared<double>
      differences() const
      {
        size_assert_intrinsic();
        af::const_ref<double> stols = stols_.const_ref();
        af::const_ref<double> target_values = target_values_.const_ref();
        af::shared<double> result((af::reserve(target_values.size())));
        for(std::size_t i=0;i<target_values.size();i++) {
          result.push_back(at_stol(stols[i]) - target_values[i]);
        }
        return result;
      }

      gaussian_fit
      apply_shifts(
        af::const_ref<double> const& shifts,
        double b_min) const
      {
        size_assert_intrinsic();
        CCTBX_ASSERT(   shifts.size() == n_ab() * 2
                     || shifts.size() == n_ab() * 2 + 1);
        af::small<float, gaussian::max_n_ab> sh_a;
        af::small<float, gaussian::max_n_ab> sh_b;
        float sh_c;
        std::size_t i=0;
        for(;i<n_ab();i++) {
          sh_a.push_back(a(i) + shifts[i]);
        }
        for(std::size_t j=0;j<n_ab();j++,i++) {
          sh_b.push_back(std::max(b_min, b(j) + shifts[i]));
        }
        if (i == shifts.size()) sh_c = 0;
        else                    sh_c = shifts[i];
        return gaussian_fit(
          stols_, target_values_, sigmas_, gaussian(sh_a,sh_b,sh_c));
      }

      af::shared<double>
      sum_of_gradients(
        int power,
        af::const_ref<double> const& differences,
        bool include_constant_term) const
      {
        CCTBX_ASSERT(differences.size() == stols_.size());
        size_assert_intrinsic();
        CCTBX_ASSERT(power == 2 || power == 4);
        af::shared<double> result(n_ab() * 2 + 1, 0);
        if (!include_constant_term) result.pop_back();
        af::const_ref<double> stols = stols_.const_ref();
        af::const_ref<double> sigmas = sigmas_.const_ref();
        af::ref<double> g = result.ref();
        for(std::size_t i_point=0;i_point<stols.size();i_point++) {
          gaussian grg = gradients_at_stol(stols[i_point]);
          double d = differences[i_point];
          double sigma_squared = sigmas[i_point] * sigmas[i_point];
          CCTBX_ASSERT(sigma_squared > 0);
          double nwpd = 2 * d / sigma_squared;
          if (power == 4) nwpd *= 2 * d * d;
          std::size_t i=0;
          for(;i<n_ab();i++) {
            g[i] += nwpd * grg.a(i);
          }
          for(std::size_t j=0;j<n_ab();j++,i++) {
            g[i] += nwpd * grg.b(j);
          }
          if (include_constant_term) {
            g[i] += nwpd * grg.c();
          }
        }
        return result;
      }

    protected:
      std::size_t size_init_;
      af::shared<double> stols_;
      af::shared<double> target_values_;
      af::shared<double> sigmas_;
  };

}}} // cctbx::eltbx::xray_scattering

#endif // CCTBX_ELTBX_GAUSSIAN_FIT_H
