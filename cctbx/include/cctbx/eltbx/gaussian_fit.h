#ifndef CCTBX_ELTBX_GAUSSIAN_FIT_H
#define CCTBX_ELTBX_GAUSSIAN_FIT_H

#include <cctbx/eltbx/xray_scattering.h>
#include <scitbx/array_family/shared.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {

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
        bool enforce_positive_b) const
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
          if (!enforce_positive_b) {
            sh_b.push_back(b(j) + shifts[i]);
          }
          else {
            CCTBX_ASSERT(b(j) >= 0);
            double sqrt_b = std::sqrt(b(j));
            sh_b.push_back(scitbx::fn::pow2(sqrt_b + shifts[i]));
          }
        }
        if (i == shifts.size()) sh_c = 0;
        else                    sh_c = shifts[i];
        return gaussian_fit(
          stols_, target_values_, sigmas_, gaussian(sh_a,sh_b,sh_c));
      }

      double
      target_function(
        int power,
        bool use_sigmas,
        af::const_ref<double> const& differences)
      {
        CCTBX_ASSERT(differences.size() == stols_.size());
        CCTBX_ASSERT(power == 2 || power == 4);
        size_assert_intrinsic();
        af::const_ref<double> sigmas = sigmas_.const_ref();
        double sigma_squared = 1;
        double result = 0;
        for(std::size_t i=0;i<differences.size();i++) {
          double diff_squared = differences[i] * differences[i];
          double term = diff_squared;
          if (use_sigmas) {
            sigma_squared = sigmas[i] * sigmas[i];
            CCTBX_ASSERT(sigma_squared > 0);
            term /= sigma_squared;
          }
          if (power == 4) term *= diff_squared;
          result += term;
        }
        return result;
      }

      af::shared<double>
      gradients_w_r_t_abc(
        int power,
        bool use_sigmas,
        af::const_ref<double> const& differences,
        bool include_constant_term) const
      {
        CCTBX_ASSERT(differences.size() == stols_.size());
        CCTBX_ASSERT(power == 2 || power == 4);
        size_assert_intrinsic();
        af::shared<double> result(n_ab() * 2 + 1, 0);
        if (!include_constant_term) result.pop_back();
        af::const_ref<double> stols = stols_.const_ref();
        af::const_ref<double> sigmas = sigmas_.const_ref();
        double sigma_squared = 1;
        af::ref<double> g = result.ref();
        for(std::size_t i_point=0;i_point<stols.size();i_point++) {
          gaussian grg = gradients_at_stol(stols[i_point]);
          double diff = differences[i_point];
          double nwpd = 2 * diff;
          if (use_sigmas) {
            sigma_squared = sigmas[i_point] * sigmas[i_point];
            CCTBX_ASSERT(sigma_squared > 0);
            nwpd /= sigma_squared;
          }
          if (power == 4) nwpd *= 2 * diff * diff;
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

      af::shared<double>
      gradients_w_r_t_shifts(
        af::const_ref<double> const& shifts,
        af::const_ref<double> const& gradients_abc) const
      {
        CCTBX_ASSERT(shifts.size() >= n_ab() * 2);
        CCTBX_ASSERT(gradients_abc.size() == shifts.size());
        af::shared<double> result(af::adapt(gradients_abc));
        af::ref<double> g = result.ref();
        std::size_t i = n_ab();
        for(std::size_t j=0;j<n_ab();j++,i++) {
          CCTBX_ASSERT(b(j) >= 0);
          double sqrt_b = std::sqrt(b(j));
          g[i] *= 2 * (sqrt_b + shifts[i]);
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
