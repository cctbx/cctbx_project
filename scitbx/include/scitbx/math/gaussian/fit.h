#ifndef SCITBX_MATH_GAUSSIAN_FIT_H
#define SCITBX_MATH_GAUSSIAN_FIT_H

#include <scitbx/math/gaussian/sum.h>

namespace scitbx { namespace math { namespace gaussian {

  template <typename FloatType=double>
  class fit : public sum<FloatType>
  {
    public:
      typedef sum<FloatType> base_t;

      fit() {}

      fit(
        af::shared<FloatType> const& table_x,
        af::shared<FloatType> const& table_y,
        af::shared<FloatType> const& table_sigmas,
        base_t const& start)
      :
        base_t(start),
        size_init_(table_x.size()),
        table_x_(table_x),
        table_y_(table_y),
        table_sigmas_(table_sigmas)
      {
        SCITBX_ASSERT(table_y.size() == table_x.size());
        SCITBX_ASSERT(table_sigmas.size() == table_x.size());
      }

      fit(
        af::shared<FloatType> const& table_x,
        base_t const& reference,
        af::shared<FloatType> const& table_sigmas,
        base_t const& start)
      :
        base_t(start),
        size_init_(table_x.size()),
        table_x_(table_x),
        table_sigmas_(table_sigmas)
      {
        SCITBX_ASSERT(table_sigmas.size() == table_x.size());
        af::const_ref<double> x = table_x_.const_ref();
        table_y_.reserve(x.size());
        for(std::size_t i=0;i<x.size();i++) {
          table_y_.push_back(reference.at_x(x[i]));
        }
      }

      af::shared<FloatType>
      table_x() const { return table_x_; }

      af::shared<FloatType>
      table_y() const { return table_y_; }

      af::shared<FloatType>
      table_sigmas() const { return table_sigmas_; }

      // Not available in Python.
      void
      size_assert_intrinsic() const
      {
        SCITBX_ASSERT(table_x_.size() == size_init_);
        SCITBX_ASSERT(table_y_.size() == size_init_);
        SCITBX_ASSERT(table_sigmas_.size() == size_init_);
      }

      af::shared<FloatType>
      fitted_values() const
      {
        size_assert_intrinsic();
        af::const_ref<FloatType> x = table_x_.const_ref();
        af::shared<FloatType> result((af::reserve(x.size())));
        for(std::size_t i=0;i<x.size();i++) {
          result.push_back(this->at_x(x[i]));
        }
        return result;
      }

      af::shared<FloatType>
      differences() const
      {
        size_assert_intrinsic();
        af::const_ref<FloatType> x = table_x_.const_ref();
        af::const_ref<FloatType> y = table_y_.const_ref();
        af::shared<FloatType> result((af::reserve(y.size())));
        for(std::size_t i=0;i<y.size();i++) {
          result.push_back(this->at_x(x[i]) - y[i]);
        }
        return result;
      }

      af::shared<FloatType>
      significant_relative_errors() const
      {
        using scitbx::fn::absolute;
        af::shared<FloatType> diffs_ = differences();
        af::const_ref<FloatType> diffs = diffs_.const_ref();
        af::const_ref<FloatType> y = table_y_.const_ref();
        af::const_ref<FloatType> sigmas = table_sigmas_.const_ref();
        FloatType zero(0);
        af::shared<FloatType> results(af::reserve(diffs.size()));
        for(std::size_t i=0;i<diffs.size();i++) {
          FloatType result = std::max(zero, absolute(diffs[i]) - sigmas[i]);
          if (result > 0) {
            SCITBX_ASSERT(absolute(y[i]) > 0 || sigmas[i] > 0);
            result /= std::max(sigmas[i], absolute(y[i]));
          }
          results.push_back(result);
        }
        return results;
      }

      fit
      apply_shifts(
        af::const_ref<FloatType> const& shifts,
        bool enforce_positive_b) const
      {
        size_assert_intrinsic();
        SCITBX_ASSERT(shifts.size() == this->n_parameters());
        af::small<term<FloatType>, base_t::max_n_terms> sh_terms;
        FloatType sh_c;
        std::size_t j=0;
        for(std::size_t i=0;i<this->n_terms();i++) {
          FloatType sh_a = this->terms_[i].a + shifts[j++];
          FloatType sh_b;
          if (!enforce_positive_b) {
            sh_b = this->terms_[i].b + shifts[j++];
          }
          else {
            FloatType b = this->terms_[i].b;
            SCITBX_ASSERT(b >= 0);
            sh_b = fn::pow2(std::sqrt(b) + shifts[j++]);
          }
          sh_terms.push_back(term<FloatType>(sh_a, sh_b));
        }
        if (this->use_c()) sh_c = this->c_ + shifts[j];
        else               sh_c = 0;
        return fit(this->table_x_, this->table_y_, this->table_sigmas_,
                   base_t(sh_terms, sh_c, this->use_c()));
      }

      FloatType
      target_function(
        int power,
        bool use_sigmas,
        af::const_ref<FloatType> const& differences)
      {
        SCITBX_ASSERT(differences.size() == table_x_.size());
        SCITBX_ASSERT(power == 2 || power == 4);
        size_assert_intrinsic();
        af::const_ref<FloatType> sigmas = table_sigmas_.const_ref();
        FloatType sigma_squared = 1;
        FloatType result = 0;
        for(std::size_t i=0;i<differences.size();i++) {
          FloatType diff_squared = differences[i] * differences[i];
          FloatType term = diff_squared;
          if (use_sigmas) {
            sigma_squared = sigmas[i] * sigmas[i];
            SCITBX_ASSERT(sigma_squared > 0);
            term /= sigma_squared;
          }
          if (power == 4) term *= diff_squared;
          result += term;
        }
        return result;
      }

      af::shared<FloatType>
      gradients_d_abc(
        int power,
        bool use_sigmas,
        af::const_ref<FloatType> const& differences)
      {
        SCITBX_ASSERT(differences.size() == this->table_x_.size());
        SCITBX_ASSERT(power == 2 || power == 4);
        size_assert_intrinsic();
        af::shared<FloatType> result(this->n_parameters());
        af::const_ref<FloatType> x = this->table_x_.const_ref();
        af::const_ref<FloatType> sigmas = this->table_sigmas_.const_ref();
        FloatType sigma_squared = 1;
        af::ref<FloatType> g = result.ref();
        for(std::size_t i_point=0;i_point<x.size();i_point++) {
          FloatType x_sq = fn::pow2(x[i_point]);
          FloatType diff = differences[i_point];
          FloatType nwpd = 2 * diff;
          if (use_sigmas) {
            sigma_squared = sigmas[i_point] * sigmas[i_point];
            SCITBX_ASSERT(sigma_squared > 0);
            nwpd /= sigma_squared;
          }
          if (power == 4) nwpd *= 2 * diff * diff;
          std::size_t j=0;
          for(std::size_t i=0;i<this->n_terms();i++) {
            term<FloatType>
              d_ab = this->terms_[i].gradients_d_ab_at_x_sq(x_sq);
            g[j++] += nwpd * d_ab.a;
            g[j++] += nwpd * d_ab.b;
          }
          if (this->use_c()) g[j] += nwpd;
        }
        return result;
      }

      af::shared<FloatType>
      gradients_d_shifts(
        af::const_ref<FloatType> const& shifts,
        af::const_ref<FloatType> const& gradients_d_abc) const
      {
        SCITBX_ASSERT(shifts.size() == this->n_parameters());
        SCITBX_ASSERT(gradients_d_abc.size() == shifts.size());
        af::shared<FloatType> result(af::adapt(gradients_d_abc));
        af::ref<FloatType> g = result.ref();
        std::size_t j=1;
        for(std::size_t i=0;i<this->n_terms();i++,j+=2) {
          FloatType b = this->terms_[i].b;
          SCITBX_ASSERT(b >= 0);
          g[j] *= 2 * (std::sqrt(b) + shifts[j]);
        }
        return result;
      }

    protected:
      std::size_t size_init_;
      af::shared<FloatType> table_x_;
      af::shared<FloatType> table_y_;
      af::shared<FloatType> table_sigmas_;
  };

}}} // scitbx::math::gaussian

#endif // SCITBX_MATH_GAUSSIAN_FIT_H
