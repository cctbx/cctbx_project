#ifndef SCITBX_MATH_GAUSSIAN_SUM_H
#define SCITBX_MATH_GAUSSIAN_SUM_H

#ifndef SCITBX_MATH_GAUSSIAN_SUM_MAX_N_TERMS
#define SCITBX_MATH_GAUSSIAN_SUM_MAX_N_TERMS 10
#endif

#include <scitbx/math/gaussian/term.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace math { namespace gaussian {

  //! Sum of Gaussian terms plus a constant.
  template <typename FloatType=double>
  class sum
  {
    public:
      //! Maximum number of terms.
      BOOST_STATIC_CONSTANT(std::size_t,
        max_n_terms=SCITBX_MATH_GAUSSIAN_SUM_MAX_N_TERMS);

      //! Default constructor. Some data members are not initialized!
      sum() {}

      //! Initialization of the constant.
      sum(FloatType const& c)
      :
        c_(c),
        use_c_(true)
      {}

      //! Initialization of the terms and optionally the constant.
      /*! If c is different from zero use_c will automatically be
          set to true.
       */
      sum(
        af::const_ref<FloatType> const& a,
        af::const_ref<FloatType> const& b,
        FloatType const& c=0,
        bool use_c=false)
      :
        c_(c),
        use_c_(use_c)
      {
        SCITBX_ASSERT(a.size() == b.size());
        SCITBX_ASSERT(a.size() <= max_n_terms);
        for(std::size_t i=0;i<a.size();i++) {
          terms_.push_back(term<FloatType>(a[i], b[i]));
        }
      }

      //! Number of terms.
      std::size_t
      n_terms() const { return terms_.size(); }

      //! Terms.
      /*! Not available in Python.
       */
      af::small<term<FloatType>, max_n_terms> const&
      terms() const { return terms_; }

      //! Array of coefficients a.
      af::shared<FloatType>
      array_of_a() const
      {
        af::shared<FloatType> result(af::reserve(terms_.size()));
        for(std::size_t i=0;i<terms_.size();i++) {
          result.push_back(terms_[i].a);
        }
        return result;
      }

      //! Array of coefficients b.
      af::shared<FloatType>
      array_of_b() const
      {
        af::shared<FloatType> result(af::reserve(terms_.size()));
        for(std::size_t i=0;i<terms_.size();i++) {
          result.push_back(terms_[i].b);
        }
        return result;
      }

      //! Coefficient c.
      FloatType const&
      c() const { return c_; }

      //! Flag.
      bool
      use_c() const { return use_c_; }

      //! Test if n_terms() == 0 and c() == 0.
      bool
      all_zero() const
      {
        return n_terms() == 0 && c() == 0;
      }

      //! Sum of Gaussian terms at the point x, given x^2.
      FloatType
      at_x_sq(FloatType const& x_sq) const
      {
        FloatType result = c_;
        for(std::size_t i=0;i<terms_.size();i++) {
          result += terms_[i].at_x_sq(x_sq);
        }
        return result;
      }

      //! Sum of Gaussian terms at the point x.
      FloatType
      at_x(FloatType const& x) const { return at_x_sq(x * x); }

      //! Gradient w.r.t. x at the point x.
      FloatType
      gradient_dx_at_x(FloatType const& x) const
      {
        FloatType result = 0;
        for(std::size_t i=0;i<terms_.size();i++) {
          result += terms_[i].gradient_dx_at_x(x);
        }
        return result;
      }

      //! Integral dx from 0 to the point x.
      FloatType
      integral_dx_at_x(
        FloatType const& x,
        FloatType const& b_min_for_erf_based_algorithm=1.e-3)
      {
        FloatType result = c_ * x;
        for(std::size_t i=0;i<terms_.size();i++) {
          result += terms_[i].integral_dx_at_x(
            x, b_min_for_erf_based_algorithm);
        }
        return result;
      }

    protected:
      af::small<term<FloatType>, max_n_terms> terms_;
      FloatType c_;
      bool use_c_;
  };

}}} // scitbx::math::gaussian

#endif // SCITBX_MATH_GAUSSIAN_SUM_H
