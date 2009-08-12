#ifndef SCITBX_MATH_GAUSSIAN_SUM_H
#define SCITBX_MATH_GAUSSIAN_SUM_H

#ifndef SCITBX_MATH_GAUSSIAN_SUM_MAX_N_TERMS
#define SCITBX_MATH_GAUSSIAN_SUM_MAX_N_TERMS 10
#endif

#include <scitbx/math/gaussian/term.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/error.h>

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
      explicit
      sum(FloatType const& c, bool use_c=true)
      :
        c_(c),
        use_c_(use_c)
      {
        SCITBX_ASSERT(use_c || c == 0);
      }

      //! Initialization of the terms and optionally the constant.
      /*! If c is different from zero use_c will automatically be
          set to true.
          <p>
          Not available in Python.
       */
      sum(
        af::small<term<FloatType>, max_n_terms> const& terms,
        FloatType const& c=0,
        bool use_c=false)
      :
        terms_(terms),
        c_(c),
        use_c_(use_c || c != 0)
      {}

      //! Initialization of the terms and optionally the constant.
      /*! If c is different from zero use_c will automatically be
          set to true.
       */
      sum(
        af::small<FloatType, max_n_terms> const& a,
        af::small<FloatType, max_n_terms> const& b,
        FloatType const& c=0,
        bool use_c=false)
      :
        c_(c),
        use_c_(use_c || c != 0)
      {
        SCITBX_ASSERT(a.size() == b.size());
        for(std::size_t i=0;i<a.size();i++) {
          terms_.push_back(term<FloatType>(a[i], b[i]));
        }
      }

      //! Initialization of the terms and optionally the constant.
      /*! If c is different from zero use_c will automatically be
          set to true. If ab contains an odd number of elements
          the last element is used to initialize c and use_c will
          bet set to true.
       */
      sum(
        af::const_ref<FloatType> const& ab,
        FloatType const& c=0,
        bool use_c=false)
      :
        c_(c),
        use_c_(use_c || c != 0)
      {
        SCITBX_ASSERT(!use_c || ab.size() % 2 == 0);
        SCITBX_ASSERT(ab.size() / 2 <= max_n_terms);
        std::size_t n = ab.size();
        if (n % 2 != 0) {
          n--;
          c_ = ab.back();
          use_c_ = true;
        }
        for(std::size_t i=0;i<n;i+=2) {
          terms_.push_back(term<FloatType>(ab[i], ab[i+1]));
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
      af::small<FloatType, max_n_terms>
      array_of_a() const
      {
        af::small<FloatType, max_n_terms> result;
        for(std::size_t i=0;i<terms_.size();i++) {
          result.push_back(terms_[i].a);
        }
        return result;
      }

      //! Array of coefficients b.
      af::small<FloatType, max_n_terms>
      array_of_b() const
      {
        af::small<FloatType, max_n_terms> result;
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

      //! Total number of parameters.
      std::size_t
      n_parameters() const
      {
        std::size_t result = n_terms() * 2;
        if (use_c()) result++;
        return result;
      }

      //! Array of parameters a0,b0,...,ai,bi[,c].
      af::shared<FloatType>
      parameters() const
      {
        af::shared<FloatType> result;
        result.reserve(n_parameters());
        for(std::size_t i=0;i<terms_.size();i++) {
          term<FloatType> const& ti = terms_[i];
          result.push_back(ti.a);
          result.push_back(ti.b);
        }
        if (use_c()) result.push_back(c_);
        return result;
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

      //! Sum of Gaussian terms at the points x, given x^2.
      af::shared<FloatType>
      at_x_sq(af::const_ref<FloatType> const& x_sq) const
      {
        af::shared<double> result(x_sq.size(),af::init_functor_null<double>());
        for(std::size_t i=0;i<x_sq.size();i++) {
          result[i] = at_x_sq(x_sq[i]);
        }
        return result;
      }

      //! Sum of Gaussian terms at the point x.
      FloatType
      at_x(FloatType const& x) const { return at_x_sq(x * x); }

      //! Sum of Gaussian terms at the points x.
      af::shared<FloatType>
      at_x(af::const_ref<FloatType> const& x) const
      {
        af::shared<double> result(x.size(),af::init_functor_null<double>());
        for(std::size_t i=0;i<x.size();i++) {
          result[i] = at_x(x[i]);
        }
        return result;
      }

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
        FloatType const& b_min_for_erf_based_algorithm=1e-3)
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
