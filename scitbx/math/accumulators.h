#ifndef SCITBX_MATH_ACCUMULATORS
#define SCITBX_MATH_ACCUMULATORS

#include <cmath>
#include <cstddef>
#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>

namespace scitbx { namespace math {

///Accumulators to compute partial series. Statistics are the main application.
/**
The classes in this namespace are designed to help separating the computation
of statistics from the iteration over the sequence of values.

See class \link basic_statistics \endlink for a typical example of their use.

*/
namespace accumulator {

template<typename DataType>
class null_accumulator
{
  public:
    null_accumulator() {}
    null_accumulator(DataType x0) {}
    null_accumulator(DataType x0, DataType x1) {}
    void operator()(DataType x) {}
    void reset() {}
};


template<typename DataType>
class enumerated_accumulator
{
  public:
    enumerated_accumulator()
      : n(0)
    {}
    enumerated_accumulator(DataType x0)
      : n(1)
    {}

    void operator()(DataType x) { n++; }

    std::size_t count() const { return n; }

    void reset() { n=0; }

  private:
    std::size_t n;

};


template<typename FloatType, class Previous=null_accumulator<FloatType> >
class min_max_accumulator : public Previous
{
  public:
    min_max_accumulator(FloatType x0)
      : Previous(x0),
        min_(x0), max_(x0), max_absolute_(std::abs(x0))
    {}

    void operator()(FloatType x) {
      Previous::operator()(x);
      if (x < min_) min_ = x;
      if (max_ < x) max_ = x;
    }

    FloatType min() const { return min_; }
    FloatType max() const { return max_; }
    FloatType max_absolute() const { return -min_ > max_ ? -min_ : max_; }

    void reset() {
      Previous::reset();
      min_ = max_ = max_absolute_ = 0;
    }

  private:
    FloatType min_, max_, max_absolute_;
};


template<typename FloatType, class Previous=null_accumulator<FloatType> >
class mean_variance_accumulator : public Previous
{
  public:
    mean_variance_accumulator()
      : Previous(),
        mean_(0), m2(0)
    {}

    mean_variance_accumulator(FloatType x0)
      : Previous(x0),
        mean_(x0), m2(0)
    {}

    void operator()(FloatType x) {
      Previous::operator()(x);
      std::size_t n = this->count();
      FloatType delta = x - mean_;
      mean_ += delta/n;
      m2 += delta*(x - mean_);
    }

    FloatType mean() const { return mean_; }

    FloatType sum() const { return this->count()*mean_; }

    FloatType biased_variance() const {
      std::size_t n = this->count();
      return m2/n;
    }

    FloatType mean_squares() const { return biased_variance() + mean()*mean(); }

    FloatType unbiased_variance() const {
      std::size_t n = this->count();
      return biased_variance()/(1 - 1./n);
    }

    FloatType biased_standard_deviation() const {
      return std::sqrt(biased_variance());
    }

    FloatType unbiased_standard_deviation() const {
      return std::sqrt(unbiased_variance());
    }

    void reset() {
      Previous::reset();
      mean_ = m2 = 0;
    }

  private:
    FloatType mean_, m2;
};


template<typename FloatType>
class deviation_accumulator
{
  public:
    deviation_accumulator(FloatType about)
      : about_(about), n(0)
    {}

    void operator()(FloatType x) {
      n++;
      delta = x - about_;
    }

    std::size_t count() const { return n; }
    FloatType about() const { return about_; }
    FloatType deviation() const { return delta; }

    void reset() {
      n = 0;
    }

  private:
    FloatType about_, width_, delta;
    std::size_t n;
};

template<typename FloatType>
class normalised_deviation_accumulator : public deviation_accumulator<FloatType>
{
  public:
    normalised_deviation_accumulator(FloatType about, FloatType width)
      : deviation_accumulator<FloatType>(about), width_(width)
    {}

    void operator()(FloatType x) {
      deviation_accumulator<FloatType>::operator()(x);
      u = this->deviation()/width_;
    }

    FloatType width() const { return width_; }
    FloatType normalised_deviation() const { return u; }

  private:
    FloatType width_, u;
    std::size_t n;
};


template<typename FloatType, class Previous=null_accumulator<FloatType> >
class skewness_accumulator : public Previous
{
  public:
    skewness_accumulator(FloatType mean, FloatType standard_deviation)
      : Previous(mean, standard_deviation),
        skewness_(0)
    {}

    void operator()(FloatType x) {
      Previous::operator()(x);
      std::size_t n = this->count();
      FloatType u = this->normalised_deviation();
      skewness_ += (u*u*u - skewness_)/n;
    }

    FloatType skewness() const { return skewness_; }

  private:
    FloatType skewness_;
};

template<typename FloatType, class Previous=null_accumulator<FloatType> >
class kurtosis_accumulator : public Previous
{
  public:
    kurtosis_accumulator(FloatType mean, FloatType standard_deviation)
      : Previous(mean, standard_deviation),
        kurtosis_(0)
    {}

    void operator()(FloatType x) {
      Previous::operator()(x);
      std::size_t n = this->count();
      FloatType u = this->normalised_deviation();
      kurtosis_ += (u*u*u*u - kurtosis_)/n;
    }

    FloatType kurtosis() const { return kurtosis_; }
    FloatType kurtosis_excess() const { return kurtosis() - 3; }

  private:
    FloatType kurtosis_;
};

template<typename FloatType, class Previous=null_accumulator<FloatType> >
class mean_absolute_deviation_accumulator : public Previous
{
  public:
    mean_absolute_deviation_accumulator(FloatType mean, FloatType width)
      : Previous(mean, width),
        mean_absolute_deviation_(0)
    {}

    mean_absolute_deviation_accumulator(FloatType mean)
      : Previous(mean),
        mean_absolute_deviation_(0)
    {}

    void operator()(FloatType x) {
      Previous::operator()(x);
      std::size_t n = this->count();
      FloatType delta = this->deviation();
      mean_absolute_deviation_ += (std::abs(delta) - mean_absolute_deviation_)/n;
    }

    FloatType mean_absolute_deviation() const {
      return mean_absolute_deviation_;
    }

    void reset() {
      Previous::reset();
      mean_absolute_deviation_ = 0;
    }

  private:
    FloatType mean_absolute_deviation_;

};

/// LAPACK-style norm of a vector: overflow- and underflow-resilient
/** Reference: DLASSQ (LAPACK) or DNRM2 (BLAS shipped with LAPACK) */
template <typename FloatType, class Previous=null_accumulator<FloatType> >
class norm_accumulator : public Previous
{
  public:
    norm_accumulator() : ssq(1), scale(0) {}

    void operator()(FloatType x) {
      Previous::operator()(x);
      if (x == 0) return;
      FloatType absx = std::abs(x);
      if (scale < absx) {
        FloatType t = scale/absx;
        ssq = 1. + ssq*t*t;
        scale = absx;
      }
      else {
        FloatType t = absx/scale;
        ssq += t*t;
      }
    }

    /// Overflow- and underflow-resilient
    FloatType norm() { return scale*std::sqrt(ssq); }

    /// Not overflow- or underflow-resilient
    FloatType sum_sq() { return scale*scale*ssq; }

  private:
    FloatType ssq, scale;
};


template<typename FloatType>
class inertia_accumulator
{
  public:
    inertia_accumulator()
      : sum_weights_(0),
        center_of_mass_(0.),
        m2(0.)
    {}

    inertia_accumulator(vec3<FloatType> const& x0, FloatType weight=1.0)
      : sum_weights_(weight),
        center_of_mass_(x0),
        m2(0.)
    {}

    void operator()(vec3<FloatType> const& x, FloatType weight=1.0) {
      sum_weights_ += weight;
      vec3<FloatType> delta = x - center_of_mass_;
      center_of_mass_ += weight * delta/sum_weights_;
      vec3<FloatType> new_delta = x - center_of_mass_;
      m2[0] += weight * delta[0] * new_delta[0];
      m2[1] += weight * delta[1] * new_delta[1];
      m2[2] += weight * delta[2] * new_delta[2];
      m2[3] += weight * delta[0] * new_delta[1];
      m2[4] += weight * delta[0] * new_delta[2];
      m2[5] += weight * delta[1] * new_delta[2];
    }

    vec3<FloatType> const& center_of_mass() const { return center_of_mass_; }

    /*! inertia_tensor = sum_weights * (identity * trace(covariance) - covariance)
          or
        inertia_tensor = identity * trace(m2) - m2

        see also: http://en.wikipedia.org/wiki/Variance#Moment_of_inertia
     */
    sym_mat3<FloatType> inertia_tensor() const {
      FloatType trace = m2.trace();
      sym_mat3<FloatType> result(trace,trace,trace,0,0,0);
      result -= m2;
      return result;
    }

    //*! The inertia tensor for about the centre of mass, with a correction
    //    applied for inertia around the given pivot point, by taking into
    //    account the Parallel Axis Theorem.
    //    http://en.wikipedia.org/wiki/Moment_of_inertia#Parallel_axis_theorem_2
    // */
    sym_mat3<FloatType> inertia_tensor(vec3<FloatType> const& pivot) const {
      if (sum_weights_ == 0) return sym_mat3<FloatType>(0.);
      vec3<FloatType> p = center_of_mass() - pivot;
      sym_mat3<FloatType> result = inertia_tensor();
      result += sum_weights_ * sym_mat3<FloatType>(
        p[1]*p[1] + p[2]*p[2],
        p[0]*p[0] + p[2]*p[2],
        p[0]*p[0] + p[1]*p[1],
        -p[0]*p[1],
        -p[0]*p[2],
        -p[1]*p[2]);
      return result;
    }

    //! covariance_matrix = 1/sum_weights * m2
    sym_mat3<FloatType> covariance_matrix() const {
      if (sum_weights_ == 0) return sym_mat3<FloatType>(0.);
      return m2 * (1/sum_weights_);
    }

    FloatType sum_weights() const { return sum_weights_; }

  private:
    FloatType sum_weights_;
    vec3<FloatType> center_of_mass_;
    //! m2 = sum(weights) * covariance matrix
    sym_mat3<FloatType> m2;
};

}}} // namespace scitbx::math::accumulator

#endif // GUARD
