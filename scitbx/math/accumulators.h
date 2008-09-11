#ifndef SCITBX_MATH_ACCUMULATORS
#define SCITBX_MATH_ACCUMULATORS

#include <cmath>
#include <cstddef>

namespace scitbx { namespace math {

///Accumulators to compute statistics.
/**
The classes in this namespace are designed to help separating the computation
of statistics from the iteration over the sequence of values.

See class \link basic_statistics \endlink for a typical example of their use.

*/
namespace accumulator {

template<typename FloatType>
class null_accumulator
{
  public:
    null_accumulator(FloatType x0) {}
    void operator()(FloatType x) {}
};


template<typename FloatType>
class enumerated_accumulator
{
  public:
    enumerated_accumulator(FloatType x0)
      : n(1)
    {}

    void operator()(FloatType x) { n++; }

    std::size_t count() const { return n; }

  private:
    std::size_t n;

};


template<typename FloatType, class Previous>
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

  private:
    FloatType min_, max_, max_absolute_;
};


template<typename FloatType, class Previous>
class mean_variance_accumulator : public Previous
{
  public:
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


template<typename FloatType, class Previous>
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

template<typename FloatType, class Previous>
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

template<typename FloatType, class Previous>
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

  private:
    FloatType mean_absolute_deviation_;
};


}}} // namespace scitbx::math::accumulator

#endif // GUARD
