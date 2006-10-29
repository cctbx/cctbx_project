#ifndef SCITBX_FFTPACK_FACTORIZATION_H
#define SCITBX_FFTPACK_FACTORIZATION_H

#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace fftpack {

  /*! \brief Determination of prime factors for both complex-to-complex
       and real-to-complex transforms.
   */
  class factorization
  {
    public:
      //! Default constructor.
      factorization() : n_(0) {}
      //! Determination of the %factorization of N.
      /*! Computation of the prime factors of N with the specialization
          that one length-4 transform is computed instead of two length-2
          transforms.
       */
      factorization(std::size_t n, bool real_to_complex);
      //! Access the n that was passed to the constructor.
      std::size_t n() const { return n_; }
      //! Access the factors of n.
      af::shared<int> factors() const { return factors_; }
    protected:
      std::size_t n_;
      af::shared<int> factors_;
  };

  namespace detail {

    template <typename IntegerType>
    IntegerType
    count_reduce(IntegerType& red_n, const IntegerType& factor)
    {
      IntegerType result = 0;
      while (red_n % factor == 0) {
        red_n /= factor;
        result++;
      }
      return result;
    }

  } // namespace detail

  inline // Potential NOINLINE
  factorization::factorization(std::size_t n, bool real_to_complex)
    : n_(n)
  {
    // Based on the first parts of FFTPACK41 cffti1.f and rffti1.f.
    const af::int3 opt_factors(3, 4, 2);
    af::int3 perm(2, 0, 1);
    if (real_to_complex) { // XXX does this make a difference?
      perm[1] = 1;
      perm[2] = 0;
    }
    af::int3 count;
    count.fill(0);
    int red_n = n_;
    int i;
    for (i = 0; red_n > 1 && i < opt_factors.size(); i++) {
      count[i] = detail::count_reduce(red_n, opt_factors[i]);
    }
    for (i = 0; i < opt_factors.size(); i++) {
      factors_.insert(factors_.end(), count[perm[i]], opt_factors[perm[i]]);
    }
    for (int factor = 5; red_n > 1; factor += 2) {
      factors_.insert(
        factors_.end(), detail::count_reduce(red_n, factor), factor);
    }
  }

}} // namespace scitbx::fftpack

#endif // SCITBX_FFTPACK_FACTORIZATION_H
