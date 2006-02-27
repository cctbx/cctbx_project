#ifndef SCITBX_RANDOM_H
#define SCITBX_RANDOM_H

// Simplified copy of boost/boost/random/mersenne_twister.hpp
//   The main motivation for the copy is to get access to
//   the state without the iostream convolution.
/*
 * Copyright Jens Maurer 2000-2001
 * Distributed under the Boost Software License, Version 1.0. (See
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#include <scitbx/array_family/shared.h>
#include <boost/cstdint.hpp>
#include <stdexcept>

namespace scitbx {
namespace boost_random {

// http://www.math.keio.ac.jp/matumoto/emt.html
template<class UIntType, int w, int n, int m, int r, UIntType a, int u,
  int s, UIntType b, int t, UIntType c, int l, UIntType val>
class mersenne_twister
{
public:
  typedef UIntType result_type;
  static const int word_size = w;
  static const int state_size = n;
  static const int shift_size = m;
  static const int mask_bits = r;
  static const UIntType parameter_a = a;
  static const int output_u = u;
  static const int output_s = s;
  static const UIntType output_b = b;
  static const int output_t = t;
  static const UIntType output_c = c;
  static const int output_l = l;

  static const bool has_fixed_range = false;

  mersenne_twister() { seed(); }

  explicit mersenne_twister(const UIntType& value)
  { seed(value); }

  void seed() { seed(UIntType(5489)); }

  void seed(const UIntType& value)
  {
    // New seeding algorithm from
    // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
    // In the previous versions, MSBs of the seed affected only MSBs of the
    // state x[].
    const UIntType mask = ~0u;
    x[0] = value & mask;
    for (i = 1; i < n; i++) {
      // See Knuth "The Art of Computer Programming" Vol. 2, 3rd ed., page 106
      x[i] = (1812433253UL * (x[i-1] ^ (x[i-1] >> (w-2))) + i) & mask;
    }
  }

  result_type min_value() const { return 0; }
  result_type max_value() const
  {
    // avoid "left shift count >= with of type" warning
    result_type res = 0;
    for(int i = 0; i < w; ++i)
      res |= (1u << i);
    return res;
  }

  result_type operator()();

  af::shared<std::size_t>
  getstate() const
  {
    af::shared<std::size_t> result;
    result.reserve(n);
    for(unsigned j=0;j<n;j++) {
      result.push_back(compute(j));
    }
    return result;
  }

  void
  setstate(af::const_ref<std::size_t> const& state)
  {
    if (state.size() != n) {
      throw std::runtime_error(
        "mersenne_twister::setstate: improper state.size()");
    }
    for(unsigned j=0;j<n;j++) {
      x[j] = state[j];
    }
    i = n;
  }

private:
  // returns x(i-n+index), where index is in 0..n-1
  UIntType compute(unsigned int index) const
  {
    // equivalent to (i-n+index) % 2n, but doesn't produce negative numbers
    return x[ (i + n + index) % (2*n) ];
  }
  void twist(int block);

  // state representation: next output is o(x(i))
  //   x[0]  ... x[k] x[k+1] ... x[n-1]     x[n]     ... x[2*n-1]   represents
  //  x(i-k) ... x(i) x(i+1) ... x(i-k+n-1) x(i-k-n) ... x[i(i-k-1)]
  // The goal is to always have x(i-n) ... x(i-1) available for
  // operator== and save/restore.

  UIntType x[2*n];
  int i;
};

template<class UIntType, int w, int n, int m, int r, UIntType a, int u,
  int s, UIntType b, int t, UIntType c, int l, UIntType val>
void mersenne_twister<UIntType,w,n,m,r,a,u,s,b,t,c,l,val>::twist(int block)
{
  const UIntType upper_mask = (~0u) << r;
  const UIntType lower_mask = ~upper_mask;

  if(block == 0) {
    for(int j = n; j < 2*n; j++) {
      UIntType y = (x[j-n] & upper_mask) | (x[j-(n-1)] & lower_mask);
      x[j] = x[j-(n-m)] ^ (y >> 1) ^ (y&1 ? a : 0);
    }
  } else if (block == 1) {
    // split loop to avoid costly modulo operations
    {  // extra scope for MSVC brokenness w.r.t. for scope
      for(int j = 0; j < n-m; j++) {
        UIntType y = (x[j+n] & upper_mask) | (x[j+n+1] & lower_mask);
        x[j] = x[j+n+m] ^ (y >> 1) ^ (y&1 ? a : 0);
      }
    }

    for(int j = n-m; j < n-1; j++) {
      UIntType y = (x[j+n] & upper_mask) | (x[j+n+1] & lower_mask);
      x[j] = x[j-(n-m)] ^ (y >> 1) ^ (y&1 ? a : 0);
    }
    // last iteration
    UIntType y = (x[2*n-1] & upper_mask) | (x[0] & lower_mask);
    x[n-1] = x[m-1] ^ (y >> 1) ^ (y&1 ? a : 0);
    i = 0;
  }
}

template<class UIntType, int w, int n, int m, int r, UIntType a, int u,
  int s, UIntType b, int t, UIntType c, int l, UIntType val>
inline typename mersenne_twister<UIntType,w,n,m,r,a,u,s,b,t,c,l,val>::result_type
mersenne_twister<UIntType,w,n,m,r,a,u,s,b,t,c,l,val>::operator()()
{
  if(i == n)
    twist(0);
  else if(i >= 2*n)
    twist(1);
  // Step 4
  UIntType z = x[i];
  ++i;
  z ^= (z >> u);
  z ^= ((z << s) & b);
  z ^= ((z << t) & c);
  z ^= (z >> l);
  return z;
}

// validation by experiment from mt19937.c
typedef mersenne_twister<uint32_t,32,624,397,31,0x9908b0df,11,
  7,0x9d2c5680,15,0xefc60000,18, 3346425566U> mt19937;

} // namespace boost_random
} // namespace scitbx

namespace scitbx {
//! Easy access to Boost.Random.
/*! See also: http://www.boost.org/libs/random/
 */
namespace random {

  //! Wrapper for boost/random/mersenne_twister.hpp
  /*! See also: http://www.boost.org/libs/random/
   */
  class mersenne_twister
  {
    public:
      //! Initialization with given seed.
      mersenne_twister(unsigned seed=0)
      :
        generator_(seed+1)
      {
        //init();
      }

      //! Re-initialization with given seed.
      void
      seed(unsigned value=0) { generator_.seed(value+1); }

      //! Smallest value returned by random_size_t().
      std::size_t
      random_size_t_min()
      {
        return static_cast<std::size_t>(generator_.min_value());
      }

      //! Largest value returned by random_size_t().
      std::size_t
      random_size_t_max()
      {
        return static_cast<std::size_t>(generator_.max_value());
      }

      /*! \brief Uniformly distributed random integer in the range
          [random_size_t_min(), random_size_t_max()].
       */
      std::size_t
      random_size_t()
      {
        return static_cast<std::size_t>(generator_());
      }

      /*! \brief Array of uniformly distributed random integers in the
          range [random_size_t_min(), random_size_t_max()].
       */
      af::shared<std::size_t>
      random_size_t(std::size_t size)
      {
        af::shared<std::size_t> result(
          size, af::init_functor_null<std::size_t>());
        for(std::size_t i=0;i<size;i++) result[i] = random_size_t();
        return result;
      }

      /*! \brief Array of uniformly distributed random integers in the
          range [0, modulus).
       */
      af::shared<std::size_t>
      random_size_t(std::size_t size, std::size_t modulus)
      {
        af::shared<std::size_t> result(
          size, af::init_functor_null<std::size_t>());
        for(std::size_t i=0;i<size;i++) {
          result[i] = random_size_t() % modulus;
        }
        return result;
      }

      /*! \brief Uniformly distributed random double in the range
          [0, 1).
       */
      double
      random_double()
      {
// From Python-2.4.1/Modules/_randommodule.c
/* genrand_res53 in the original code;
 * generates a random number on [0,1) with 53-bit resolution; note that
 * 9007199254740992 == 2**53; 67108864 is 2**26.  In
 * effect, a contains 27 random bits shifted left 26, and b fills in the
 * lower 26 bits of the 53-bit numerator.
 * The orginal code credited Isaku Wada for this algorithm, 2002/01/09.
 */
        std::size_t a = random_size_t() >> 5;
        std::size_t b = random_size_t() >> 6;
        static const double c = 1.0/9007199254740992.0;
        return (a*67108864.0+b)*c;
      }

      /*! \brief Array of uniformly distributed random doubles in the
          range [0, 1).
       */
      af::shared<double>
      random_double(std::size_t size)
      {
        af::shared<double> result(size, af::init_functor_null<double>());
        for(std::size_t i=0;i<size;i++) {
          result[i] = random_double();
        }
        return result;
      }

      /*! \brief Array of uniformly distributed random doubles in the
          range [0, factor).
       */
      af::shared<double>
      random_double(std::size_t size, double factor)
      {
        af::shared<double> result(size, af::init_functor_null<double>());
        for(std::size_t i=0;i<size;i++) {
          result[i] = random_double() * factor;
        }
        return result;
      }

      //! Random permutation of integers in the range [0, size).
      af::shared<std::size_t>
      random_permutation(std::size_t size)
      {
        af::shared<std::size_t> result(
          size, af::init_functor_null<std::size_t>());
        for(std::size_t i=0;i<size;i++) {
          result[i] = i;
        }
        for(std::size_t i=0;i<size;i++) {
          std::size_t j = static_cast<std::size_t>(generator_()) % size;
          std::swap(result[i], result[j]);
        }
        return result;
      }

      af::shared<std::size_t>
      getstate() const { return generator_.getstate(); }

      void
      setstate(af::const_ref<std::size_t> const& state)
      {
        generator_.setstate(state);
      }

    private:
      boost_random::mt19937 generator_;
  };

}} // namespace scitbx::random

#endif // SCITBX_RANDOM_H
