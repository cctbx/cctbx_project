#ifndef SCITBX_RANDOM_H
#define SCITBX_RANDOM_H

#include <boost/random/mersenne_twister.hpp>
#include <scitbx/array_family/shared.h>
#include <sstream>

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
        return static_cast<std::size_t>(generator_.min());
      }

      //! Largest value returned by random_size_t().
      std::size_t
      random_size_t_max()
      {
        return static_cast<std::size_t>(generator_.max());
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

      std::string
      getstate() const
      {
        std::stringstream s;
        s << generator_;
        return s.str();
      }

      void
      setstate(std::string const& state)
      {
        std::stringstream s(state);
        s >> generator_;
      }

    private:
      boost::mt19937 generator_;
  };

}} // namespace scitbx::random
#endif // SCITBX_RANDOM_H
