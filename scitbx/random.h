#ifndef SCITBX_RANDOM_H
#define SCITBX_RANDOM_H

#include <scitbx/random/mersenne_twister.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/math/utils.h>
#include <scitbx/array_family/shared.h>
#include <stdexcept>

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
      //! Initialization with generator instance
      mersenne_twister(boost_random::mt19937 &generator)
      :
        generator_(generator) {}

      //! Initialization with given seed.
      mersenne_twister(unsigned seed=0)
      :
        generator_(seed+1) {}

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
        std::size_t* r = result.begin();
        for(std::size_t i=0;i<size;i++) *r++ = random_size_t();
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
        std::size_t* r = result.begin();
        for(std::size_t i=0;i<size;i++) *r++ = random_size_t() % modulus;
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
        double* r = result.begin();
        for(std::size_t i=0;i<size;i++) *r++ = random_double();
        return result;
      }

      /*! \brief Array of uniformly distributed random doubles in the
          range [0, factor).
       */
      af::shared<double>
      random_double(std::size_t size, double factor)
      {
        af::shared<double> result(size, af::init_functor_null<double>());
        double* r = result.begin();
        for(std::size_t i=0;i<size;i++) *r++ = random_double() * factor;
        return result;
      }

      //! \brief Array of results: random_double() < threshold
      af::shared<bool>
      random_bool(std::size_t size, double threshold)
      {
        af::shared<bool> result(size, af::init_functor_null<bool>());
        bool *r = result.begin();
        for(std::size_t i=0;i<size;i++) *r++ = random_double() < threshold;
        return result;
      }

      //! Random permutation of integers in the range [0, size).
      af::shared<std::size_t>
      random_permutation(std::size_t size)
      {
        af::shared<std::size_t> result(
          size, af::init_functor_null<std::size_t>());
        std::size_t* r = result.begin();
        for(std::size_t i=0;i<size;i++) *r++ = i;
        r = result.begin();
        for(std::size_t i=0;i<size;i++) {
          std::size_t j = static_cast<std::size_t>(generator_()) % size;
          std::swap(r[i], r[j]);
        }
        return result;
      }

      //! Uniform random points on 3D sphere.
      /*! http://cgafaq.info/wiki/Random_Points_On_Sphere (2006_08_30)
          (probably by Colas Schretter)
          Trig. method
            This method works only in 3-space, but it is very fast. It
            depends on the slightly counterintuitive fact (see proof
            below) that each of the three coordinates is uniformly
            distributed on [-1,1] (but the three are not independent,
            obviously). Therefore, it suffices to choose one axis
            (Z, say) and generate a uniformly distributed value on
            that axis. This constrains the chosen point to lie on a
            circle parallel to the X-Y plane, and the obvious trig
            method may be used to obtain the remaining coordinates.

              1. Choose z uniformly distributed in [-1,1].
              2. Choose t uniformly distributed on [0, 2 p).
              3. Let r = sqrt(1-z**2).
              4. Let x = r * cos(t).
              5. Let y = r * sin(t).
       */
      vec3<double>
      random_double_point_on_sphere()
      {
        vec3<double> result;
        double z = 2 * random_double() - 1;
        double t = constants::two_pi * random_double();
        double r = std::sqrt(1-z*z);
        result[0] = r * std::cos(t);
        result[1] = r * std::sin(t);
        result[2] = z;
        return result;
      }

      af::tiny<double, 4>
      random_double_unit_quaternion()
      {
        // pick a uniformly distributed, random point on 3-sphere
        // http://mathworld.wolfram.com/HyperspherePointPicking.html
        af::tiny<double, 4> result;
        double r = 0.0;
        do {
          for (int i=0; i<4; i++) {
            result[i] = random_double();
            r += result[i] * result[i];
          }
          r = std::sqrt(r);
        }
        while (r == 0.0);
        for (int i=0; i<4; i++) {
          result[i] /= r;
        }
        return result;
      }

      //! Uniformly distributed random 3D rotation matrix using quaternions
      mat3<double>
      random_double_r3_rotation_matrix_quaternion()
      {
        return math::r3_rotation::unit_quaternion_as_matrix(
          random_double_unit_quaternion());
      }

      mat3<double>
      random_double_r3_rotation_matrix()
      {
        vec3<double> axis = random_double_point_on_sphere();
        double angle = constants::two_pi * random_double();
        return math::r3_rotation::axis_and_angle_as_matrix(
          axis, angle, /* deg */ false);
      }

      //! Uniformly distributed random 3D rotation matrix using Arvo's method.
      /*! See also: scitbx::math::r3_rotation::random_matrix_arvo_1992
       */
      mat3<double>
      random_double_r3_rotation_matrix_arvo_1992()
      {
        /* Results are not predictable if the calls of random_double()
           are inlined into the function call. This is because the compiler
           is free to generate code that evaluates the arguments in any order.
         */
        double x0 = random_double();
        double x1 = random_double();
        double x2 = random_double();
        return math::r3_rotation::random_matrix_arvo_1992(x0, x1, x2);
      }

      //! Integer array with Gaussian distribution.
      /*! Algorithm copied from the Python source code.
       */
      af::shared<int>
      random_int_gaussian_distribution(
        std::size_t size,
        double const& mean,
        double const& sigma)
      {
        af::shared<int> result(size, af::init_functor_null<int>());
        int* r = result.begin();
        for(std::size_t i=0;i<size;i++) {
          double x2pi = random_double() * 2.0 * constants::pi;
          double g2rad = std::sqrt(-2.0 * std::log(1.-random_double()));
          *r++ = math::nearest_integer(mean + std::cos(x2pi) * g2rad * sigma);
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
