#ifndef SCITBX_RANDOM_LEGACY_BOOST_1_63_H
#define SCITBX_RANDOM_LEGACY_BOOST_1_63_H

#include <scitbx/random.h>

namespace boost_1_63 { namespace random {

template<class RealType = double>
class exponential_distribution
{
public:
    typedef RealType result_type;
    /**
     * Constructs an exponential_distribution with a given lambda.
     *
     * Requires: lambda > 0
     */
    explicit exponential_distribution(RealType lambda_arg = RealType(1.0))
      : _lambda(lambda_arg) { BOOST_ASSERT(_lambda > RealType(0)); }
    /**
     * Returns a random variate distributed according to the
     * exponential distribution.
     */
    template<class Engine>
    result_type operator()(Engine& eng) const
    {
        using std::log;
        return -result_type(1) /
            _lambda * log(result_type(1)-boost::random::uniform_01<RealType>()(eng));
    }
private:
    result_type _lambda;
};
} // namespace random
} // namespace boost_1_63

namespace boost_1_70 { namespace random {
namespace detail {

template<class RealType = double>
struct unit_normal_distribution: public boost::random::detail::unit_normal_distribution<RealType>
{
    template<class Engine>
    RealType operator()(Engine& eng) {
        const double * const table_x = boost::random::detail::normal_table<double>::table_x;
        const double * const table_y = boost::random::detail::normal_table<double>::table_y;
        for(;;) {
            std::pair<RealType, int> vals = boost::random::detail::generate_int_float_pair<RealType, 8>(eng);
            int i = vals.second;
            int sign = (i & 1) * 2 - 1;
            i = i >> 1;
            RealType x = vals.first * RealType(table_x[i]);
            if(x < table_x[i + 1]) { return x * sign; }
            if(i == 0) { return generate_tail(eng) * sign; }

            RealType y01 = boost::random::uniform_01<RealType>()(eng);
            RealType y = RealType(table_y[i]) + y01 * RealType(table_y[i + 1] - table_y[i]);

            // These store the value y - bound, or something proportional to that difference:
            RealType y_above_ubound, y_above_lbound;

            // There are three cases to consider:
            // - convex regions (where x[i] > x[j] >= 1)
            // - concave regions (where 1 <= x[i] < x[j])
            // - region containing the inflection point (where x[i] > 1 > x[j])
            // For convex (concave), exp^(-x^2/2) is bounded below (above) by the tangent at
            // (x[i],y[i]) and is bounded above (below) by the diagonal line from (x[i+1],y[i+1]) to
            // (x[i],y[i]).
            //
            // *If* the inflection point region satisfies slope(x[i+1]) < slope(diagonal), then we
            // can treat the inflection region as a convex region: this condition is necessary and
            // sufficient to ensure that the curve lies entirely below the diagonal (that, in turn,
            // also implies that it will be above the tangent at x[i]).
            //
            // For the current table size (128), this is satisfied: slope(x[i+1]) = -0.60653 <
            // slope(diag) = -0.60649, and so we have only two cases below instead of three.

            if (table_x[i] >= 1) { // convex (incl. inflection)
                y_above_ubound = RealType(table_x[i] - table_x[i+1]) * y01 - (RealType(table_x[i]) - x);
                y_above_lbound = y - (RealType(table_y[i]) + (RealType(table_x[i]) - x) * RealType(table_y[i]) * RealType(table_x[i]));
            }
            else { // concave
                y_above_lbound = RealType(table_x[i] - table_x[i+1]) * y01 - (RealType(table_x[i]) - x);
                y_above_ubound = y - (RealType(table_y[i]) + (RealType(table_x[i]) - x) * RealType(table_y[i]) * RealType(table_x[i]));
            }

            if (y_above_ubound < 0 // if above the upper bound reject immediately
                    &&
                    (
                     y_above_lbound < 0 // If below the lower bound accept immediately
                     ||
                     y < f(x) // Otherwise it's between the bounds and we need a full check
                    )
               ) {return x * sign;
            }
        }
    }
    static RealType f(RealType x) { using std::exp; return exp(-(x*x/2));
    }
    template<class Engine>
    RealType generate_tail(Engine& eng) {
        const RealType tail_start = RealType(boost::random::detail::normal_table<double>::table_x[1]);
        boost_1_63::random::exponential_distribution<RealType> exp_x(tail_start);
        boost_1_63::random::exponential_distribution<RealType> exp_y;
        for(;;) {
            RealType x = exp_x(eng);
            RealType y = exp_y(eng);
            // If we were doing non-transformed rejection sampling, this condition would be:
            // if (unif01 < exp(-.5*x*x)) return x + tail_start;
            if(2*y > x*x) return x + tail_start;
        }
    }
};
} // namespace detail

template<class RealType = double>
class normal_distribution
{
public:
    typedef RealType result_type;

    explicit normal_distribution(const RealType& mean_arg = RealType(0.0),
                                 const RealType& sigma_arg = RealType(1.0))
      : _mean(mean_arg), _sigma(sigma_arg)
    {
        BOOST_ASSERT(_sigma >= RealType(0));
    }

    /**  Returns a normal variate. */
    template<class Engine>
    result_type operator()(Engine& eng)
    {
        detail::unit_normal_distribution<RealType> impl;
        return impl(eng) * _sigma + _mean;
    }

private:
    RealType _mean, _sigma;
};
} // namespace random
} // namespace boost_1_70

namespace scitbx {

//! Easy access to Boost.Random.
/*! See also: http://www.boost.org/libs/random/
 */
namespace random_legacy_boost_1_63 {

  //! Wrapper for boost/random/mersenne_twister.hpp
  /*! See also: http://www.boost.org/libs/random/
   */
  class mersenne_twister: public scitbx::random::mersenne_twister
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

      /*! \brief Random quaternion giving identical results
           as boost 1.63 algorithm
       */
      af::tiny<double, 4>
      random_double_unit_quaternion()
      {
        // pick a uniformly distributed, random point on 3-sphere
        // http://mathworld.wolfram.com/HyperspherePointPicking.html
        af::tiny<double, 4> result;
        boost_1_70::random::normal_distribution<double> distribution(0.0, 1.0);
        boost::variate_generator<
          boost_random::mt19937&,
          boost_1_70::random::normal_distribution<double> >
            normal_distribution_random_double(generator_, distribution);
        double r = 0.0;
        do {
          for (int i=0; i<4; i++) {
            result[i] = normal_distribution_random_double();
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

      /*! \brief Uniformly distributed random 3D rotation matrix
           giving identical results as using boost 1.63 algorithm
       */
      mat3<double>
      random_double_r3_rotation_matrix()
      {
        return math::r3_rotation::unit_quaternion_as_matrix(
          random_double_unit_quaternion());
      }

    private:
      boost_random::mt19937 generator_;
  };

}} // namespace scitbx::random_legacy_boost_1_63

#endif // SCITBX_RANDOM_LEGACY_BOOST_1_63_H
