#ifndef SCITBX_RANDOM_VARIATE_GENERATOR_H
#define SCITBX_RANDOM_VARIATE_GENERATOR_H

#include <scitbx/array_family/shared.h>
#include <boost/random/variate_generator.hpp>

namespace scitbx { namespace random {

  /*! An extension of boost::variate_generator.
      See also: http://www.boost.org/libs/random/random-variate.html
      Operator() overloaded to generate array of random numbers
      in addition to single numbers.
   */
  template<class Engine, class Distribution>
  struct variate_generator
        : public boost::variate_generator<Engine, Distribution>
  {
  private:
    typedef boost::variate_generator<Engine, Distribution> base_class;

  public:
    typedef typename base_class::result_type result_type;

    variate_generator(Engine e, Distribution d) : base_class(e, d) {}

    result_type operator()() {
      return base_class::operator()();
    }

    af::shared<result_type> operator()(std::size_t size) {
      af::reserve reserved(size);
      af::shared<result_type> result(reserved);
      for(std::size_t i=0; i<size; ++i) result.push_back((*this)());
      return result;
    }
  };

}} // scitbx::random

#endif // GUARD
