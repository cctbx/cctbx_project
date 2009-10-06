#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/random/variate_generator.h>
#include <scitbx/random/mersenne_twister.h>
#include <boost/random/normal_distribution.hpp>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace random { namespace boost_python {

namespace {

  template<class Engine, class Distribution>
  struct variate_generator_wrappers
  {
    typedef scitbx::random::variate_generator<Engine, Distribution> wt;
    typedef typename wt::result_type FloatType;

    static void
    wrap(char const *name)
    {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<Engine, Distribution>(args("engine", "distribution")))
        .def("__call__", (FloatType(wt::*)()) &wt::operator())
        .def("__call__", (scitbx::af::shared<FloatType>(wt::*)(std::size_t))
          &wt::operator(), arg("size"))
      ;
    }
  };

  template<typename FloatType>
  struct normal_distribution_wrappers
  {
    typedef boost::normal_distribution<FloatType> wt;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<wt>("normal_distribution", no_init)
        .def(init<FloatType, FloatType>((arg("mean")=FloatType(0), arg("sigma")=FloatType(1))))
        .def("mean", &wt::mean)
        .def("sigma", &wt::sigma)
        .def("reset", &wt::reset)
      ;
    }
  };

  /*! This is intentionally kept useless as a generator
      in order to discourage use directly from python.
      Instead use variate_generators.
   */
  struct mt19937_wrappers
  {
    typedef boost_random::mt19937 wt;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<wt>("mt19937", no_init)
        .def(init<wt::result_type>(arg("value")))
        .def("seed", (void(wt::*)()) &wt::seed)
        .def("seed", (void(wt::*)(const wt::result_type&)) &wt::seed)
      ;
    }
  };

} // namespace <anonymous>

  void wrap_random()
  {
    mt19937_wrappers::wrap();
    normal_distribution_wrappers<double>::wrap();
    variate_generator_wrappers<
      boost_random::mt19937 &,
      boost::normal_distribution<double>
    >::wrap("normal_variate_generator");
  }

}}} // namespace scitbx::random::boost_python
