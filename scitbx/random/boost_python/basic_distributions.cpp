#include <scitbx/random/boost_python/random.h>

#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>


namespace scitbx { namespace random { namespace boost_python {

namespace {

  struct uniform
  {
    typedef boost::uniform_real<double> wt;

    static std::string name() { return "uniform"; }

    static void wrap_specific(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      klass
        .def(init<double, double>((arg("min"), arg("max"))))
        .add_property("min", &wt::min)
        .add_property("min", &wt::min)
      ;
    }
  };

  struct normal
  {
    typedef boost::normal_distribution<double> wt;

    static std::string name() { return "normal"; }

    static void wrap_specific(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      klass
        .def(init<double, double>((arg("mean") =0.,
                                   arg("sigma")=1.)))
        .add_property("mean", &wt::mean)
        .add_property("sigma", &wt::sigma)
        ;
    }
  };

  struct bernoulli
  {
    typedef boost::bernoulli_distribution<double> wt;

    static std::string name() { return "bernoulli"; }

    static void wrap_specific(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      klass
        .def(init<double>())
        .add_property("p", &wt::p)
        ;
    }
  };

  struct gamma
  {
    typedef boost::gamma_distribution<double> wt;

    static std::string name() { return "gamma"; }

    static void wrap_specific(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      klass
        .def(init<double, double>((arg("alpha")=1.0,
                                   arg("beta")=1.0)))
        .add_property("alpha", &wt::alpha)
        .add_property("beta", &wt::beta)
        ;
    }
  };


} // namespace <anonymous>

  void wrap_random()
  {
    wrap_distribution_and_variate<uniform>();
    wrap_distribution_and_variate<normal>();
    wrap_distribution_and_variate<bernoulli>();
    wrap_distribution_and_variate<gamma>();
  }

}}} // namespace scitbx::random::boost_python
