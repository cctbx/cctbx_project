#include <scitbx/random/mersenne_twister.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace random { namespace boost_python {

  void wrap_random();

  /*! This is intentionally kept useless as a generator
   in order to discourage use directly from python.
   Instead use variate_generators.
   */
  struct mt19937_wrappers
  {
    typedef boost_random::mt19937 wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt>("mersenne_twister_19937", no_init)
        .def(init<wt::result_type>(arg("value")))
        .def("seed", (void(wt::*)()) &wt::seed)
        .def("seed", (void(wt::*)(const wt::result_type&)) &wt::seed)
        ;
    }
  };

  namespace {
    void init_module() {
      mt19937_wrappers::wrap();
      wrap_random();
    }
  }
}}}

BOOST_PYTHON_MODULE(scitbx_random_ext)
{
  scitbx::random::boost_python::init_module();
}
