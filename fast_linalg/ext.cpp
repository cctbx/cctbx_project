#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/object.hpp>

#include <fast_linalg/environment.h>
#include <fast_linalg/lapacke.h>

namespace fast_linalg {

namespace {

  void init_module() {
    using namespace boost::python;
    typedef environment e;
    class_<environment>("environment")
      .add_property("threads", &e::threads, &e::set_threads,
                    "Number of threads to be used in parallelisation")
      .add_property("physical_cores", &e::physical_cores,
                    "Number of physical cores on the machine\n\n"
                    "E.g. a machine with 2 hexacore processors will "
                    "reports 12 cores\n(with or without hyperthreading).")
      .add_property("cpu_family", &e::cpu_family,
                    "The family of the CPU\n\n"
                    "E.g. \"Nehalem\" for an Intel processor.")
      .add_property("build_config", &e::build_config)
      .add_static_property("initialised", &fast_linalg::is_initialised)
      .def("initialise", &fast_linalg::initialise, (arg("lib_name")))
      .def("finalise", &fast_linalg::finalise)
      .staticmethod("initialise")
      ;
    scope().attr("env") = object(new environment());
  }

}}

BOOST_PYTHON_MODULE(fast_linalg_ext) {
  fast_linalg::init_module();
}
