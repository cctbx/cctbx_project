#include <boost/python/module.hpp>
#include <boost/python/class.hpp>

#include <scitbx/minpack/levenberg_marquardt.h>

namespace scitbx { namespace minpack {
namespace {

  struct levenberg_marquardt_wrappers
  {
    typedef levenberg_marquardt w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("levenberg_marquardt", no_init)
        .def(init<
          int,
          af::shared<double>,
          optional<
            double,
            double,
            double,
            int,
            double,
            bool> >((
          arg_("m"),
          arg_("x"),
          arg_("ftol")=-1,
          arg_("xtol")=-1,
          arg_("gtol")=0,
          arg_("maxfev")=0,
          arg_("factor")=1.0e2,
          arg_("call_back_after_iteration")=false)))
        .def_readonly("m", &w_t::m)
        .def_readonly("ftol", &w_t::ftol)
        .def_readonly("xtol", &w_t::xtol)
        .def_readonly("gtol", &w_t::gtol)
        .def_readonly("maxfev", &w_t::maxfev)
        .def_readonly("factor", &w_t::factor)
        .def("has_terminated", &w_t::has_terminated)
        .def("requests_fvec", &w_t::requests_fvec)
        .def("requests_fjac", &w_t::requests_fjac)
        .def("calls_back_after_iteration", &w_t::calls_back_after_iteration)
        .def("process_fvec", &w_t::process_fvec, (arg_("fvec")))
        .def("process_fjac", &w_t::process_fjac, (arg_("fjac")))
        .def("continue_after_call_back_after_iteration",
          &w_t::continue_after_call_back_after_iteration)
        .def_readonly("info", &w_t::info)
        .def_readonly("nfev", &w_t::nfev)
        .def_readonly("njev", &w_t::njev)
      ;
    }
  };

  void init_module()
  {
    levenberg_marquardt_wrappers::wrap();
  }

}}} // namespace scitbx::minpack::<anonymous>

BOOST_PYTHON_MODULE(scitbx_minpack_ext)
{
  scitbx::minpack::init_module();
}
