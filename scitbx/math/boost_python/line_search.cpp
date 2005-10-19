#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <scitbx/line_search/more_thuente_1994.h>

namespace scitbx { namespace line_search { namespace boost_python {

namespace {

  struct more_thuente_1994_wrappers
  {
    typedef more_thuente_1994<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("line_search_more_thuente_1994")
        .def_readwrite("xtol", &w_t::xtol)
        .def_readwrite("ftol", &w_t::ftol)
        .def_readwrite("gtol", &w_t::gtol)
        .def_readwrite("stpmin", &w_t::stpmin)
        .def_readwrite("stpmax", &w_t::stpmax)
        .def_readwrite("maxfev", &w_t::maxfev)
        .def_readonly("info_code", &w_t::info_code)
        .def_readonly("info_meaning", &w_t::info_meaning)
        .def_readonly("stp", &w_t::stp)
        .def_readonly("nfev", &w_t::nfev)
        .def("start", &w_t::start, (
          arg_("x"),
          arg_("functional"),
          arg_("gradients"),
          arg_("search_direction"),
          arg_("initial_estimate_of_satisfactory_step_length")))
        .def("next", &w_t::next, (
          arg_("x"),
          arg_("functional"),
          arg_("gradients")))
      ;
    }
  };

  void
  wrap()
  {
    more_thuente_1994_wrappers::wrap();
  }

}}} // namespace line_search::boost_python::<anonymous>

namespace math { namespace boost_python {

  void wrap_line_search()
  {
    line_search::boost_python::wrap();
  }

}}} // namespace scitbx::math::boost_python
