#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/match_indices.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct match_indices_wrappers
  {
    typedef match_indices w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("match_indices", no_init)
        .def(init<af::shared<index<> > const&,
                  af::shared<index<> > const&>())
        .def(init<af::shared<index<> > const&>())
        .def("match_cached", &w_t::match_cached)
        .def("match_cached_fast", &w_t::match_cached_fast)
        .def("pairs", &w_t::pairs)
        .def("singles", &w_t::singles)
        .def("have_singles", &w_t::have_singles)
        .def("pair_selection", &w_t::pair_selection)
        .def("single_selection", &w_t::single_selection)
        .def("paired_miller_indices", &w_t::paired_miller_indices)
        .def("permutation", &w_t::permutation)
#define CCTBX_DEF(function_name) \
        .def(# function_name, \
          (af::shared<double>(w_t::*)(af::const_ref<double> const&, \
                                      af::const_ref<double> const&) const) \
          &w_t::function_name) \
        .def(# function_name, \
          (af::shared<std::complex<double> >(w_t::*)( \
            af::const_ref<std::complex<double> > const&, \
            af::const_ref<std::complex<double> > const&) const) \
          &w_t::function_name)
        CCTBX_DEF(plus)
        CCTBX_DEF(minus)
        CCTBX_DEF(divides)
        CCTBX_DEF(multiplies)
        CCTBX_DEF(additive_sigmas)
#undef CCTBX_DEF
      ;
    }
  };

} // namespace <anoymous>

  void wrap_match_indices()
  {
    match_indices_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
