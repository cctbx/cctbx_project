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
        .def("match_cached", &w_t::match_cached,
           "To be used when the match_indices instance was constructed from \n \
           one Miller array h0. Upon calling match_indices.match_cached(h1),\n \
           the matches between h0 and h1 are found and stored as tuples\n \
           (i0, i1) where i0 and i1 are indices into the arrays h0 and h1.\n \
           Retrieve the array of (i0, i1) tuples by calling\n \
           match_indices.pairs()."
            )
        .def("match_cached_fast", &w_t::match_cached_fast,
           "See match_indices.match_cached. This is faster when only \n \
           the pairs (i0, i1) are needed. Other member functions such as \n \
           match_cached.singles(i) are not available when the matching was \n \
           done via match_cached_fast. Note that match_cached_fast \n \
           populates the pairs array in ascending order of indices i1 \n \
           instead of i0."
            )
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
