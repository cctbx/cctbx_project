#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/match_bijvoet_mates.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct match_bijvoet_mates_wrappers
  {
    typedef match_bijvoet_mates w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("match_bijvoet_mates", no_init)
        .def(init<sgtbx::space_group_type const&,
                  af::shared<index<> > const&,
                  bool>((
                    arg("sg_type"), arg("indices"),
                    arg("assert_is_unique_set_under_symmetry")=true)))
        .def(init<sgtbx::reciprocal_space::asu const&,
                  af::shared<index<> > const&,
                  bool>((
                    arg("asu"), arg("indices"),
                    arg("assert_is_unique_set_under_symmetry")=true)))
        .def(init<af::shared<index<> > const&,
                  bool>((
                    arg("indices"),
                    arg("assert_is_unique_set_under_symmetry")=true)))
        .def("pairs", &w_t::pairs)
        .def("singles", &w_t::singles)
        .def("n_singles", &w_t::n_singles)
        .def("pairs_hemisphere_selection", &w_t::pairs_hemisphere_selection)
        .def("singles_hemisphere_selection",
          &w_t::singles_hemisphere_selection, ccr())
        .def("miller_indices_in_hemisphere",
          &w_t::miller_indices_in_hemisphere)
#define CCTBX_DEF(function_name) \
        .def(# function_name, \
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const) \
          &w_t::function_name)
        CCTBX_DEF(minus)
        CCTBX_DEF(additive_sigmas)
        CCTBX_DEF(average)
#undef CCTBX_DEF
      ;
    }
  };

} // namespace <anoymous>

  void wrap_match_bijvoet_mates()
  {
    match_bijvoet_mates_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
