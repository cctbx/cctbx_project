#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/crystal/close_packing.h>

namespace cctbx { namespace crystal { namespace close_packing {

namespace {

  struct hexagonal_sampling_wrappers
  {
    typedef hexagonal_sampling<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("close_packing_hexagonal_sampling", no_init)
        .def(init<direct_space_asu::float_asu<double>,
                  af::tiny<bool, 3> const&,
                  double const&,
                  optional<double const&,
                           bool> >(
          (arg_("float_asu"),
           arg_("continuous_shift_flags"),
           arg_("point_distance"),
           arg_("buffer_thickness")=-1,
           arg_("all_twelve_neighbors")=false)))
        .def("float_asu", &w_t::float_asu, rir())
        .def("continuous_shift_flags", &w_t::continuous_shift_flags, ccr())
        .def("point_distance", &w_t::point_distance, ccr())
        .def("buffer_thickness", &w_t::buffer_thickness, ccr())
        .def("all_twelve_neighbors", &w_t::all_twelve_neighbors)
        .def("box_lower", &w_t::box_lower, ccr())
        .def("box_upper", &w_t::box_upper, ccr())
        .def("all_sites_frac", &w_t::all_sites_frac)
      ;
    }
  };

}} // namespace close_packing::<anoymous>

namespace boost_python {

  void wrap_close_packing()
  {
    close_packing::hexagonal_sampling_wrappers::wrap();
  }

}}} // namespace cctbx::crystal::boost_python
