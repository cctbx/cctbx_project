#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <cctbx/crystal/close_packing.h>

namespace cctbx { namespace crystal { namespace close_packing {

namespace {

  struct hexagonal_sampling_wrappers
  {
    typedef hexagonal_sampling<> w_t;

    static fractional<>
    next_site_frac(w_t& o)
    {
      if (o.at_end()) {
        PyErr_SetString(PyExc_StopIteration, "Sampling sites are exhausted.");
        boost::python::throw_error_already_set();
      }
      return o.next_site_frac();
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("close_packing_hexagonal_sampling", no_init)
        .def(init<sgtbx::change_of_basis_op const&,
                  direct_space_asu::float_asu<double>,
                  af::tiny<bool, 3> const&,
                  double const&,
                  optional<double const&,
                           bool> >(
          (arg_("cb_op_original_to_sampling"),
           arg_("float_asu"),
           arg_("continuous_shift_flags"),
           arg_("point_distance"),
           arg_("buffer_thickness")=-1,
           arg_("all_twelve_neighbors")=false)))
        .def("cb_op_original_to_sampling",
          &w_t::cb_op_original_to_sampling, rir())
        .def("float_asu", &w_t::float_asu, rir())
        .def("continuous_shift_flags", &w_t::continuous_shift_flags, ccr())
        .def("point_distance", &w_t::point_distance, ccr())
        .def("buffer_thickness", &w_t::buffer_thickness, ccr())
        .def("all_twelve_neighbors", &w_t::all_twelve_neighbors)
        .def("box_lower", &w_t::box_lower, ccr())
        .def("box_upper", &w_t::box_upper, ccr())
        .def("at_end", &w_t::at_end)
        .def("next_site_frac", next_site_frac)
        .def("next", next_site_frac)
        .def("__iter__", scitbx::boost_python::pass_through)
        .def("all_sites_frac", &w_t::all_sites_frac)
        .def("restart", &w_t::restart)
        .def("count_sites", &w_t::count_sites)
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
