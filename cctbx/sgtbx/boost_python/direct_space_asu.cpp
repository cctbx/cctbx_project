#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/boost_python/container_conversions.h>
#include <cctbx/sgtbx/direct_space_asu.h>

namespace cctbx { namespace sgtbx { namespace direct_space_asu {

namespace {

  struct float_cut_plane_wrappers
  {
    typedef float_cut_plane<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("direct_space_asu_float_cut_plane", no_init)
        .def(init<fractional<double> const&, double>((arg("n"), arg("c"))))
        .add_property("n",
          make_getter(&w_t::n, rbv()),
          make_setter(&w_t::n, dcp()))
        .def_readwrite("c", &w_t::c)
        .def("evaluate", &w_t::evaluate)
        .def("is_inside", &w_t::is_inside)
        .def("get_point_in_plane", &w_t::get_point_in_plane)
        .def("add_buffer", &w_t::add_buffer,
          (arg("unit_cell"), arg("thickness")))
      ;
    }
  };

  struct float_asu_wrappers
  {
    typedef float_asu<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("direct_space_asu_float_asu", no_init)
        .def(init<uctbx::unit_cell const&, w_t::facets_t const&>(
          (arg("unit_cell"), arg("facets"))))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("facets", &w_t::facets, ccr())
        .def("is_inside", &w_t::is_inside)
        .def("_add_buffer", &w_t::add_buffer)
      ;
    }
  };

  void register_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping<float_asu_wrappers::w_t::facets_t, fixed_capacity_policy>();
  }

}} // namespace direct_space_asu::<anoymous>

namespace boost_python {

  void wrap_direct_space_asu()
  {
    direct_space_asu::float_cut_plane_wrappers::wrap();
    direct_space_asu::float_asu_wrappers::wrap();
    direct_space_asu::register_tuple_mappings();
  }

}}} // namespace cctbx::sgtbx::boost_python
