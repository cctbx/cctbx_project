#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
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
        .def(init<scitbx::vec3<double> const&, double>())
        .add_property("n",
          make_getter(&w_t::n, rbv()),
          make_setter(&w_t::n, dcp()))
        .def_readwrite("c", &w_t::c)
        .def("evaluate", &w_t::evaluate)
        .def("is_inside", &w_t::is_inside)
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
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("direct_space_asu_float_asu", no_init)
        .def(init<w_t::facets_t const&>())
        .add_property("facets", make_getter(&w_t::facets, rbv()))
        .def("is_inside", &w_t::is_inside)
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
