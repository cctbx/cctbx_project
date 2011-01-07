#include <cctbx/boost_python/flex_fwd.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include "cut.h"
#include "direct_space_asu.h"

namespace cctbx { namespace sgtbx { namespace asu { namespace {

  void wrap_cut()
  {
    typedef cut w_t;

    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;

    bool (w_t::*const cut_is_inside1)( const rvector3_t &) const =
      &w_t::is_inside;
    short (w_t::*const cut_where_is1)( const scitbx::af::int3 &,
      const scitbx::af::int3 & ) const = &w_t::where_is;

    class_<w_t>("cut", no_init)
      .def(init<
        int3_t const&,
        rational_t ,
        optional< bool > >((
          arg("n"),
          arg("c"),
          arg("inclusive") )))
      .def(init<
        int_type, int_type, int_type,
        rational_t ,
        optional< bool > >((
          arg("x"),arg("y"),arg("z"),
          arg("c"),
          arg("inclusive") )))
      .add_property("n", make_getter(&w_t::n, rbv()))
      .def_readonly("c", &w_t::c)
      .def_readonly("inclusive", &w_t::inclusive)
      .def("__pos__", &w_t::operator+)
      .def("__neg__", &w_t::operator-)
      .def("__inv__", &w_t::operator~)
      .def("__mul__", &w_t::operator*)
      .def("__div__", &w_t::operator/)
      .def("__repr__", &w_t::as_string)
      .def("one", &w_t::one)
      .def("is_inside",cut_is_inside1)
      .def("where_is", cut_where_is1)
      .def("get_point_in_plane", &w_t::get_point_in_plane)
      .def("change_basis", &w_t::change_basis)
      .def("evaluate", &w_t::evaluate)
    ;
  }

  direct_space_asu asu_copy(const direct_space_asu &a)
  {
    return direct_space_asu(a);
  }

  scitbx::af::shared<rvector3_t> asu_shape_vertices(const direct_space_asu &a)
  {
    set_rvector3_t v;
    a.shape_vertices(v);
    scitbx::af::shared<rvector3_t> result;
    for(set_rvector3_t::const_iterator it=v.begin(); it!=v.end(); ++it)
      result.push_back(*it);
    return result;
  }

  scitbx::af::long3 asu_grid_limits(const direct_space_asu &a)
  {
    scitbx::af::long3 result;
    a.get_optimized_grid_limits(result);
    return result;
  }

  void wrap_direct_space_asu()
  {
    typedef direct_space_asu w_t;

    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;

    bool (w_t::*const is_inside1)( const rvector3_t &) const = &w_t::is_inside;
    bool (w_t::*const is_inside2)( const rvector3_t &, bool ) const =
      &w_t::is_inside;
    bool (w_t::*const is_inside_shape_only1)(const rvector3_t &point) const =
      &w_t::is_inside_shape_only;
    bool (w_t::*const is_inside_shape_only2)(const scitbx::af::double3 &point,
      double tol) const = &w_t::is_inside_shape_only;
    short (w_t::*const where_is1)(const scitbx::int3 &num,
      const scitbx::int3 &den) const = &w_t::where_is;
    short (w_t::*const where_is2)(const scitbx::int3 &num) const =
      &w_t::where_is;

    class_<w_t>("direct_space_asu", no_init)
      .def(init< const std::string& >(( arg("group_symbol") )))
      .def(init< const space_group_type& >(( arg("group_type") )))
      .def_readonly("hall_symbol", &w_t::hall_symbol)
      .def("is_inside", is_inside1)
      .def("is_inside", is_inside2)
      .def("is_inside_shape_only", is_inside_shape_only1)
      .def("is_inside_shape_only", is_inside_shape_only2)
      .def("where_is", where_is1)
      .def("where_is", where_is2)
      .def("change_basis", &w_t::change_basis)
      .def("get_nth_plane", &w_t::get_nth_plane)
      .def("shape_only", &w_t::shape_only)
      .def("in_which_planes", &w_t::in_which_planes)
      .def("in_which_facets", &w_t::in_which_facets, rbv())
      .def("n_faces", &w_t::n_faces)
      .def("shape_vertices", asu_shape_vertices)
      .def("box_max", &w_t::box_max)
      .def("box_min", &w_t::box_min)
      .def("as_string", &w_t::as_string)
      .def("as_float_asu", &w_t::as_float_asu)
      .def("__copy__", asu_copy)
      .def("optimize_for_grid", &w_t::optimize_for_grid)
      .def("is_optimized", &w_t::is_optimized)
      .def("get_optimized_grid_limits", asu_grid_limits)
    ;
  }

  void init_module()
  {
    wrap_cut();
    wrap_direct_space_asu();
    scitbx::af::boost_python::shared_wrapper< cut >::wrap("cut_shared_array");
    scitbx::af::boost_python::shared_wrapper< rvector3_t >::wrap(
        "rational3_shared_array");
  }

}}}} // namespace cctbx::sgtbx::asu

BOOST_PYTHON_MODULE(cctbx_sgtbx_asu_ext)
{
  cctbx::sgtbx::asu::init_module();
}
