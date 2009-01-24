#include <cctbx/boost_python/flex_fwd.h>

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

    class_<w_t>("cut", no_init)
      .def(init<
        ivector3_t const&,
        rational_t ,
        optional< bool > >((
          arg_("n"),
          arg_("c"),
          arg_("inclusive") )))
      .def_readonly("n", &w_t::n)
      .def_readonly("c", &w_t::c)
      .def_readonly("inclusive", &w_t::inclusive)
      .def("__pos__", &w_t::operator+)
      .def("__neg__", &w_t::operator-)
      .def("__inv__", &w_t::operator~)
      .def("__mul__", &w_t::operator*)
      .def("__div__", &w_t::operator/)
      // .def("is_inside", &w_t::is_inside)
      .def("get_point_in_plane", &w_t::get_point_in_plane)
      .def("change_basis", &w_t::change_basis)
      .def("evaluate", &w_t::evaluate)
    ;
  }

  void wrap_direct_space_asu()
  {
    typedef direct_space_asu w_t;

    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;

    class_<w_t>("direct_space_asu", no_init)
      .def(init< unsigned short >(( arg_("spgr") )))
      .def_readonly("hall_symbol", &w_t::hall_symbol)
      .def("is_inside", &w_t::is_inside)
      .def("change_basis", &w_t::change_basis)
      .def("get_nth_plane", &w_t::get_nth_plane)
      .def("volume_only", &w_t::volume_only)
      .def("in_which_planes", &w_t::in_which_planes)
      .def("n_faces", &w_t::n_faces)
      .def("volume_vertices", &w_t::volume_vertices)
      .def("box_corners", &w_t::box_corners)
    ;
  }


  void init_module()
  {
    wrap_cut();
    wrap_direct_space_asu();
  }

}}}} // namespace cctbx::sgtbx::asu

BOOST_PYTHON_MODULE(cctbx_sgtbx_asu_ext)
{
  cctbx::sgtbx::asu::init_module();
}
