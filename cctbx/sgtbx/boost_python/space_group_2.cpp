/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/space_group_type.h>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct space_group_wrappers
  {
    typedef space_group w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_valid_phase_overloads, is_valid_phase, 2, 4)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_compatible_unit_cell_overloads, is_compatible_unit_cell, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      all_ops_overloads, all_ops, 0, 2)

    static void
    wrap_more(boost::python::class_<w_t>& cl)
    {
      using namespace boost::python;
      cl
        .def("is_sys_absent",
          (bool(w_t::*)(miller::index<> const&) const)
          &w_t::is_sys_absent)
        .def("is_sys_absent",
          (bool(w_t::*)(miller::index<> const&) const)
          &w_t::is_sys_absent)
        .def("is_sys_absent",
          (af::shared<bool>(w_t::*)
            (af::const_ref<miller::index<> > const&) const)
          &w_t::is_sys_absent)
        .def("is_centric",
          (bool(w_t::*)(miller::index<> const&) const)
          &w_t::is_centric)
        .def("is_centric",
          (af::shared<bool>(w_t::*)
            (af::const_ref<miller::index<> > const&) const)
          &w_t::is_centric)
        .def("phase_restriction", &w_t::phase_restriction)
        .def("is_valid_phase",
          &w_t::is_valid_phase,
          is_valid_phase_overloads())
        .def("multiplicity",
          (int(w_t::*)(miller::index<> const&, bool) const)
          &w_t::multiplicity)
        .def("multiplicity",
          (af::shared<int>(w_t::*)
            (af::const_ref<miller::index<> > const&, bool) const)
          &w_t::multiplicity)
        .def("epsilon",
          (int(w_t::*)(miller::index<> const&) const)
          &w_t::epsilon)
        .def("epsilon",
          (af::shared<int>(w_t::*)
            (af::const_ref<miller::index<> > const&) const)
          &w_t::epsilon)
        .def("is_compatible_unit_cell",
          &w_t::is_compatible_unit_cell,
          is_compatible_unit_cell_overloads())
        .def("build_derived_patterson_group",
          &w_t::build_derived_patterson_group)
        .def("build_derived_point_group",
          &w_t::build_derived_point_group)
        .def("build_derived_laue_group",
          &w_t::build_derived_laue_group)
        .def("point_group_type", &w_t::point_group_type)
        .def("laue_group_type", &w_t::laue_group_type)
        .def("crystal_system", &w_t::crystal_system)
        .def("match_tabulated_settings", &w_t::match_tabulated_settings)
        .def("gridding", &w_t::gridding)
        .def("refine_gridding",
          (sg_vec3(w_t::*)(sg_vec3 const&) const)
          &w_t::refine_gridding)
        .def("all_ops", &w_t::all_ops, all_ops_overloads())
        .def("unique", &w_t::unique)
        .def("type", &w_t::type)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_space_group_2(boost::python::class_<space_group>& cl)
  {
    space_group_wrappers::wrap_more(cl);
  }

}}} // namespace cctbx::sgtbx::boost_python
