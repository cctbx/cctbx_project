/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/sgtbx/site_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct site_symmetry_wrappers
  {
    typedef site_symmetry w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_compatible_u_star_overloads, is_compatible_u_star, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("site_symmetry", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  fractional<double> const&,
                  optional<double, bool> >())
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("space_group", &w_t::space_group, rir())
        .def("original_site", &w_t::original_site, ccr())
        .def("min_distance_sym_equiv", &w_t::min_distance_sym_equiv)
        .def("exact_site", &w_t::exact_site, ccr())
        .def("distance_moved", &w_t::distance_moved)
        .def("shortest_distance", &w_t::shortest_distance)
        .def("check_min_distance_sym_equiv",&w_t::check_min_distance_sym_equiv)
        .def("multiplicity", &w_t::multiplicity)
        .def("special_op", &w_t::special_op, ccr())
        .def("is_point_group_1", &w_t::is_point_group_1)
        .def("point_group_type", &w_t::point_group_type)
        .def("unique_ops", &w_t::unique_ops)
        .def("is_compatible_u_star",
           (bool(w_t::*)(scitbx::sym_mat3<double> const&, double) const)0,
           is_compatible_u_star_overloads())
        .def("average_u_star",
          (scitbx::sym_mat3<double>
            (w_t::*)(scitbx::sym_mat3<double> const&) const)
          &w_t::average_u_star)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_site_symmetry()
  {
    site_symmetry_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
