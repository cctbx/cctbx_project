/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct sym_equiv_sites_wrappers
  {
    typedef sym_equiv_sites<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("sym_equiv_sites", no_init)
        .def(init<site_symmetry const&>())
        .def(init<wyckoff::mapping const&>())
        .def(init<sgtbx::space_group const&,
                  fractional<> const&,
                  optional<uctbx::unit_cell const&> >())
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  fractional<> const&,
                  rt_mx const&>())
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  fractional<> const&,
                  optional<double, double> >())
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("space_group", &w_t::space_group, rir())
        .def("original_site", &w_t::original_site, ccr())
        .def("special_op", &w_t::special_op, ccr())
        .def("max_accepted_tolerance", &w_t::max_accepted_tolerance)
        .def("coordinates", &w_t::coordinates)
        .def("sym_op", &w_t::sym_op)
        .def("sym_op_indices", &w_t::sym_op_indices)
      ;
    }
  };

  struct min_sym_equiv_distance_info_wrappers
  {
    typedef min_sym_equiv_distance_info<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("min_sym_equiv_distance_info", no_init)
        .def(init<sym_equiv_sites<> const&,
                  fractional<> const&,
                  optional<af::tiny<bool, 3> const&> >())
        .def("sym_op", &w_t::sym_op, ccr())
        .def("continuous_shifts", &w_t::continuous_shifts, ccr())
        .def("diff", &w_t::diff, ccr())
        .def("dist", &w_t::dist)
        .def("apply", &w_t::apply)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_sym_equiv_sites()
  {
    sym_equiv_sites_wrappers::wrap();
    min_sym_equiv_distance_info_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
