/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/space_group_type.h>
#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/boost_python/utils.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct space_group_wrappers : boost::python::pickle_suite
  {
    typedef space_group w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      reset_overloads, reset, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      parse_hall_symbol_overloads, parse_hall_symbol, 1, 3)

    static rt_mx
    getitem(w_t const& o, std::size_t i_op)
    {
      if (i_op >= o.order_z()) scitbx::boost_python::raise_index_error();
      return o(i_op);
    }

    static rt_mx
    call_3(w_t const& o,
           std::size_t i_ltr, std::size_t i_inv, std::size_t i_smx)
    {
      return o(i_ltr, i_inv, i_smx);
    }

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      z2p_op_overloads, z2p_op, 0, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      construct_z2p_op_overloads, construct_z2p_op, 0, 2)

    static boost::python::tuple
    getinitargs(w_t const& o)
    {
      return boost::python::make_tuple(o.type().hall_symbol());
    }

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_valid_phase_overloads, is_valid_phase, 2, 4)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_compatible_unit_cell_overloads, is_compatible_unit_cell, 1, 3)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      all_ops_overloads, all_ops, 0, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("space_group")
        .def(init<parse_string&, optional<bool, bool, bool> >())
        .def(init<std::string const&, optional<bool, bool, bool> >())
        .def(init<space_group_symbols const&>())
        .def(init<space_group const&>())
        .def("reset", &w_t::reset, reset_overloads())
        .def("expand_ltr", &w_t::expand_ltr)
        .def("expand_inv", &w_t::expand_inv)
        .def("expand_smx", (void(w_t::*)(rt_mx const&)) &w_t::expand_smx)
        .def("expand_conventional_centring_type",
          &w_t::expand_conventional_centring_type)
        .def("parse_hall_symbol",
          &w_t::parse_hall_symbol,
          parse_hall_symbol_overloads())
        .def("change_basis", &w_t::change_basis)
        .def("r_den", &w_t::r_den)
        .def("t_den", &w_t::t_den)
        .def("order_p", &w_t::order_p)
        .def("order_z", &w_t::order_z)
        .def("__len__", &w_t::order_z)
        .def("n_ltr", &w_t::n_ltr)
        .def("is_centric", (bool(w_t::*)() const) &w_t::is_centric)
        .def("is_origin_centric", &w_t::is_origin_centric)
        .def("f_inv", &w_t::f_inv)
        .def("n_smx", &w_t::n_smx)
        .def("__call__", call_3)
        .def("__call__", getitem)
        .def("__getitem__", getitem)
        .def("make_tidy", &w_t::make_tidy)
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("conventional_centring_type_symbol",
          &w_t::conventional_centring_type_symbol)
        .def("z2p_op", &w_t::z2p_op, z2p_op_overloads())
        .def("construct_z2p_op",
          &w_t::construct_z2p_op,
          construct_z2p_op_overloads())
        .def("is_chiral", &w_t::is_chiral)
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
        .def("average_unit_cell", &w_t::average_unit_cell)
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
        .def_pickle(space_group_wrappers())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_space_group()
  {
    space_group_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
