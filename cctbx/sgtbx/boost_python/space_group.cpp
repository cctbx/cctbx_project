#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/space_group_type.h>
#include <cctbx/sgtbx/space_group_hash.h>
#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <scitbx/boost_python/utils.h>
#include <boost_adaptbx/hash.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct space_group_wrappers : boost::python::pickle_suite,
                                boost_adaptbx::py_hashable<space_group>
  {
    typedef space_group w_t;

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

    static boost::python::tuple
    getinitargs(w_t const& o)
    {
      return boost::python::make_tuple(o.type().hall_symbol());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("space_group")
        .def(init<parse_string&, optional<bool, bool, bool, int> >((
          arg("hall_symbol"),
          arg("pedantic")=false,
          arg("no_centring_type_symbol")=false,
          arg("no_expand")=false,
          arg("t_den")=sg_t_den)))
        .def(init<std::string const&, optional<bool, bool, bool, int> >((
          arg("hall_symbol"),
          arg("pedantic")=false,
          arg("no_centring_type_symbol")=false,
          arg("no_expand")=false,
          arg("t_den")=sg_t_den)))
        .def(init<space_group_symbols const&, optional<int> >((
          arg("space_group_symbols"),
          arg("t_den")=sg_t_den)))
        .def(init<space_group const&>((arg("other"))))
        .def("reset", &w_t::reset, (arg("t_den")=sg_t_den))
        .def("expand_ltr", &w_t::expand_ltr, return_self<>(), (arg("new_t")))
        .def("expand_inv", &w_t::expand_inv,
          return_self<>(), (arg("new_inv_t")))
        .def("expand_smx", (space_group&(w_t::*)(rt_mx const&))
          &w_t::expand_smx, return_self<>(), (arg("new_smx")))
        .def("expand_smx", (space_group&(w_t::*)(std::string const&))
          &w_t::expand_smx, return_self<>(), (arg("smx_symbol")))
        .def("expand_conventional_centring_type",
          &w_t::expand_conventional_centring_type, (arg("symbol")))
        .def("parse_hall_symbol", &w_t::parse_hall_symbol, (
          arg("hall_symbol"),
          arg("pedantic")=false,
          arg("no_centring_type_symbol")=false))
        .def("change_basis", &w_t::change_basis, (arg("cb_op")))
        .def("change_of_origin_realising_origin_centricity",
             &w_t::change_of_origin_realising_origin_centricity)
        .def("r_den", &w_t::r_den)
        .def("t_den", &w_t::t_den)
        .def("order_p", &w_t::order_p)
        .def("order_z", &w_t::order_z)
        .def("__len__", &w_t::order_z)
        .def("n_equivalent_positions", &w_t::n_equivalent_positions)
        .def("n_ltr", &w_t::n_ltr)
        .def("is_centric", (bool(w_t::*)() const) &w_t::is_centric)
        .def("is_origin_centric", &w_t::is_origin_centric)
        .def("f_inv", &w_t::f_inv)
        .def("n_smx", &w_t::n_smx)
        .def("__call__", call_3, (
          arg("i_ltr"),
          arg("i_inv"),
          arg("i_smx")))
        .def("__call__", getitem, (arg("i_op")))
        .def("__getitem__", getitem)
        .def("make_tidy", &w_t::make_tidy, return_self<>())
        .def("is_tidy", &w_t::is_tidy)
        .def("contains", &w_t::contains, (arg("smx")))
        .def("__contains__", &w_t::contains)
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("__hash__", py_hash)
        .def("conventional_centring_type_symbol",
          &w_t::conventional_centring_type_symbol)
        .def("z2p_op", &w_t::z2p_op, (
          arg("r_den")=cb_r_den,
          arg("t_den")=cb_t_den))
        .def("construct_z2p_op", &w_t::construct_z2p_op, (
          arg("r_den")=cb_r_den,
          arg("t_den")=cb_t_den))
        .def("is_chiral", &w_t::is_chiral)
        .def("is_sys_absent",
          (bool(w_t::*)(miller::index<> const&) const)
          &w_t::is_sys_absent, (arg("miller_index")))
        .def("is_sys_absent",
          (af::shared<bool>(w_t::*)
            (af::const_ref<miller::index<> > const&) const)
          &w_t::is_sys_absent, (arg("miller_indices")))
        .def("is_centric",
          (bool(w_t::*)(miller::index<> const&) const)
          &w_t::is_centric, (arg("miller_index")))
        .def("is_centric",
          (af::shared<bool>(w_t::*)
            (af::const_ref<miller::index<> > const&) const)
          &w_t::is_centric, (arg("miller_indices")))
        .def("phase_restriction", &w_t::phase_restriction, (
          arg("miller_index")))
        .def("is_valid_phase", &w_t::is_valid_phase, (
          arg("miller_index"),
          arg("phi"),
          arg("deg")=false,
          arg("tolerance")=1e-5))
        .def("nearest_valid_phases", &w_t::nearest_valid_phases, (
          arg("miller_indices"),
          arg("phases"),
          arg("deg")=false))
        .def("multiplicity",
          (int(w_t::*)(miller::index<> const&, bool) const)
          &w_t::multiplicity, (arg("miller_index")))
        .def("multiplicity",
          (af::shared<int>(w_t::*)
            (af::const_ref<miller::index<> > const&, bool) const)
          &w_t::multiplicity, (arg("miller_indices")))
        .def("epsilon",
          (int(w_t::*)(miller::index<> const&) const)
          &w_t::epsilon, (arg("miller_index")))
        .def("epsilon",
          (af::shared<int>(w_t::*)
            (af::const_ref<miller::index<> > const&) const)
          &w_t::epsilon, (arg("miller_indices")))
        .def("multiplicity",
          (int(w_t::*)(vec3_rat const&) const)
            &w_t::multiplicity,
              (arg("site")))
        .def("average_unit_cell", &w_t::average_unit_cell, (
          arg("unit_cell")))
        .def("is_compatible_unit_cell", &w_t::is_compatible_unit_cell, (
          arg("unit_cell"),
          arg("relative_length_tolerance")=0.01,
          arg("absolute_angle_tolerance")=1.))
        .def("average_u_star",
          (scitbx::sym_mat3<double>(w_t::*)(
            scitbx::sym_mat3<double> const&) const) &w_t::average_u_star, (
          arg("u_star")))
        .def("build_derived_acentric_group",
          &w_t::build_derived_acentric_group)
        .def("build_derived_group",
           &w_t::build_derived_group)
        .def("build_derived_reflection_intensity_group",
          &w_t::build_derived_reflection_intensity_group, (
           arg("anomalous_flag")))
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
          &w_t::refine_gridding, (arg("grid")))
        .def("all_ops", &w_t::all_ops, (arg("mod")=0, arg("cancel")=false))
        .def("unique", &w_t::unique, (arg("special_op")))
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
