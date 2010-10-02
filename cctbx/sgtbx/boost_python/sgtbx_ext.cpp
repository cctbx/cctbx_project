#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/site_symmetry.h>
#include <cctbx/sgtbx/operator_from_axis_direction.h>
#include <scitbx/array_family/shared.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <scitbx/boost_python/container_conversions.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  void wrap_brick();
  void wrap_change_of_basis_op();
  void wrap_find_affine();
  void wrap_lattice_symmetry();
  void wrap_phase_info();
  void wrap_reciprocal_space_asu();
  void wrap_rot_mx();
  void wrap_rt_mx();
  void wrap_search_symmetry();
  void wrap_seminvariant();
  void wrap_site_symmetry();
  void wrap_space_group();
  void wrap_space_group_type();
  void wrap_sym_equiv_sites();
  void wrap_symbols();
  void wrap_tensor_rank_2();
  void wrap_tr_vec();
  void wrap_wyckoff();
  void wrap_select_generators();

namespace {

  struct parse_string_wrappers
  {
    typedef parse_string w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("parse_string", no_init)
        .def(init<std::string const&>((arg("str"))))
        .def("string", &w_t::string)
        .def("where", &w_t::where)
      ;
    }
  };

  struct crystal_system_code_to_string
  {
    static PyObject* convert(crystal_system::code const& c)
    {
      using namespace boost::python;
      return incref(object(crystal_system::label(c)).ptr());
    }
  };

  struct matrix_group_code_to_string
  {
    static PyObject* convert(matrix_group::code const& c)
    {
      using namespace boost::python;
      return incref(object(c.label()).ptr());
    }
  };

  fractional<>
  fractional_mod_positive(
    fractional<> const& site)
  {
    return site.mod_positive();
  }

  fractional<>
  fractional_mod_short(
    fractional<> const& site)
  {
    return site.mod_short();
  }

  void register_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping_variable_capacity<af::shared<site_symmetry_ops> >();
    tuple_mapping_fixed_capacity<af::small<ss_vec_mod, 3> >();
  }

  void init_module()
  {
    using namespace boost::python;

    sanity_check();

    scope s;
    s.attr("sg_t_den") = sg_t_den;
    s.attr("cb_r_den") = cb_r_den;
    s.attr("cb_t_den") = cb_t_den;

    register_tuple_mappings();

    parse_string_wrappers::wrap();

    to_python_converter<crystal_system::code, crystal_system_code_to_string>();
    to_python_converter<matrix_group::code, matrix_group_code_to_string>();

    def("fractional_mod_positive", fractional_mod_positive, (arg("site")));
    def("fractional_mod_short", fractional_mod_short, (arg("site")));

    wrap_brick();
    wrap_change_of_basis_op();
    wrap_find_affine();
    wrap_lattice_symmetry();
    wrap_phase_info();
    wrap_reciprocal_space_asu();
    wrap_rot_mx();
    wrap_rt_mx();
    wrap_search_symmetry();
    wrap_seminvariant();
    wrap_site_symmetry();
    wrap_space_group();
    wrap_space_group_type();
    wrap_sym_equiv_sites();
    wrap_symbols();
    wrap_tensor_rank_2();
    wrap_tr_vec();
    wrap_wyckoff();
    wrap_select_generators();

    def("n_fold_operator_from_axis_direction",
      n_fold_operator_from_axis_direction, (
        arg("ev_cart"), arg("n"), arg("sense")=1));
  }

} // namespace <anonymous>
}}} // namespace cctbx::sgtbx::boost_python

BOOST_PYTHON_MODULE(cctbx_sgtbx_ext)
{
  cctbx::sgtbx::boost_python::init_module();
}
