#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/adptbx.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace adptbx { namespace boost_python {

  void wrap_hirshfeld();

namespace {

  sym_mat3<double>
  debye_waller_factor_u_star_gradient_coefficients_double(
    miller::index<> const& h)
  {
    return debye_waller_factor_u_star_gradient_coefficients(
      h, scitbx::type_holder<double>());
  }

  af::shared<double>
  debye_waller_factor_u_star_curvature_coefficients_double(
    miller::index<> const& h)
  {
    return debye_waller_factor_u_star_curvature_coefficients(
      h, scitbx::type_holder<double>());
  }

  struct eigensystem_wrappers
  {
    typedef eigensystem<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("eigensystem", no_init)
        .def(init<sym_mat3<double> const&>((arg("sym_mat3"))))
        .def("vectors", &w_t::vectors, ccr())
        .def("values", &w_t::values, ccr())
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;

    def("u_as_b", (double(*)(double)) u_as_b);
    def("b_as_u", (double(*)(double)) b_as_u);
    def("u_as_b", (sym_mat3<double>(*)(sym_mat3<double> const&)) u_as_b);
    def("b_as_u", (sym_mat3<double>(*)(sym_mat3<double> const&)) b_as_u);

#define CCTBX_DEF(function_name) \
    def(# function_name, \
      (sym_mat3<double>(*)(uctbx::unit_cell const&, sym_mat3<double> const&)) \
      function_name);
    CCTBX_DEF(u_cif_as_u_star)
    CCTBX_DEF(u_star_as_u_cif)
    CCTBX_DEF(u_cart_as_u_star)
    CCTBX_DEF(u_star_as_u_cart)
    CCTBX_DEF(u_cart_as_u_cif)
    CCTBX_DEF(u_cif_as_u_cart)
    def("u_star_as_beta",
      (sym_mat3<double>(*)(sym_mat3<double> const&)) u_star_as_beta);
    def("beta_as_u_star",
      (sym_mat3<double>(*)(sym_mat3<double> const&)) beta_as_u_star);
    CCTBX_DEF(u_cart_as_beta)
    CCTBX_DEF(beta_as_u_cart)
    CCTBX_DEF(u_cif_as_beta)
    CCTBX_DEF(beta_as_u_cif)
#undef CCTBX_DEF

    def("u_cart_as_u_iso",
      (double(*)(sym_mat3<double> const&)) u_cart_as_u_iso);
    def("u_iso_as_u_cart",
      (sym_mat3<double>(*)(double const&)) u_iso_as_u_cart);

#define CCTBX_DEF(function_name) \
    def(# function_name, \
      (double(*)(uctbx::unit_cell const&, sym_mat3<double> const&)) \
      function_name);
    CCTBX_DEF(u_star_as_u_iso)
    CCTBX_DEF(u_cif_as_u_iso)
    CCTBX_DEF(beta_as_u_iso)
#undef CCTBX_DEF

#define CCTBX_DEF(function_name) \
    def(# function_name, \
      (sym_mat3<double>(*)(uctbx::unit_cell const&, double const&)) \
      function_name);
    CCTBX_DEF(u_iso_as_u_star)
    CCTBX_DEF(u_iso_as_u_cif)
    CCTBX_DEF(u_iso_as_beta)
#undef CCTBX_DEF

    class_<factor_u_cart_u_iso<> >("factor_u_cart_u_iso", no_init)
      .def(init<sym_mat3<double> const&>((arg("u_cart"))))
      .def_readonly("u_iso", &factor_u_cart_u_iso<>::u_iso)
      .add_property("u_cart_minus_u_iso",
        make_getter(&factor_u_cart_u_iso<>::u_cart_minus_u_iso, rbv()))
    ;
    class_<factor_u_star_u_iso<> >("factor_u_star_u_iso", no_init)
      .def(init<
        uctbx::unit_cell const&, sym_mat3<double> const&>((
          arg("unit_cell"), arg("u_star"))))
      .def_readonly("u_iso", &factor_u_star_u_iso<>::u_iso)
      .add_property("u_star_minus_u_iso",
        make_getter(&factor_u_star_u_iso<>::u_star_minus_u_iso, rbv()))
    ;
    class_<factor_u_cif_u_iso<> >("factor_u_cif_u_iso", no_init)
      .def(init<
        uctbx::unit_cell const&, sym_mat3<double> const&>((
          arg("unit_cell"), arg("u_cif"))))
      .def_readonly("u_iso", &factor_u_cif_u_iso<>::u_iso)
      .add_property("u_cif_minus_u_iso",
        make_getter(&factor_u_cif_u_iso<>::u_cif_minus_u_iso, rbv()))
    ;
    class_<factor_beta_u_iso<> >("factor_beta_u_iso", no_init)
      .def(init<
        uctbx::unit_cell const&, sym_mat3<double> const&>((
          arg("unit_cell"), arg("beta"))))
      .def_readonly("u_iso", &factor_beta_u_iso<>::u_iso)
      .add_property("beta_minus_u_iso",
        make_getter(&factor_beta_u_iso<>::beta_minus_u_iso, rbv()))
    ;

    def("debye_waller_factor_b_iso",
      (double(*)(
        double,
        double, double, bool))
          debye_waller_factor_b_iso, (
            arg("stol_sq"),
            arg("b_iso"),
            arg("exp_arg_limit")=50,
            arg("truncate_exp_arg")=false));
    def("debye_waller_factor_b_iso",
      (af::shared<double>(*)(
        af::const_ref<double> const&,
        double const&, double const&, bool))
          debye_waller_factor_b_iso, (
            arg("stol_sq"),
            arg("b_iso"),
            arg("exp_arg_limit")=50,
            arg("truncate_exp_arg")=false));
    def("debye_waller_factor_u_iso",
      (double(*)(double, double)) debye_waller_factor_u_iso);
    def("debye_waller_factor_b_iso",
      (double(*)(uctbx::unit_cell const&, miller::index<> const&, double))
      debye_waller_factor_b_iso);
    def("debye_waller_factor_u_iso",
      (double(*)(uctbx::unit_cell const&, miller::index<> const&, double))
      debye_waller_factor_u_iso);
    def("debye_waller_factor_u_star",
      (double(*)(
        miller::index<> const&,
        sym_mat3<double> const&, double const&, bool))
          debye_waller_factor_u_star, (
            arg("h"),
            arg("u_star"),
            arg("exp_arg_limit")=50,
            arg("truncate_exp_arg")=false));
    def("debye_waller_factor_u_star",
      (af::shared<double>(*)(
        af::const_ref<miller::index<> > const&,
        sym_mat3<double> const&, double const&, bool))
          debye_waller_factor_u_star, (
            arg("miller_indices"),
            arg("u_star"),
            arg("exp_arg_limit")=50,
            arg("truncate_exp_arg")=false));
    def("debye_waller_factor_u_star_gradient_coefficients",
      (scitbx::sym_mat3<double>(*)(
        miller::index<> const&))
      debye_waller_factor_u_star_gradient_coefficients_double);
    def("debye_waller_factor_u_star_gradient_coefficients",
      (af::shared<scitbx::sym_mat3<double> >(*)(
        af::const_ref<miller::index<> > const&))
      debye_waller_factor_u_star_gradient_coefficients);
    def("debye_waller_factor_u_star_curvature_coefficients",
      debye_waller_factor_u_star_curvature_coefficients_double, (
        arg("h")));
    def("debye_waller_factor_beta",
      (double(*)(miller::index<> const&, sym_mat3<double> const&))
      debye_waller_factor_beta);
    def("debye_waller_factor_u_cif",
      (double(*)(uctbx::unit_cell const&,
                 miller::index<> const&,
                 sym_mat3<double> const&))
      debye_waller_factor_u_cif);
    def("debye_waller_factor_u_cart",
      (double(*)(uctbx::unit_cell const&,
                 miller::index<> const&,
                 sym_mat3<double> const&))
      debye_waller_factor_u_cart);

    def("grad_u_star_as_u_cart",
      (sym_mat3<double>(*)(
        uctbx::unit_cell const&,
        sym_mat3<double> const&)) grad_u_star_as_u_cart);
    def("grad_u_star_as_u_cart",
      (af::shared<sym_mat3<double> >(*)(
        uctbx::unit_cell const&,
        af::const_ref<sym_mat3<double> > const&)) grad_u_star_as_u_cart);
    def("grad_u_cart_as_u_star",
      (sym_mat3<double>(*)(
        uctbx::unit_cell const&,
        sym_mat3<double> const&)) grad_u_cart_as_u_star);
    def("grad_u_cart_as_u_star",
      (af::shared<sym_mat3<double> >(*)(
        uctbx::unit_cell const&,
        af::const_ref<sym_mat3<double> > const&)) grad_u_cart_as_u_star);

    def("eigenvalues", (vec3<double>(*)(sym_mat3<double> const&)) eigenvalues);
    def("is_positive_definite",
      (bool(*)(vec3<double> const&)) is_positive_definite);
    def("is_positive_definite",
      (bool(*)(vec3<double> const&, double const&)) is_positive_definite);
    def("is_positive_definite",
      (bool(*)(sym_mat3<double> const&)) is_positive_definite);
    def("is_positive_definite",
      (bool(*)(sym_mat3<double> const&, double const&)) is_positive_definite);

    def("is_positive_definite",
      (af::shared<bool>(*)(af::const_ref<sym_mat3<double> > const&,
                                         double const&)) is_positive_definite);

    def("eigenvalue_filtering",
      (sym_mat3<double>(*)(
        sym_mat3<double> const&, double const&, double const&))
          eigenvalue_filtering, (
            arg("u_cart"),
            arg("u_min")=0,
            arg("u_max")=0));

    eigensystem_wrappers::wrap();

    def("c_u_c_transpose",
      (sym_mat3<double>(*)(mat3<double> const&, sym_mat3<double> const&))
      c_u_c_transpose);

    wrap_hirshfeld();
}

} // namespace <anonymous>
}}} // namespace cctbx::adptbx::boost_python

BOOST_PYTHON_MODULE(cctbx_adptbx_ext)
{
  cctbx::adptbx::boost_python::init_module();
}
