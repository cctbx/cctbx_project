/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Refactored (rwgk)
     2001 Oct: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/adptbx.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace adptbx { namespace boost_python {
namespace {

  struct eigensystem_wrappers
  {
    typedef eigensystem<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("eigensystem", no_init)
        .def(init<sym_mat3<double> const&, optional<double> >())
        .def("vectors", &w_t::vectors, ccr())
        .def("values", &w_t::values, ccr())
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

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

    def("debye_waller_factor_b_iso",
      (double(*)(double, double)) debye_waller_factor_b_iso);
    def("debye_waller_factor_u_iso",
      (double(*)(double, double)) debye_waller_factor_u_iso);
    def("debye_waller_factor_b_iso",
      (double(*)(uctbx::unit_cell const&, miller::index<> const&, double))
      debye_waller_factor_b_iso);
    def("debye_waller_factor_u_iso",
      (double(*)(uctbx::unit_cell const&, miller::index<> const&, double))
      debye_waller_factor_u_iso);
    def("debye_waller_factor_u_star",
      (double(*)(miller::index<> const&, sym_mat3<double> const&))
      debye_waller_factor_u_star);
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
      (sym_mat3<double>(*)(uctbx::unit_cell const&,
                           sym_mat3<double> const&)) grad_u_star_as_u_cart);
    def("grad_u_cart_as_u_star",
      (sym_mat3<double>(*)(uctbx::unit_cell const&,
                           sym_mat3<double> const&)) grad_u_cart_as_u_star);

    def("eigenvalues", (vec3<double>(*)(sym_mat3<double> const&)) eigenvalues);
    def("is_positive_definite",
      (bool(*)(vec3<double> const&)) is_positive_definite);
    def("is_positive_definite",
      (bool(*)(sym_mat3<double> const&)) is_positive_definite);
    def("eigenvalue_filtering",
      (sym_mat3<double>(*)(sym_mat3<double> const&)) eigenvalue_filtering);

    eigensystem_wrappers::wrap();

    def("c_u_c_transpose",
      (sym_mat3<double>(*)(mat3<double> const&, sym_mat3<double> const&))
      c_u_c_transpose);
}

} // namespace <anonymous>
}}} // namespace cctbx::adptbx::boost_python

BOOST_PYTHON_MODULE(adptbx_ext)
{
  cctbx::adptbx::boost_python::init_module();
}
