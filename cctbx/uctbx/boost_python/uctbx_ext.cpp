/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Ported from cctbx/uxtbx/uctbxmodule.cpp (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/uctbx.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace uctbx { namespace boost_python {

  void wrap_fast_minimal_reduction();

namespace {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    d_star_sq_as_two_theta_overloads, d_star_sq_as_two_theta, 2, 3)

  struct unit_cell_wrappers : boost::python::pickle_suite
  {
    typedef unit_cell w_t;
    typedef cartesian<> cart_t;
    typedef fractional<> frac_t;
    typedef miller::index<> mix_t;
    typedef af::const_ref<mix_t> cr_mix_t;
    typedef af::shared<double> sh_dbl_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_similar_to_overloads, is_similar_to, 1, 3)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      change_basis_overloads, change_basis, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      max_miller_indices_overloads, max_miller_indices, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      two_theta_overloads, two_theta, 2, 3)

    static boost::python::tuple
    getinitargs(w_t const& ucell)
    {
      return boost::python::make_tuple(ucell.parameters());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("unit_cell", no_init)
        .def(init<af::small<double, 6> const&, optional<bool> >())
        .def("parameters", &w_t::parameters, ccr())
        .def("reciprocal_parameters", &w_t::reciprocal_parameters, ccr())
        .def("metrical_matrix", &w_t::metrical_matrix, ccr())
        .def("reciprocal_metrical_matrix",
          &w_t::reciprocal_metrical_matrix, ccr())
        .def("volume", &w_t::volume)
        .def("reciprocal", &w_t::reciprocal)
        .def("longest_vector_sq", &w_t::longest_vector_sq)
        .def("is_similar_to",
          &w_t::is_similar_to,
          is_similar_to_overloads())
        .def("fractionalization_matrix", &w_t::fractionalization_matrix, ccr())
        .def("orthogonalization_matrix", &w_t::orthogonalization_matrix, ccr())
        .def("fractionalize",
          (frac_t(w_t::*)(cart_t const&) const)
          &w_t::fractionalize)
        .def("orthogonalize",
          (cart_t(w_t::*)(frac_t const& xf) const)
          &w_t::orthogonalize)
        .def("length",
          (double(w_t::*)(frac_t const& xf) const)
          &w_t::length)
        .def("distance",
          (double(w_t::*)(frac_t const& xf, frac_t const& yf) const)
          &w_t::distance)
        .def("mod_short_length",
          (double(w_t::*)(frac_t const& xf) const)
          &w_t::mod_short_length)
        .def("mod_short_distance",
          (double(w_t::*)(frac_t const& xf, frac_t const& yf) const)
          &w_t::mod_short_distance)
        .def("min_mod_short_distance",
          (double(w_t::*)
            (af::const_ref<scitbx::vec3<double> > const& xf,
             frac_t const& yf) const)
          &w_t::min_mod_short_distance)
        .def("change_basis",
          (w_t(w_t::*)(uc_mat3 const&, double) const) 0,
          change_basis_overloads())
        .def("max_miller_indices",
          (mix_t(w_t::*)(double, double) const) 0,
          max_miller_indices_overloads())
        .def("d_star_sq",
          (double(w_t::*)(mix_t const&) const)
          &w_t::d_star_sq)
        .def("d_star_sq",
          (sh_dbl_t(w_t::*)(cr_mix_t const&) const)
          &w_t::d_star_sq)
        .def("max_d_star_sq",
          (double(w_t::*)(cr_mix_t const& h) const)
          &w_t::max_d_star_sq)
        .def("min_max_d_star_sq",
          (af::double2(w_t::*)(cr_mix_t const& h) const)
          &w_t::min_max_d_star_sq)
        .def("stol_sq",
          (double(w_t::*)(mix_t  const& h) const)
          &w_t::stol_sq)
        .def("stol_sq",
          (sh_dbl_t(w_t::*)(cr_mix_t const& h) const)
          &w_t::stol_sq)
        .def("two_stol",
          (double(w_t::*)(mix_t const& h) const)
          &w_t::two_stol)
        .def("two_stol",
          (sh_dbl_t(w_t::*)(cr_mix_t const& h) const)
          &w_t::two_stol)
        .def("stol",
          (double(w_t::*)(mix_t const& h) const)
          &w_t::stol)
        .def("stol",
          (sh_dbl_t(w_t::*)(cr_mix_t const& h) const)
          &w_t::stol)
        .def("d",
          (double(w_t::*)(mix_t const& h) const)
          &w_t::d)
        .def("d",
          (sh_dbl_t(w_t::*)(cr_mix_t const& h) const)
          &w_t::d)
        .def("two_theta",
          (double(w_t::*)(mix_t const&, double, bool) const) 0,
          two_theta_overloads())
        .def("two_theta",
          (sh_dbl_t(w_t::*)(cr_mix_t const&, double, bool) const) 0,
          two_theta_overloads())
        .def_pickle(unit_cell_wrappers())
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    def("d_star_sq_as_stol_sq", d_star_sq_as_stol_sq);
    def("d_star_sq_as_two_stol", d_star_sq_as_two_stol);
    def("d_star_sq_as_stol", d_star_sq_as_stol);
    def("d_star_sq_as_two_theta", d_star_sq_as_two_theta,
      d_star_sq_as_two_theta_overloads());
    def("d_star_sq_as_d", d_star_sq_as_d);

    unit_cell_wrappers::wrap();
    wrap_fast_minimal_reduction();
  }

} // namespace <anonymous>
}}} // namespace cctbx::uctbx::boost_python

BOOST_PYTHON_MODULE(uctbx_ext)
{
  cctbx::uctbx::boost_python::init_module();
}
