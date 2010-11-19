#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/change_of_basis_op.h>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <scitbx/vec3.h>

namespace cctbx { namespace uctbx { namespace boost_python {

  void wrap_fast_minimum_reduction();

namespace {

  struct unit_cell_wrappers : boost::python::pickle_suite
  {
    typedef unit_cell w_t;
    typedef cartesian<> cart_t;
    typedef fractional<> frac_t;
    typedef miller::index<> mix_t;
    typedef scitbx::vec3<double> frac_mix_t;
    typedef af::const_ref<mix_t> cr_mix_t;
    typedef af::shared<double> sh_dbl_t;

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
      typedef return_internal_reference<> rir;
      class_<w_t>("unit_cell", no_init)
        .def(init<scitbx::mat3<double> const&>((
          arg("orthogonalization_matrix"))))
        .def(init<scitbx::sym_mat3<double> const&>((
          arg("metrical_matrix"))))
        .def(init<af::small<double, 6> const&>((
          arg("parameters"))))
        .def("parameters", &w_t::parameters, ccr())
        .def("reciprocal_parameters", &w_t::reciprocal_parameters, ccr())
        .def("metrical_matrix", &w_t::metrical_matrix, ccr())
        .def("reciprocal_metrical_matrix",
          &w_t::reciprocal_metrical_matrix, ccr())
        .def("volume", &w_t::volume)
        .def("reciprocal", &w_t::reciprocal)
        .def("longest_vector_sq", &w_t::longest_vector_sq)
        .def("shortest_vector_sq", &w_t::shortest_vector_sq)
        .def("is_degenerate",
          &w_t::is_degenerate, (
            arg("min_min_length_over_max_length")=1e-10,
            arg("min_volume_over_min_length")=1e-5))
        .def("is_similar_to",
          &w_t::is_similar_to, (
            arg("other"),
            arg("relative_length_tolerance")=0.01,
            arg("absolute_angle_tolerance")=1.))
        .def("similarity_transformations",
          &w_t::similarity_transformations, (
            arg("other"),
            arg("relative_length_tolerance")=0.02,
            arg("absolute_angle_tolerance")=2,
            arg("unimodular_generator_range")=1))
        .def("fractionalization_matrix", &w_t::fractionalization_matrix, ccr())
        .def("orthogonalization_matrix", &w_t::orthogonalization_matrix, ccr())
        .def("grid_index_as_site_cart_matrix",
          (uc_mat3 (w_t::*)(scitbx::vec3<int> const&) const)
            &w_t::grid_index_as_site_cart_matrix,
              (arg("gridding")))
        .def("fractionalize",
          (scitbx::vec3<double>(w_t::*)(scitbx::vec3<double> const&) const)
          &w_t::fractionalize, (
            arg("site_cart")))
        .def("orthogonalize",
          (scitbx::vec3<double>(w_t::*)(scitbx::vec3<double> const&) const)
          &w_t::orthogonalize, (
            arg("site_frac")))
        .def("fractionalize",
          (af::shared<scitbx::vec3<double> >(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&) const)
              &w_t::fractionalize, (
                arg("sites_cart")))
        .def("orthogonalize",
          (af::shared<scitbx::vec3<double> >(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&) const)
              &w_t::orthogonalize, (
                arg("sites_frac")))
        .def("fractionalize_gradient",
          (scitbx::vec3<double>(w_t::*)(scitbx::vec3<double> const&) const)
          &w_t::fractionalize_gradient, (
            arg("site_cart")))
        .def("u_star_to_u_iso_linear_form",
             &w_t::u_star_to_u_iso_linear_form, ccr())
        .def("u_star_to_u_cif_linear_map",
             &w_t::u_star_to_u_cif_linear_map, ccr())
        .def("u_star_to_u_cart_linear_map",
             &w_t::u_star_to_u_cart_linear_map, ccr())
        .def("length",
          (double(w_t::*)(frac_t const&) const)
          &w_t::length, (
            arg("site_frac")))
        .def("distance",
          (double(w_t::*)(frac_t const&, frac_t const&) const)
          &w_t::distance, (
            arg("site_frac_1"), arg("site_frac_2")))
        .def("angle",
          (boost::optional<double>(w_t::*)(
            frac_t const&, frac_t const&, frac_t const&) const)
          &w_t::angle, (
            arg("site_frac_1"), arg("site_frac_2"), arg("site_frac_3")))
        .def("dihedral",
          (boost::optional<double>(w_t::*)(
            frac_t const&, frac_t const&, frac_t const&, frac_t const&) const)
          &w_t::dihedral, (
            arg("site_frac_1"), arg("site_frac_2"),
            arg("site_frac_3"), arg("site_frac_4")))
        .def("mod_short_length",
          (double(w_t::*)(frac_t const&) const)
          &w_t::mod_short_length, (
            arg("site_frac")))
        .def("mod_short_distance",
          (double(w_t::*)(frac_t const&, frac_t const&) const)
          &w_t::mod_short_distance, (
            arg("site_frac_1"), arg("site_frac_2")))
        .def("min_mod_short_distance",
          (double(w_t::*)
            (af::const_ref<scitbx::vec3<double> > const&,
             frac_t const&) const)
          &w_t::min_mod_short_distance, (
            arg("site_frac_1"), arg("site_frac_2")))
        .def("change_basis",
          (w_t(w_t::*)(uc_mat3 const&, double) const)
            &w_t::change_basis, (
              arg("c_inv_r"), arg("r_den")=1.))
        .def("change_basis",
          (w_t(w_t::*)(sgtbx::change_of_basis_op const&) const)
            &w_t::change_basis, (
              arg("cb_op")))
        .def("max_miller_indices",
          (mix_t(w_t::*)(double, double) const)
            &w_t::max_miller_indices, (
              arg("d_min"), arg("tolerance")=1e-4))
        .def("d_star_sq",
          (double(w_t::*)(mix_t const&) const)
          &w_t::d_star_sq, (
            arg("miller_index")))
        .def("d_star_sq",
          (sh_dbl_t(w_t::*)(cr_mix_t const&) const)
          &w_t::d_star_sq, (
            arg("miller_indices")))
        .def("max_d_star_sq",
          (double(w_t::*)(cr_mix_t const&) const)
          &w_t::max_d_star_sq, (
            arg("miller_indices")))
        .def("min_max_d_star_sq",
          (af::double2(w_t::*)(cr_mix_t const&) const)
          &w_t::min_max_d_star_sq, (
            arg("miller_indices")))
        .def("stol_sq",
          (double(w_t::*)(mix_t  const&) const)
          &w_t::stol_sq, (
            arg("miller_index")))
        .def("stol_sq",
          (sh_dbl_t(w_t::*)(cr_mix_t const&) const)
          &w_t::stol_sq, (
            arg("miller_indices")))
        .def("two_stol",
          (double(w_t::*)(mix_t const&) const)
          &w_t::two_stol, (
            arg("miller_index")))
        .def("two_stol",
          (sh_dbl_t(w_t::*)(cr_mix_t const&) const)
          &w_t::two_stol, (
            arg("miller_indices")))
        .def("stol",
          (double(w_t::*)(mix_t const&) const)
          &w_t::stol, (
            arg("miller_index")))
        .def("stol",
          (sh_dbl_t(w_t::*)(cr_mix_t const&) const)
          &w_t::stol, (
            arg("miller_indices")))
        .def("d",
          (double(w_t::*)(mix_t const&) const)
          &w_t::d, (
            arg("miller_index")))
        .def("d_frac",
          (double(w_t::*)(frac_mix_t const&) const)
          &w_t::d_frac, (
            arg("miller_index")))
        .def("d",
          (sh_dbl_t(w_t::*)(cr_mix_t const&) const)
          &w_t::d, (
            arg("miller_indices")))
        .def("two_theta",
          (double(w_t::*)(mix_t const&, double, bool) const)
            &w_t::two_theta, (
              arg("miller_index"), arg("wavelength"), arg("deg")=false))
        .def("two_theta",
          (sh_dbl_t(w_t::*)(cr_mix_t const&, double, bool) const)
            &w_t::two_theta, (
              arg("miller_indices"), arg("wavelength"), arg("deg")=false))
        .def("bases_mean_square_difference",
          &w_t::bases_mean_square_difference,
            (arg("other")))
        .def("compare_orthorhombic", &w_t::compare_orthorhombic,
          (arg("other")))
        .def("compare_monoclinic", &w_t::compare_monoclinic,
          (arg("other"), arg("unique_axis"), arg("angular_tolerance")))
        .def("change_of_basis_op_for_best_monoclinic_beta",
          &w_t::change_of_basis_op_for_best_monoclinic_beta, rir())
        .def_pickle(unit_cell_wrappers())
      ;
    }
  };

  inline
  scitbx::vec3<int>
  fractional_unit_shifts_d(fractional<> const& distance_frac)
  {
    return distance_frac.unit_shifts();
  }

  inline
  scitbx::vec3<int>
  fractional_unit_shifts_s_s(
    fractional<> const& site_frac_1,
    fractional<> const& site_frac_2)
  {
    return fractional<>(site_frac_1-site_frac_2).unit_shifts();
  }

  struct distance_mod_1_wrappers
  {
    typedef distance_mod_1 w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("distance_mod_1", no_init)
        .def(init<
          unit_cell const&,
          fractional<> const&,
          fractional<> const&>((
            arg("unit_cell"),
            arg("site_frac_1"),
            arg("site_frac_2"))))
        .add_property("diff_raw", make_getter(&w_t::diff_raw, rbv()))
        .add_property("diff_mod", make_getter(&w_t::diff_mod, rbv()))
        .def_readonly("dist_sq", &w_t::dist_sq)
        .def("unit_shifts", &w_t::unit_shifts)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    //! Forward conversions
    def("d_star_sq_as_stol_sq", (
      double(*)(double)) d_star_sq_as_stol_sq, (arg("d_star_sq")));
    def("d_star_sq_as_stol_sq", (
      af::shared<double>(*)(af::const_ref<double> const &))
      d_star_sq_as_stol_sq, (arg("d_star_sq")));
    def("d_star_sq_as_two_stol", (
      double(*)(double)) d_star_sq_as_two_stol, (arg("d_star_sq")));
    def("d_star_sq_as_two_stol", (
      af::shared<double>(*)(af::const_ref<double> const &))
      d_star_sq_as_two_stol, (arg("d_star_sq")));
    def("d_star_sq_as_stol", (
      double(*)(double)) d_star_sq_as_stol, (arg("d_star_sq")));
    def("d_star_sq_as_stol", (
      af::shared<double>(*)(af::const_ref<double> const &))
      d_star_sq_as_stol, (arg("d_star_sq")));
    def("d_star_sq_as_d", (
      double(*)(double)) d_star_sq_as_d, (arg("d_star_sq")));
    def("d_star_sq_as_d", (
      af::shared<double>(*)(af::const_ref<double> const &))
      d_star_sq_as_d, (arg("d_star_sq")));
    def("d_star_sq_as_two_theta", (
      double(*)(double, double, bool)) d_star_sq_as_two_theta,
      (arg("d_star_sq"), arg("wavelength"), arg("deg")=false));
    def("d_star_sq_as_two_theta", (
      af::shared<double>(*)(af::const_ref<double> const &, double, bool))
      d_star_sq_as_two_theta,
      (arg("d_star_sq"), arg("wavelength"), arg("deg")=false));
    //! Reverse conversions
    def("stol_sq_as_d_star_sq", (
      double(*)(double)) stol_sq_as_d_star_sq, (arg("stol_sq")));
    def("stol_sq_as_d_star_sq", (
      af::shared<double>(*)(af::const_ref<double> const &))
      stol_sq_as_d_star_sq, (arg("stol_sq")));
    def("two_stol_as_d_star_sq", (
      double(*)(double)) two_stol_as_d_star_sq, (arg("two_stol")));
    def("two_stol_as_d_star_sq", (
      af::shared<double>(*)(af::const_ref<double> const &))
      two_stol_as_d_star_sq, (arg("two_stol")));
    def("stol_as_d_star_sq", (
      double(*)(double)) stol_as_d_star_sq, (arg("stol")));
    def("stol_as_d_star_sq", (
      af::shared<double>(*)(af::const_ref<double> const &))
      stol_as_d_star_sq, (arg("stol")));
    def("d_as_d_star_sq", (
      double(*)(double)) d_as_d_star_sq, (arg("d")));
    def("d_as_d_star_sq", (
      af::shared<double>(*)(af::const_ref<double> const &))
      d_as_d_star_sq, (arg("d")));
    def("two_theta_as_d_star_sq", (
      double(*)(double, double, bool)) two_theta_as_d_star_sq,
      (arg("two_theta"), arg("wavelength"), arg("deg")=false));
    def("two_theta_as_d_star_sq", (
      af::shared<double>(*)(af::const_ref<double> const &, double, bool))
      two_theta_as_d_star_sq,
      (arg("two_theta"), arg("wavelength"), arg("deg")=false));
    def("two_theta_as_d", (
      double(*)(double, double, bool)) two_theta_as_d,
      (arg("two_theta"), arg("wavelength"), arg("deg")=false));
    def("two_theta_as_d", (
      af::shared<double>(*)(af::const_ref<double> const &, double, bool))
      two_theta_as_d,
      (arg("two_theta"), arg("wavelength"), arg("deg")=false));

    unit_cell_wrappers::wrap();
    wrap_fast_minimum_reduction();

    def("fractional_unit_shifts", fractional_unit_shifts_d, (
      arg("distance_frac")));
    def("fractional_unit_shifts", fractional_unit_shifts_s_s, (
      arg("site_frac_1"), arg("site_frac_2")));

    distance_mod_1_wrappers::wrap();
  }

} // namespace <anonymous>
}}} // namespace cctbx::uctbx::boost_python

BOOST_PYTHON_MODULE(cctbx_uctbx_ext)
{
  cctbx::uctbx::boost_python::init_module();
}
