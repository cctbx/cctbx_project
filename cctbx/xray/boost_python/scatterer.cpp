#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scatterer_utils.h>
#include <cctbx/crystal/direct_space_asu.h>
#include <cctbx/eltbx/henke.h>
#include <cctbx/eltbx/sasaki.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scatterer_wrappers
  {
    typedef scatterer<> w_t;
    typedef w_t::float_type flt_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("scatterer", no_init)
        .def(init<w_t const&>(arg("other")))
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  flt_t const&,
                  flt_t const&,
                  std::string const&,
                  flt_t const&,
                  flt_t const&>((
          arg("label"),
          arg("site"),
          arg("u_iso"),
          arg("occupancy"),
          arg("scattering_type"),
          arg("fp"),
          arg("fdp"))))
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  scitbx::sym_mat3<flt_t> const&,
                  flt_t const&,
                  std::string const&,
                  flt_t const&,
                  flt_t const&>((
          arg("label"),
          arg("site"),
          arg("u_star"),
          arg("occupancy"),
          arg("scattering_type"),
          arg("fp"),
          arg("fdp"))))
        .add_property("label", make_getter(&w_t::label, rbv()),
                               make_setter(&w_t::label, dcp()))
        .add_property("scattering_type",
          make_getter(&w_t::scattering_type, rbv()),
          make_setter(&w_t::scattering_type, dcp()))
        .add_property("fp", make_getter(&w_t::fp, rbv()),
                            make_setter(&w_t::fp, dcp()))
        .add_property("fdp", make_getter(&w_t::fdp, rbv()),
                             make_setter(&w_t::fdp, dcp()))
        .add_property("site", make_getter(&w_t::site, rbv()),
                              make_setter(&w_t::site, dcp()))
        .add_property("occupancy", make_getter(&w_t::occupancy, rbv()),
                                   make_setter(&w_t::occupancy, dcp()))
        .add_property("u_iso", make_getter(&w_t::u_iso, rbv()),
                               make_setter(&w_t::u_iso, dcp()))
        .add_property("u_star", make_getter(&w_t::u_star, rbv()),
                                make_setter(&w_t::u_star, dcp()))
        .def_readwrite("flags", &w_t::flags)
        .def("set_use_u", &w_t::set_use_u, (arg("iso"),arg("aniso")))
        .def("set_use_u_iso_only", &w_t::set_use_u_iso_only)
        .def("set_use_u_aniso_only", &w_t::set_use_u_aniso_only)
        .def("convert_to_isotropic", &w_t::convert_to_isotropic, (
          arg("unit_cell")))
        .def("convert_to_anisotropic", &w_t::convert_to_anisotropic, (
          arg("unit_cell")))
        .def("is_positive_definite_u",
          (bool(w_t::*)(uctbx::unit_cell const&) const)
            &w_t::is_positive_definite_u, (
          arg("unit_cell")))
        .def("is_positive_definite_u",
          (bool(w_t::*)(uctbx::unit_cell const&, double const&) const)
            &w_t::is_positive_definite_u, (
          arg("unit_cell"),
          arg("u_cart_tolerance")))
        .def("u_iso_or_equiv", &w_t::u_iso_or_equiv, (arg("unit_cell")))
        .def("u_cart_plus_u_iso", &w_t::u_cart_plus_u_iso, (arg("unit_cell")))
        .def("tidy_u",
          (void(w_t::*)(
            uctbx::unit_cell const&,
            sgtbx::site_symmetry_ops const&,
            double const&,
            double const&,
            double const&)) &w_t::tidy_u, (
          arg("unit_cell"),
          arg("site_symmetry_ops"),
          arg("u_min"),
          arg("u_max"),
          arg("anisotropy_min")))
        .def("shift_u",
          (void(w_t::*)(
            uctbx::unit_cell const&,
            double const&)) &w_t::shift_u, (
          arg("unit_cell"),
          arg("u_shift")))
        .def("shift_occupancy",
          (void(w_t::*)(
            double const&)) &w_t::shift_occupancy, (
          arg("q_shift")))
        .def("apply_symmetry",
          (sgtbx::site_symmetry(w_t::*)(
            uctbx::unit_cell const&,
            sgtbx::space_group const&,
            double const&,
            double const&,
            bool)) &w_t::apply_symmetry, (
              arg("unit_cell"),
              arg("space_group"),
              arg("min_distance_sym_equiv")=0.5,
              arg("u_star_tolerance")=0,
              arg("assert_min_distance_sym_equiv")=true))
        .def("apply_symmetry",
          (void(w_t::*)(
            sgtbx::site_symmetry_ops const&,
            double const&)) &w_t::apply_symmetry, (
              arg("site_symmetry_ops"),
              arg("u_star_tolerance")=0))
        .def("apply_symmetry_site", &w_t::apply_symmetry_site, (
          arg("site_symmetry_ops")))
        .def("apply_symmetry_u_star", &w_t::apply_symmetry_u_star, (
          arg("site_symmetry_ops"),
          arg("u_star_tolerance")=0))
        .def("multiplicity", &w_t::multiplicity)
        .def("weight_without_occupancy", &w_t::weight_without_occupancy)
        .def("weight", &w_t::weight)
        .def("report_details", &w_t::report_details, (
          arg("unit_cell"),
          arg("prefix")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_scatterer()
  {
    using namespace boost::python;

    scatterer_wrappers::wrap();

    def("is_positive_definite_u",
      (af::shared<bool>(*)(
        af::const_ref<scatterer<> > const&,
        uctbx::unit_cell const&)) is_positive_definite_u, (
          arg("scatterers"),
          arg("unit_cell")));

    def("is_positive_definite_u",
      (af::shared<bool>(*)(
        af::const_ref<scatterer<> > const&,
        uctbx::unit_cell const&,
        double)) is_positive_definite_u, (
          arg("scatterers"),
          arg("unit_cell"),
          arg("u_cart_tolerance")));

    def("tidy_us",
      (void(*)(
        af::ref<scatterer<> > const&,
        uctbx::unit_cell const&,
        sgtbx::site_symmetry_table const&,
        double u_min,
        double u_max,
        double anisotropy_min)) tidy_us, (
          arg("scatterers"),
          arg("unit_cell"),
          arg("site_symmetry_table"),
          arg("u_min"),
          arg("u_max"),
          arg("anisotropy_min")));

    def("u_star_plus_u_iso",
      (void(*)(
        af::ref<scatterer<> > const&,
        uctbx::unit_cell const&)) u_star_plus_u_iso, (
          arg("scatterers"),
          arg("unit_cell")));

    def("shift_us",
      (void(*)(
        af::ref<scatterer<> > const&,
        uctbx::unit_cell const&,
        double u_min)) shift_us, (
          arg("scatterers"),
          arg("unit_cell"),
          arg("u_shift")));

    def("shift_us",
      (void(*)(
        af::ref<scatterer<> > const&,
        uctbx::unit_cell const&,
        double u_min,
        af::const_ref<std::size_t> const&)) shift_us, (
          arg("scatterers"),
          arg("unit_cell"),
          arg("u_shift"),
          arg("selection")));

    def("shift_occupancies",
      (void(*)(
        af::ref<scatterer<> > const&,
        double,
        af::const_ref<std::size_t> const&)) shift_occupancies, (
          arg("scatterers"),
          arg("q_shift"),
          arg("selection")));

    def("shift_occupancies",
      (void(*)(
        af::ref<scatterer<> > const&,
        double)) shift_occupancies, (
          arg("scatterers"),
          arg("q_shift")));

    def("apply_symmetry_sites",
      (void(*)(
        sgtbx::site_symmetry_table const&,
        af::ref<scatterer<> > const&)) apply_symmetry_sites, (
          arg("site_symmetry_table"),
          arg("scatterers")));

    def("apply_symmetry_u_stars",
      (void(*)(
        sgtbx::site_symmetry_table const&,
        af::ref<scatterer<> > const&,
        double)) apply_symmetry_u_stars, (
          arg("site_symmetry_table"),
          arg("scatterers"),
          arg("u_star_tolerance")=0));

    def("add_scatterers_ext",
      (void(*)(
        uctbx::unit_cell const&,
        sgtbx::space_group const&,
        af::ref<scatterer<> > const&,
        sgtbx::site_symmetry_table&,
        sgtbx::site_symmetry_table const&,
        double,
        double,
        bool,
        bool)) add_scatterers_ext, (
          arg("unit_cell"),
          arg("space_group"),
          arg("scatterers"),
          arg("site_symmetry_table"),
          arg("site_symmetry_table_for_new"),
          arg("min_distance_sym_equiv"),
          arg("u_star_tolerance"),
          arg("assert_min_distance_sym_equiv"),
          arg("non_unit_occupancy_implies_min_distance_sym_equiv_zero")));

    def("change_basis",
      (af::shared<scatterer<> >(*)(
        af::const_ref<scatterer<> > const&,
        sgtbx::change_of_basis_op const&)) change_basis, (
      arg("scatterers"), arg("cb_op")));

    def("expand_to_p1",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        sgtbx::space_group const&,
        af::const_ref<scatterer<> > const&,
        sgtbx::site_symmetry_table const&,
        bool)) expand_to_p1, (
          arg("unit_cell"),
          arg("space_group"),
          arg("scatterers"),
          arg("site_symmetry_table"),
          arg("append_number_to_labels")));

    def("n_undefined_multiplicities",
      (std::size_t(*)(
        af::const_ref<scatterer<> > const&)) n_undefined_multiplicities, (
      arg("scatterers")));

    def("asu_mappings_process",
      (void(*)(
        crystal::direct_space_asu::asu_mappings<>&,
        af::const_ref<scatterer<> > const&,
        sgtbx::site_symmetry_table const&)) asu_mappings_process, (
      arg("asu_mappings"), arg("scatterers"), arg("site_symmetry_table")));

    def("rotate",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        scitbx::mat3<double> const&,
        af::const_ref<scatterer<> > const&)) rotate, (
          arg("unit_cell"),
          arg("rotation_matrix"),
          arg("scatterers")));

    typedef return_value_policy<return_by_value> rbv;
    class_<apply_rigid_body_shift<> >("apply_rigid_body_shift")
      .def(init<af::shared<scitbx::vec3<double> > const&,
                af::shared<scitbx::vec3<double> > const&,
                scitbx::mat3<double> const&,
                scitbx::vec3<double> const&,
                af::const_ref<double> const&,
                uctbx::unit_cell const&,
                af::const_ref<std::size_t> const& >((arg("sites_cart"),
                                              arg("sites_frac"),
                                              arg("rot"),
                                              arg("trans"),
                                              arg("atomic_weights"),
                                              arg("unit_cell"),
                                              arg("selection"))))
      .add_property("sites_frac",  make_getter(
                             &apply_rigid_body_shift<>::sites_frac, rbv()))
      .add_property("sites_cart",  make_getter(
                             &apply_rigid_body_shift<>::sites_cart, rbv()))
      .add_property("center_of_mass",  make_getter(
                             &apply_rigid_body_shift<>::center_of_mass, rbv()))
    ;

    def("set_inelastic_form_factors_from_henke",
      (void(*)(af::ref<scatterer<> > const &,
               eltbx::wavelengths::characteristic,
               bool))
        inelastic_form_factors<eltbx::henke::table>::set,
        (arg("scatterers"), arg("photon"), arg("set_use_fp_fdp")=true));
    def("set_inelastic_form_factors_from_sasaki",
      (void(*)(af::ref<scatterer<> > const &,
               eltbx::wavelengths::characteristic,
               bool))
        inelastic_form_factors<eltbx::henke::table>::set,
        (arg("scatterers"), arg("photon"), arg("set_use_fp_fdp")=true));
    def("set_inelastic_form_factors_from_henke",
      (void(*)(af::ref<scatterer<> > const &, float, bool))
        inelastic_form_factors<eltbx::henke::table>::set,
        (arg("scatterers"), arg("wavelength"), arg("set_use_fp_fdp")=true));
    def("set_inelastic_form_factors_from_sasaki",
      (void(*)(af::ref<scatterer<> > const &, float, bool))
        inelastic_form_factors<eltbx::henke::table>::set,
        (arg("scatterers"), arg("wavelength"), arg("set_use_fp_fdp")=true));
  }

}}} // namespace cctbx::xray::boost_python
