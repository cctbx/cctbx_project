#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/bulk_solvent/bulk_solvent.h>
#include <mmtbx/f_model/f_model.h>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace mmtbx { namespace bulk_solvent {
namespace {

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;

    class_<ls_kbp_sol_u_star<> >(
      "ls_kbp_sol_u_star")
      .def(init<
           f_model::core<double, std::complex<double> > const&,
           af::const_ref<double> const&,
           double,
           bool const&,
           bool const&,
           bool const&,
           bool const&,
           bool const& >(
             (arg("f_model"),
              arg("f_obs"),
              arg("scale"),
              arg("kb_sol_grad"),
              arg("p_sol_grad"),
              arg("u_star_grad"),
              arg("kb_sol_curv"),
              arg("p_sol_curv"))))
      .def("target",      &ls_kbp_sol_u_star<>::target)
      .def("grad_u_star", &ls_kbp_sol_u_star<>::grad_u_star)
      .def("grad_k_sols", &ls_kbp_sol_u_star<>::grad_k_sols)
      .def("grad_b_sols", &ls_kbp_sol_u_star<>::grad_b_sols)
      .def("grad_p_sols", &ls_kbp_sol_u_star<>::grad_p_sols)
      .def("curv_k_sols", &ls_kbp_sol_u_star<>::curv_k_sols)
      .def("curv_b_sols", &ls_kbp_sol_u_star<>::curv_b_sols)
      .def("curv_p_sols", &ls_kbp_sol_u_star<>::curv_p_sols)
   ;

   class_<ls_u_star<> >(
      "ls_u_star")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const&,
           af::const_ref<double> const& >(
             (arg("f_model_abs_no_k_total"),
              arg("f_obs"),
              arg("miller_indices"),
              arg("k_anisotropic"))))
      .def("target",      &ls_u_star<>::target)
      .def("grad_u_star", &ls_u_star<>::grad_u_star)
   ;

   class_<f_kb_scaled<> >(
      "f_kb_scaled")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<double> const& >(
             (arg("f1"),
              arg("f2"),
              arg("b_range"),
              arg("ss"))))
      .def("k", &f_kb_scaled<>::k)
      .def("b", &f_kb_scaled<>::b)
      .def("scaled", &f_kb_scaled<>::scaled)
   ;

   class_<add_complex_f_kb_scaled<> >(
      "add_complex_f_kb_scaled")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<double> const& >(
             (arg("f0"),
              arg("f1"),
              arg("f2"),
              arg("k_range"),
              arg("b_range"),
              arg("ss"))))
      .def("k",      &add_complex_f_kb_scaled<>::k)
      .def("b",      &add_complex_f_kb_scaled<>::b)
      .def("r",      &add_complex_f_kb_scaled<>::r)
      .def("scaled", &add_complex_f_kb_scaled<>::scaled)
   ;

   class_<complex_f_kb_scaled<> >(
      "complex_f_kb_scaled")
      .def(init<
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<double> const&,
           af::const_ref<double> const& >(
             (arg("f1"),
              arg("f2"),
              arg("b_range"),
              arg("ss"))))
      .def("k", &complex_f_kb_scaled<>::k)
      .def("b", &complex_f_kb_scaled<>::b)
      .def("scaled", &complex_f_kb_scaled<>::scaled)
   ;

   class_<k_sol_b_sol_k_anisotropic_scaler_twin<> >(
      "k_sol_b_sol_k_anisotropic_scaler_twin")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<double> const&,
           double const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const&,
           cctbx::uctbx::unit_cell const&,
           double const&>(
             (arg("f_obs"),
              arg("f_calc_1"),
              arg("f_calc_2"),
              arg("f_mask_1"),
              arg("f_mask_2"),
              arg("ss"),
              arg("twin_fraction"),
              arg("k_sol_range"),
              arg("b_sol_range"),
              arg("miller_indices"),
              arg("unit_cell"),
              arg("r_ref"))))
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           double const&>(
             (arg("f_obs"),
              arg("f_calc"),
              arg("f_mask"),
              arg("k_total"),
              arg("ss"),
              arg("k_sol_range"),
              arg("b_sol_range"),
              arg("r_ref"))))
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const&,
           double const&>(
             (arg("f_obs"),
              arg("f_calc"),
              arg("f_mask"),
              arg("ss"),
              arg("k_sol_range"),
              arg("b_sol_range"),
              arg("miller_indices"),
              arg("r_ref"))))
      .def("r", &k_sol_b_sol_k_anisotropic_scaler_twin<>::r)
      .def("k_sol", &k_sol_b_sol_k_anisotropic_scaler_twin<>::k_sol)
      .def("b_sol", &k_sol_b_sol_k_anisotropic_scaler_twin<>::b_sol)
      .def("k_mask", &k_sol_b_sol_k_anisotropic_scaler_twin<>::k_mask)
      .def("k_anisotropic", &k_sol_b_sol_k_anisotropic_scaler_twin<>::k_anisotropic)
      .def("updated", &k_sol_b_sol_k_anisotropic_scaler_twin<>::updated)
      .def("u_star", &k_sol_b_sol_k_anisotropic_scaler_twin<>::u_star)
   ;

   class_<overall_and_bulk_solvent_scale_coefficients_analytical<> >(
      "overall_and_bulk_solvent_scale_coefficients_analytical")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<bool> const& >(
             (arg("f_obs"),
              arg("f_calc"),
              arg("f_mask"),
              arg("selection"))))
      .add_property("x", make_getter(
        &overall_and_bulk_solvent_scale_coefficients_analytical<>::x, rbv()))
      .add_property("r", make_getter(
        &overall_and_bulk_solvent_scale_coefficients_analytical<>::r, rbv()))
      .add_property("x_best", make_getter(
        &overall_and_bulk_solvent_scale_coefficients_analytical<>::x_best, rbv()))
      .add_property("r_best", make_getter(
        &overall_and_bulk_solvent_scale_coefficients_analytical<>::r_best, rbv()))
   ;

   class_<bulk_solvent_scale_coefficients_analytical<> >(
      "bulk_solvent_scale_coefficients_analytical")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<std::complex<double> > const&,
           af::const_ref<bool> const& >(
             (arg("f_obs"),
              arg("f_calc"),
              arg("f_mask"),
              arg("selection"))))
      .add_property("x", make_getter(
        &bulk_solvent_scale_coefficients_analytical<>::x, rbv()))
      .add_property("r", make_getter(
        &bulk_solvent_scale_coefficients_analytical<>::r, rbv()))
      .add_property("x_best", make_getter(
        &bulk_solvent_scale_coefficients_analytical<>::x_best, rbv()))
      .add_property("r_best", make_getter(
        &bulk_solvent_scale_coefficients_analytical<>::r_best, rbv()))
   ;

   class_<aniso_u_scaler<> >("aniso_u_scaler")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const& >(
             (arg("f_model_abs"),
              arg("f_obs"),
              arg("miller_indices"))))
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const&,
           af::const_ref<double, af::mat_grid> const& >(
             (arg("f_model_abs"),
              arg("f_obs"),
              arg("miller_indices"),
              arg("adp_constraint_matrix"))))
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const&,
           cctbx::uctbx::unit_cell const& >(
             (arg("f_model_abs"),
              arg("f_obs"),
              arg("miller_indices"),
              arg("unit_cell"))))
      .add_property("u_star_independent",
        make_getter(&aniso_u_scaler<>::u_star_independent, rbv()))
      .add_property("a", make_getter(&aniso_u_scaler<>::a, rbv()))
      .add_property("u_star", make_getter(&aniso_u_scaler<>::u_star, rbv()))
   ;

    def("k_mask_and_k_overall_grid_search",
      (af::tiny<double, 2>(*)
        (af::const_ref<double>                const& f_obs,
         af::const_ref<std::complex<double> > const& f_calc,
         af::const_ref<std::complex<double> > const& f_mask,
         af::const_ref<double>                const& k_mask_range,
         af::const_ref<bool>                  const& selection
         )) k_mask_and_k_overall_grid_search);
   ;

   def("fit_k_exp_b_to_k_total",
      (af::tiny<double, 3>(*)
        (af::const_ref<double> const& data,
         af::const_ref<double> const& ss,
         double                       k_start,
         double                       b_start
         )) fit_k_exp_b_to_k_total);
   ;

   def("ksol_bsol_grid_search",
      (af::shared<double>(*)
        (af::const_ref<double>   const&,
         af::const_ref<std::complex<double> > const&,
         af::const_ref<std::complex<double> > const&,
         af::const_ref<double>   const&,
         af::const_ref<double>   const&,
         af::const_ref<double>   const&,
         double                  const&,
         af::const_ref<double>   const&,
         af::const_ref<double>   const&,
         double                  const&)) ksol_bsol_grid_search);
   ;

   def("complex_f_minus_f_kb_scaled",
      (af::shared<std::complex<double> >(*)
        (af::const_ref<std::complex<double> > const&,
         af::const_ref<std::complex<double> > const&,
         af::const_ref<double>   const&,
         af::const_ref<double>   const&)) complex_f_minus_f_kb_scaled);
   ;

   def("set_k_mask_to_cubic_polynom",
      (af::shared<double>(*)
        (af::const_ref<double> const& ss,
         double                const& ss_cutoff,
         af::tiny<double, 4>   const& coeffs)) set_k_mask_to_cubic_polynom);
   ;

   def("set_to_linear_interpolated",
      (af::shared<double>(*)
        (af::const_ref<double> const& ss,
         double                const& k,
         double                const& b,
         af::const_ref<bool>   const& selection,
         af::shared<double>           data)) set_to_linear_interpolated);
   ;

    def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&)) r_factor);
   ;
   def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&,
         af::const_ref<bool> const&)) r_factor);
   ;
   def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&,
         af::const_ref< std::complex<double> > const&,
         double const&)) r_factor);
   ;
   def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&,
         double const&)) r_factor);
   ;
   def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&,
         af::const_ref<bool> const&,
         double const&)) r_factor);
   ;
   def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref<double> const&,
         double const&)) r_factor);
   ;
   def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref<double> const&)) r_factor);
   ;
   def("r_factor",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&,
         af::const_ref< std::complex<double> > const&,
         double const&,
         double const&)) r_factor);
   ;
   def("scale",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&)) scale);
   ;
   def("scale",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&,
         af::const_ref<bool> const&)) scale);
   ;
   def("scale",
      (double(*)
        (af::const_ref<double> const&,
         af::const_ref< std::complex<double> > const&,
         af::const_ref< std::complex<double> > const&,
         double const&)) scale);
   ;
    def("fb_cart",fb_cart)
   ;
   def("symmetrize_mask",
      (void(*)
        (af::ref<int, af::c_grid<3> > const&,
         af::const_ref<long, af::c_grid<3> > const&)) symmetrize_mask);
  }

} // namespace <anonymous>
}} // namespace mmtbx::bulk_solvent

BOOST_PYTHON_MODULE(mmtbx_bulk_solvent_ext)
{
  mmtbx::bulk_solvent::init_module();
}
