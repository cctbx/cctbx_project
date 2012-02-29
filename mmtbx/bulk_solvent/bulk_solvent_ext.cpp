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

    class_<target_gradients_aniso_ml>("target_gradients_aniso_ml",
                             init<af::const_ref<double> const&,
                                  af::const_ref< std::complex<double> > const&,
                                  af::const_ref< std::complex<double> > const&,
                                  sym_mat3<double> const&,
                                  double const&,
                                  double const&,
                                  af::const_ref<cctbx::miller::index<> > const&,
                                  cctbx::uctbx::unit_cell const&,
                                  cctbx::sgtbx::space_group const&,
                                  af::const_ref<bool> const&,
                                  af::const_ref<double> const&,
                                  af::const_ref<double> const&,
                                  double const&>())
      .def("target", &target_gradients_aniso_ml::target)
      .def("grad_b_cart", &target_gradients_aniso_ml::grad_b_cart)
      .def("grad_ksol", &target_gradients_aniso_ml::grad_ksol)
      .def("grad_bsol", &target_gradients_aniso_ml::grad_bsol)
      .def("grad_k", &target_gradients_aniso_ml::grad_k)
    ;
    //
    class_<bulk_solvent_and_aniso_scale_target_and_grads_ls<> >(
      "bulk_solvent_and_aniso_scale_target_and_grads_ls")
      .def(init<
           f_model::core<double, std::complex<double> > const&,
           f_model::core<double, std::complex<double> > const&,
           double const&,
           af::const_ref<double> const&,
           bool const&,
           bool const&,
           bool const& >(
             (arg("fm1"),
              arg("fm2"),
              arg("twin_fraction"),
              arg("fo"),
              arg("compute_k_sol_grad"),
              arg("compute_b_sol_grad"),
              arg("compute_u_star_grad"))))
      .def(init<
           f_model::core<double, std::complex<double> > const&,
           af::const_ref<double> const&,
           bool const&,
           bool const&,
           bool const& >(
             (arg("fm"),
              arg("fo"),
              arg("compute_k_sol_grad"),
              arg("compute_b_sol_grad"),
              arg("compute_u_star_grad"))))
      .def("target", &bulk_solvent_and_aniso_scale_target_and_grads_ls<>::target)
      .def("grad_u_star", &bulk_solvent_and_aniso_scale_target_and_grads_ls<>::grad_u_star)
      .def("grad_k_sols", &bulk_solvent_and_aniso_scale_target_and_grads_ls<>::grad_k_sols)
      .def("grad_b_sol", &bulk_solvent_and_aniso_scale_target_and_grads_ls<>::grad_b_sol)
   ;
   //
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
           af::const_ref<std::complex<double> > const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const&,
           af::const_ref<double, af::mat_grid> const& >(
             (arg("f_model"),
              arg("f_obs"),
              arg("miller_indices"),
              arg("adp_constraint_matrix"))))
      .def(init<
           af::const_ref<std::complex<double> > const&,
           af::const_ref<double> const&,
           af::const_ref<cctbx::miller::index<> > const&,
           cctbx::uctbx::unit_cell const& >(
             (arg("f_model"),
              arg("f_obs"),
              arg("miller_indices"),
              arg("unit_cell"))))
      .add_property("u_star_independent",
        make_getter(&aniso_u_scaler<>::u_star_independent, rbv()))
      .add_property("a", make_getter(&aniso_u_scaler<>::a, rbv()))
   ;

    def("k_mask_and_k_overall_grid_search",
      (af::tiny<double, 2>(*)
        (af::const_ref<double>                const& f_obs,
         af::const_ref<std::complex<double> > const& f_calc,
         af::const_ref<std::complex<double> > const& f_mask,
         af::const_ref<double>                const& k_mask_range
         )) k_mask_and_k_overall_grid_search);
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

   def("set_k_mask_to_cubic_polynom",
      (af::shared<double>(*)
        (af::const_ref<double> const& ss,
         double                const& ss_cutoff,
         af::tiny<double, 4>   const& coeffs)) set_k_mask_to_cubic_polynom);
   ;

   def("set_to_liear_interpolated",
      (af::shared<double>(*)
        (af::const_ref<double> const& ss,
         double                const& k,
         double                const& b,
         af::const_ref<bool>   const& selection,
         af::shared<double>           data)) set_to_liear_interpolated);
   ;
   //
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
