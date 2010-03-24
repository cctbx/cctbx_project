#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/utils/utils.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace mmtbx { namespace utils {
namespace {

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;

    class_<fit_hoh<> >("fit_hoh")
       .def(init<
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            double const&,
            cctbx::uctbx::unit_cell const& >((arg("site_frac_o"),
                                              arg("site_frac_h1"),
                                              arg("site_frac_h2"),
                                              arg("site_frac_peak1"),
                                              arg("site_frac_peak2"),
                                              arg("angular_shift"),
                                              arg("unit_cell"))))
       .add_property("site_cart_o_fitted",
          make_getter(&fit_hoh<>::site_cart_o_fitted, rbv()))
       .add_property("site_cart_h1_fitted",
          make_getter(&fit_hoh<>::site_cart_h1_fitted, rbv()))
       .add_property("site_cart_h2_fitted",
          make_getter(&fit_hoh<>::site_cart_h2_fitted, rbv()))
       .def("dist_best", &fit_hoh<>::dist_best)
    ;
    //
    class_<density_distribution_per_atom<> >("density_distribution_per_atom")
       .def(init<
            af::ref<vec3<double> > const&,
            af::const_ref<vec3<double> > const&,
            af::const_ref<double> const&,
            cctbx::uctbx::unit_cell const& >((
                                              arg("sites_frac_atoms"),
                                              arg("sites_frac_peaks"),
                                              arg("density_values"),
                                              arg("unit_cell"))))
       .def("distances", &density_distribution_per_atom<>::distances)
       .def("map_values", &density_distribution_per_atom<>::map_values)
    ;
    //

    def("select_water_by_distance",
         (af::shared<std::size_t>(*)
               (af::shared<vec3<double> > const&,
                af::shared<std::string> const&,
                af::shared<std::size_t> const&,
                double const&,
                double const&,
                double const&,
                cctbx::uctbx::unit_cell const&)) select_water_by_distance,
                                                          (arg("sites_frac_all"),
                                                           arg("element_symbols_all"),
                                                           arg("water_selection_o"),
                                                           arg("dist_max"),
                                                           arg("dist_min_mac"),
                                                           arg("dist_min_sol"),
                                                           arg("unit_cell")))
   ;
       def("correct_drifted_waters",
         (void(*)
               (af::ref<vec3<double> > const&,
                af::const_ref<vec3<double> > const&,
                af::const_ref<bool> const&,
                cctbx::uctbx::unit_cell const&)) correct_drifted_waters,
                                                          (arg("sites_frac_all"),
                                                           arg("sites_frac_peaks"),
                                                           arg("water_selection"),
                                                           arg("unit_cell")))
   ;

   def("create_twin_mate",
         (af::shared<cctbx::miller::index<> >(*)
               (af::const_ref<cctbx::miller::index<> > const&,
                scitbx::mat3<double>)) create_twin_mate,
                  (arg("miller_indices"),
                   arg("twin_law_matrix")))
   ;

   def("apply_twin_fraction",
         (af::shared<double>(*)
               (af::const_ref<double> const&,
                af::const_ref<double> const&,
                double const&)) apply_twin_fraction,
                  (arg("amplitude_data_part_one"),
                   arg("amplitude_data_part_two"),
                   arg("twin_fraction")))
   ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::utils

BOOST_PYTHON_MODULE(mmtbx_utils_ext)
{
  mmtbx::utils::init_module();
}
