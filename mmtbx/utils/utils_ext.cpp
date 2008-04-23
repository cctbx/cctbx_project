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
    typedef boost::python::arg arg_;
    typedef return_value_policy<return_by_value> rbv;

    class_<fit_hoh<> >("fit_hoh")
       .def(init<
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            cctbx::fractional<> const&,
            double const&,
            cctbx::uctbx::unit_cell const& >((arg_("site_frac_o"),
                                              arg_("site_frac_h1"),
                                              arg_("site_frac_h2"),
                                              arg_("site_frac_peak1"),
                                              arg_("site_frac_peak2"),
                                              arg_("angular_shift"),
                                              arg_("unit_cell"))))
       .add_property("site_cart_o_fitted",
          make_getter(&fit_hoh<>::site_cart_o_fitted, rbv()))
       .add_property("site_cart_h1_fitted",
          make_getter(&fit_hoh<>::site_cart_h1_fitted, rbv()))
       .add_property("site_cart_h2_fitted",
          make_getter(&fit_hoh<>::site_cart_h2_fitted, rbv()))
       .def("dist_best", &fit_hoh<>::dist_best)
    ;

    def("select_water_by_distance",
         (af::shared<std::size_t>(*)
               (af::shared<vec3<double> > const&,
                af::shared<std::string> const&,
                af::shared<std::size_t> const&,
                double const&,
                double const&,
                double const&,
                cctbx::uctbx::unit_cell const&)) select_water_by_distance,
                                                          (arg_("sites_frac_all"),
                                                           arg_("element_symbols_all"),
                                                           arg_("water_selection_o"),
                                                           arg_("dist_max"),
                                                           arg_("dist_min_mac"),
                                                           arg_("dist_min_sol"),
                                                           arg_("unit_cell")))
   ;
       def("correct_drifted_waters",
         (void(*)
               (af::ref<vec3<double> > const&,
                af::const_ref<vec3<double> > const&,
                af::const_ref<bool> const&,
                cctbx::uctbx::unit_cell const&)) correct_drifted_waters,
                                                          (arg_("sites_frac_all"),
                                                           arg_("sites_frac_peaks"),
                                                           arg_("water_selection"),
                                                           arg_("unit_cell")))
   ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::utils

BOOST_PYTHON_MODULE(mmtbx_utils_ext)
{
  mmtbx::utils::init_module();
}
