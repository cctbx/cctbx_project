#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <rstbx/indexing_api/indexing_api.h>


using namespace boost::python;
using namespace rstbx::indexing_api;

namespace indexing_api{

struct find_green_bar {
  find_green_bar(const scitbx::af::shared<double> rayleigh_cdf_x,
                 const scitbx::af::shared<double> rayleigh_cdf,
                 const scitbx::af::shared<double> dr,
                 const scitbx::af::shared<double> x,
                 const double& sd ):
                 is_set(false){
    const double* ptr_rayleigh_cdf_x = rayleigh_cdf_x.begin();
    const double* ptr_rayleigh_cdf = rayleigh_cdf.begin();
    const double* ptr_x = x.begin();
    const double* ptr_dr = dr.begin();

    for (std::size_t i=0; i < rayleigh_cdf_x.size(); ++i){
      double mx = ptr_rayleigh_cdf_x[i];
      double my = ptr_rayleigh_cdf[i];
      for (std::size_t j=1; j < dr.size(); ++j){
        double upper_x = ptr_dr[j];
        double upper_y = ptr_x[j];
        double lower_x = ptr_dr[j-1];
        double lower_y = ptr_x[j-1];
        if ((my >= lower_y) && (my < upper_y)){
          if ((sd <= (upper_x - mx)) && ((lower_x - mx) > 0.0)){
            //sd_data = ((mx,my),(lower_x,lower_y))
            sd_mx = mx;
            sd_my = my;
            sd_lower_x = lower_x;
            sd_lower_y = lower_y;
            is_set = true;
            radius_outlier_index = j-1;
            limit_outlier = lower_x;
            break;
          }
        }
        if (is_set){
          break;
        }
      }
    }
  }

  bool is_set;
  double sd_mx,sd_my,limit_outlier,sd_lower_x,sd_lower_y;
  int radius_outlier_index;

};

}

BOOST_PYTHON_MODULE(rstbx_indexing_api_ext)
{

   def("cpp_absence_test",cpp_absence_test);

   class_<dps_extended, bases<rstbx::dps_core> >("dps_extended",init< >())
     .def("getData",&dps_extended::getData)
     .def("setData",&dps_extended::setData)
     .def("refine_direction",&dps_extended::refine_direction,
          (arg("candidate"),arg("current_grid"),
           arg("target_grid")))
   ;

   def("raw_spot_positions_mm_to_reciprocal_space_xyz",
     ( scitbx::af::shared< scitbx::vec3<double> > (*) (
       rstbx::pointlist,dxtbx::model::Detector const&, double const&,
       scitbx::vec3<double> const& , scitbx::vec3<double> const&, scitbx::af::shared<int>) )
     raw_spot_positions_mm_to_reciprocal_space_xyz);
   def("raw_spot_positions_mm_to_reciprocal_space_xyz",
     ( scitbx::af::shared< scitbx::vec3<double> > (*) (
       rstbx::pointlist,dxtbx::model::Detector const&, double const&,
       scitbx::vec3<double> const& , scitbx::af::shared<int>) )
     raw_spot_positions_mm_to_reciprocal_space_xyz);

  typedef return_value_policy<return_by_value> rbv;
  class_<indexing_api::find_green_bar>("find_green_bar",
    init<const scitbx::af::shared<double>, const scitbx::af::shared<double>,
              const scitbx::af::shared<double>, const scitbx::af::shared<double>,
              const double& > ((
              arg("rayleigh_cdf_x"),arg("rayleigh_cdf"),arg("dr"),arg("x"),arg("sd"))))
    .add_property("is_set",make_getter(&indexing_api::find_green_bar::is_set, rbv()))
    .add_property("sd_mx",make_getter(&indexing_api::find_green_bar::sd_mx, rbv()))
    .add_property("sd_my",make_getter(&indexing_api::find_green_bar::sd_my, rbv()))
    .add_property("limit_outlier",make_getter(&indexing_api::find_green_bar::limit_outlier, rbv()))
    .add_property("sd_lower_x",make_getter(&indexing_api::find_green_bar::sd_lower_x, rbv()))
    .add_property("sd_lower_y",make_getter(&indexing_api::find_green_bar::sd_lower_y, rbv()))
    .add_property("radius_outlier_index",make_getter(&indexing_api::find_green_bar::radius_outlier_index, rbv()))
  ;
}
