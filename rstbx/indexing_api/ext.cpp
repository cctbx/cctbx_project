#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <rstbx/indexing_api/indexing_api.h>


using namespace boost::python;
using namespace rstbx::indexing_api;

namespace indexing_api{

struct rayleigh_cpp { // a fast C++ version of the Rayleigh distribution class
  /*
  =============================================================================
  Class models a 1-d Rayleigh distribution using one parameter, sigma.

              x                x^2
    pdf = --------- exp(- ------------)
           sigma^2         2 sigma^2

                        x^2
    cdf = 1 - exp(- -----------)
                     2 sigma^2

  The derivative of the cdf with respect to sigma is,

      d(cdf)          x^2               x^2             x
    ---------- = - --------- exp( - -----------) = - ------- pdf
     d(sigma)       sigma^3          2 sigma^2        sigma

  Methods:
    set_parameters
    get_parameters
    estimate_parameters_from_cdf
    pdf
    cdf
    d_cdf_d_sigma
    d_cdf_d_sigma_finite
    cdf_gradients
  -----------------------------------------------------------------------------
  */
  rayleigh_cpp(): sigma(1.),interface("C++"){}
  rayleigh_cpp(const double& s): sigma(s){}

  void set_parameters(scitbx::af::shared<double> p) {
    SCITBX_ASSERT(p.size() == 1);
    sigma = p[0];
  }

  scitbx::af::shared<double> get_parameters(){
    return scitbx::af::shared<double>(1,sigma);
  }

  void estimate_parameters_from_cdf(scitbx::af::shared<double> x_data,scitbx::af::shared<double>y_data){
    //Function estimates the parameter values based on the data (cdf)
    // sigma is the mode of the distribution
    // approximate with the median (cdf = 0.5)
    int midpoint = 0;
    for (int i=0; i < x_data.size(); ++i){
      if (y_data[i] > 0.5){
        midpoint = i;
        break;
      }
    }
    if (midpoint == 0){
      midpoint = x_data.size() - 1;
    }
    sigma = x_data[midpoint];
  }

  double pdf(const double& x){
    //Function returns the probability density function at x
      double x_sigma = x/sigma;
      return (x_sigma/sigma)*std::exp(-0.5*x_sigma*x_sigma);
  }

  scitbx::af::shared<double> pdf(scitbx::af::shared<double> x){
    //Function returns the probability density function at x
    scitbx::af::shared<double> f;
    for (int i = 0; i < x.size(); ++i){
      double x_sigma = x[i]/sigma;
      f.push_back( (x_sigma/sigma)*std::exp(-0.5*x_sigma*x_sigma) );
    }
    return f;
  }

  double cdf(const double& x){
    //Function returns the cumulative distribution function at x
      double x_sigma = x/sigma;
      return 1.0 - std::exp(-0.5*x_sigma*x_sigma);
  }

  scitbx::af::shared<double> cdf(scitbx::af::shared<double> x){
    //Function returns the cumulative distribution function at x
    scitbx::af::shared<double> f;
    for (int i = 0; i < x.size(); ++i){
      double x_sigma = x[i]/sigma;
      f.push_back( 1.0 - std::exp(-0.5*x_sigma*x_sigma) );
    }
    return f;
  }

  double d_cdf_d_sigma(const double& x){
    //Function returns the derivative of the cdf at x with respect to the standard deviation
    double p = pdf(x);
    return -(x/sigma)*p ;
  }

  scitbx::af::shared<double> d_cdf_d_sigma(scitbx::af::shared<double> x){
    //Function returns the derivative of the cdf at x with respect to the standard deviation
    scitbx::af::shared<double> p = pdf(x);
    scitbx::af::shared<double> df;
    for (int i = 0; i < x.size(); ++i){
      df.push_back ( -(x[i]/sigma)*p[i] );
    }
    return df;
  }

  scitbx::af::shared<double> cdf_gradients(const double& x){
    //Function returns a flex.double containing all derivatives
    scitbx::af::shared<double> result;
    result.push_back( d_cdf_d_sigma(x) );
    return result;
  }
  scitbx::af::shared<double> gradients(
    scitbx::af::shared<double> x, const int& nparams, scitbx::af::shared<double>difference){
    //Convenience function to return the gradients in the context of fit_distribution.py
    scitbx::af::shared<double> gradients = scitbx::af::shared<double>(nparams);
    for (int i = 0; i < x.size(); ++i){
      scitbx::af::shared<double> g_i = cdf_gradients(x[i]);
      for (int j = 0; j < nparams; ++j){
        gradients[j] = gradients[j] + difference[i]*g_i[j];
      }
    }
    for (int i = 0; i < gradients.size(); ++i){
      gradients[i] = 2.0*gradients[i];
    }
    return gradients;
  }
  double sigma;
  std::string interface;
};

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
   class_<indexing_api::rayleigh_cpp >("rayleigh_cpp",init< >())
     .def("set_parameters",&indexing_api::rayleigh_cpp::set_parameters, (arg("p")))
     .def("get_parameters",&indexing_api::rayleigh_cpp::get_parameters)
     .def("estimate_parameters_from_cdf",&indexing_api::rayleigh_cpp::estimate_parameters_from_cdf,
         (arg("x_data"),arg("y_data")))
     .def("pdf",(double (indexing_api::rayleigh_cpp::*)(const double&))&indexing_api::rayleigh_cpp::pdf, (arg("x")))
     .def("pdf",(scitbx::af::shared<double> (indexing_api::rayleigh_cpp::*)(scitbx::af::shared<double>))&indexing_api::rayleigh_cpp::pdf, (arg("x")))
     .def("cdf",(double (indexing_api::rayleigh_cpp::*)(const double&))(&indexing_api::rayleigh_cpp::cdf), (arg("x")))
     .def("cdf",(scitbx::af::shared<double> (indexing_api::rayleigh_cpp::*)(scitbx::af::shared<double>))(&indexing_api::rayleigh_cpp::cdf), (arg("x")))
     .def("d_cdf_d_sigma",(double (indexing_api::rayleigh_cpp::*)(const double&))&indexing_api::rayleigh_cpp::d_cdf_d_sigma)
     .def("d_cdf_d_sigma",(scitbx::af::shared<double> (indexing_api::rayleigh_cpp::*)(scitbx::af::shared<double>))&indexing_api::rayleigh_cpp::d_cdf_d_sigma)
     .def("cdf_gradients",(scitbx::af::shared<double> (indexing_api::rayleigh_cpp::*)(const double&))&indexing_api::rayleigh_cpp::cdf_gradients, (arg("x")))
     .add_property("interface",make_getter(&indexing_api::rayleigh_cpp::interface, rbv()))
     .def("gradients",&indexing_api::rayleigh_cpp::gradients, (arg("x"),arg("nparams"),arg("difference")))
   ;
}
