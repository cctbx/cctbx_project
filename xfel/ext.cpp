#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/enum.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/constants.h>
#include <scitbx/math/mean_and_variance.h>

using namespace boost::python;

namespace xfel {

struct correction_vector_store {
  typedef scitbx::vec2<double> vec2;
  typedef scitbx::vec3<double> vec3;
  scitbx::af::flex_int tiles;
  scitbx::af::shared<int> tilecounts;
  scitbx::af::shared<vec3> tilecenters;

  scitbx::af::shared<double> radii;
  scitbx::af::shared<vec2> mean_cv;
  scitbx::af::shared<vec2> master_coords;
  scitbx::af::shared<vec2> master_cv;
  scitbx::af::shared<int> master_tiles;
  scitbx::af::shared<vec2> all_tile_obs_spo;
  vec2 overall_cv;
  double sum_sq_cv;

  void initialize_per_tile_sums(){
    tilecounts = scitbx::af::shared<int>(tiles.size()/4,0);
    radii = scitbx::af::shared<double>(tiles.size()/4,0.);
    mean_cv = scitbx::af::shared<vec2>(tiles.size()/4,vec2(0.,0.));
    master_tiles = scitbx::af::shared<int>();
    master_coords = scitbx::af::shared<vec2>();
    master_cv = scitbx::af::shared<vec2>();
    all_tile_obs_spo = scitbx::af::shared<vec2>();
    overall_cv = vec2();
    sum_sq_cv = 0.;
    tilecenters = scitbx::af::shared<vec3>();
    for (int x = 0; x < tiles.size()/4; ++x){
      tilecenters.push_back( vec3(
        (tiles[4*x+0] + tiles[4*x+2])/2.,
        (tiles[4*x+1] + tiles[4*x+3])/2.,
         0.) );
    }
  }

  void
  register_line(double const&a,double const&b,double const&c,double const&d,
                double const&e,double const&f,double const&g,double const&h){
    vec2 observed_center(a,b);
    vec2 refined_center(c,d);
    vec2 observed_spot(e,f);
    vec2 predicted_spot(g,h);
    vec2 prediction = predicted_spot - refined_center;

    vec2 correction_vector = predicted_spot - observed_spot;

    int itile = 0;
    for (int x = 0; x < tiles.size()/4; ++x){
      if (tiles[4*x+0]<predicted_spot[0] && predicted_spot[0]<tiles[4*x+2] &&
          tiles[4*x+1]<predicted_spot[1] && predicted_spot[1]<tiles[4*x+3]){
         itile = x;
         break;
      }
    }
    SCITBX_ASSERT(correction_vector.length() <= 10);
    tilecounts[itile]+=1;
    radii[itile]+=prediction.length();
    mean_cv[itile] = mean_cv[itile] + correction_vector;
    master_tiles.push_back(itile);
    master_cv.push_back(correction_vector);
    master_coords.push_back(prediction);
    all_tile_obs_spo.push_back(observed_spot -
                               vec2(tilecenters[itile][0],tilecenters[itile][1]));
    overall_cv += correction_vector;
    sum_sq_cv += correction_vector.length_sq();
  }

  double
  weighted_average_angle_deg_from_tile(int const& itile) const {

    scitbx::af::shared<vec2> selected_cv;
    scitbx::af::shared<vec2> selected_tile_obs_spo;
    scitbx::af::shared<vec2> translated_correction_vectors;
    scitbx::af::shared<vec2> all_tile_pred_spo;

    for (int x = 0; x < master_tiles.size(); ++x){
      if (master_tiles[x]==itile){
        selected_cv.push_back( master_cv[x] );
        selected_tile_obs_spo.push_back( all_tile_obs_spo[x] );
        translated_correction_vectors.push_back( master_cv[x] - mean_cv[itile] );
        all_tile_pred_spo.push_back( all_tile_obs_spo[x] + master_cv[x] );
      }
    }

    double numerator = 0.;
    double denominator = 0.;

    for (int x = 0; x < selected_cv.size(); ++x){
      vec2 co = selected_tile_obs_spo[x];
      vec2 cp = all_tile_pred_spo[x];
      double co_cp_norm = co.length()*cp.length();
      double co_dot_cp = co*cp;
      double co_cross_cp_coeff = (co[0]*cp[1]-cp[0]*co[1]);
      double sin_theta = co_cross_cp_coeff / co_cp_norm;
      double cos_theta = co_dot_cp / co_cp_norm;
      double angle_deg = std::atan2(sin_theta,cos_theta)/scitbx::constants::pi_180;
      double weight = std::sqrt(co_cp_norm);
      numerator += weight*angle_deg;
      denominator += weight;
    }
    return numerator/denominator;
  }
};

static boost::python::tuple
get_radial_tangential_vectors(correction_vector_store const& L, int const& itile){

    scitbx::vec2<double> radial(0,0);
    scitbx::vec2<double> tangential;
    for (int x = 0; x < L.master_tiles.size(); ++x){
      if (L.master_tiles[x]==itile){
        radial += L.master_coords[x];
      }
    }
    radial = radial.normalize();
    tangential = scitbx::vec2<double>( -radial[1], radial[0] );

    // Now consider 2D Gaussian distribution of all the observations
    scitbx::af::shared<double> radi_projection;
    scitbx::af::shared<double> tang_projection;
    for (int x = 0; x < L.master_tiles.size(); ++x){
      if (L.master_tiles[x]==itile){
        scitbx::vec2<double> recentered_cv = L.master_cv[x] - L.mean_cv[itile];
        radi_projection.push_back( recentered_cv*radial );
        tang_projection.push_back( recentered_cv*tangential );
      }
    }
    scitbx::math::mean_and_variance<double> radistats(radi_projection.const_ref());
    scitbx::math::mean_and_variance<double> tangstats(tang_projection.const_ref());

    return make_tuple(radial,tangential,radistats.mean(),tangstats.mean(),
                      radistats.unweighted_sample_standard_deviation(),
                      tangstats.unweighted_sample_standard_deviation());
}

static boost::python::tuple
get_correction_vector_xy(correction_vector_store const& L, int const& itile){

    scitbx::af::shared<double> xcv;
    scitbx::af::shared<double> ycv;
    for (int x = 0; x < L.master_tiles.size(); ++x){
      if (L.master_tiles[x]==itile){
        xcv.push_back( L.master_cv[x][0] );
        ycv.push_back( L.master_cv[x][1] );
      }
    }
    return make_tuple(xcv,ycv);
}

namespace boost_python { namespace {

  boost::python::tuple
  foo2()
  {
    return boost::python::make_tuple(1,2,3,4);
  }

  void
  init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    def("get_correction_vector_xy", &get_correction_vector_xy);
    def("get_radial_tangential_vectors", &get_radial_tangential_vectors);

    class_<correction_vector_store>("correction_vector_store",init<>())
      .add_property("tiles",
        make_getter(&correction_vector_store::tiles, rbv()),
        make_setter(&correction_vector_store::tiles, dcp()))
      .def("register_line",&correction_vector_store::register_line)
      .def("initialize_per_tile_sums",&correction_vector_store::initialize_per_tile_sums)
      .add_property("tilecounts",
        make_getter(&correction_vector_store::tilecounts, rbv()))
      .add_property("radii",
        make_getter(&correction_vector_store::radii, rbv()),
        make_setter(&correction_vector_store::radii, dcp()))
      .add_property("mean_cv",
        make_getter(&correction_vector_store::mean_cv, rbv()),
        make_setter(&correction_vector_store::mean_cv, dcp()))
      .add_property("master_tiles",
        make_getter(&correction_vector_store::master_tiles, rbv()))
      .add_property("master_cv",
        make_getter(&correction_vector_store::master_cv, rbv()))
      .add_property("overall_cv",
        make_getter(&correction_vector_store::overall_cv, rbv()),
        make_setter(&correction_vector_store::overall_cv, dcp()))
      .add_property("sum_sq_cv",
        make_getter(&correction_vector_store::sum_sq_cv, rbv()))
      .add_property("master_coords",
        make_getter(&correction_vector_store::master_coords, rbv()))
      .add_property("all_tile_obs_spo",
        make_getter(&correction_vector_store::all_tile_obs_spo, rbv()))
      .def("weighted_average_angle_deg_from_tile",
           &correction_vector_store::weighted_average_angle_deg_from_tile)
    ;

}
}}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xfel_ext)
{
  xfel::boost_python::init_module();

}
