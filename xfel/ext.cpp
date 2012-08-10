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
    all_tile_obs_spo.push_back(observed_spot - vec2(tilecenters[itile][0],tilecenters[itile][1]));
    overall_cv += correction_vector;
    sum_sq_cv += correction_vector.length_sq();

  }
};

static boost::python::tuple
bar(){
    // This is going to require some revision to assure it works in an arbitrary
    // principle value region for phi_start_rad and phi_end_rad

    scitbx::af::shared<scitbx::vec3<double> > return_indices;
    scitbx::af::shared<double> return_angles_rad;
    return make_tuple(return_indices,return_angles_rad);
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

    def("foo2", &foo2);

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
    ;

}
}}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xfel_ext)
{
  xfel::boost_python::init_module();

}
