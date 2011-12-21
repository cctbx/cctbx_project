#include <algorithm>
#include <scitbx/error.h>//SCITBX_CHECK_POINT
#include <spotfinder/core_toolbox/distl.h>
#include <spotfinder/core_toolbox/spotfilter.h>

namespace af = scitbx::af;
namespace di = spotfinder::distltbx;

typedef af::shared<int>::iterator PtrTyp;
typedef bool(*funcptr)(di::w_spot const&,di::SpotFilterAgent*);
typedef scitbx::vec3<double> point;
typedef scitbx::mat3<double> matrix;

double radius(di::w_spot const& spot,di::SpotFilterAgent* const& agent){
  double dx = (spot.max_pxl_x()*(agent->pixel_size)) - agent->xbeam;
  double dy = (spot.max_pxl_y()*(agent->pixel_size)) - agent->ybeam;
  // never permit a radius of zero as it ends up in the denominator later on.
  return std::max(agent->pixel_size,std::sqrt( dx*dx + dy*dy ));
}

bool ice_ring_test(di::w_spot const& spot,di::SpotFilterAgent* agent){
  double r = radius(spot,agent);
  for(
   af::shared<Distl::icering>::const_iterator ice_ring =
     agent->icerings.begin();
   ice_ring != agent->icerings.end();
   ++ice_ring) {
     if (
      r > agent->pixel_size*(std::sqrt(ice_ring->lowerr2)) - 0.4 &&
      r < agent->pixel_size*(std::sqrt(ice_ring->upperr2)) + 0.4
     ) {return false;}
  }
  // add a 0.4mm annulus above & below ring to allow for inaccurate center
  return true;
}

bool resolution_test(di::w_spot const& spot,di::SpotFilterAgent* agent){
  //omit spots too close to the beamstop
  double r = radius(spot,agent);
  return r > agent->arguments[0];
}

bool resolution_test_nztt(di::w_spot const& spot,di::SpotFilterAgent* agent){
  //omit spots too close to the beamstop
  return spot.resolution < agent->arguments[0];
}

bool lo_pass_resolution_test(di::w_spot const& spot,di::SpotFilterAgent* agent){
  //omit high-resolution spots
  double r = radius(spot,agent);
  return r < agent->arguments[0];
}

bool lo_pass_resolution_test_nztt(di::w_spot const& spot,di::SpotFilterAgent* agent){
  //omit high-resolution spots
  double r = radius(spot,agent);
  return spot.resolution > agent->arguments[0];
}

bool modal_test(di::w_spot const& spot,di::SpotFilterAgent* agent){
  return spot.nmaxima<=agent->arguments[0] &&
         spot.bodypixels.size()>agent->arguments[1];
}

bool low_skew_test(di::w_spot const& spot,di::SpotFilterAgent* agent){
  return spot.skewness() <= agent->arguments[0];
}

bool intensity_test(di::w_spot const& spot,di::SpotFilterAgent* agent){
  return (spot.intensity() > agent->arguments[0]) &&
         (spot.intensity() < agent->arguments[1]) &&
         (spot.intensity() < agent->arguments[2]);
}

funcptr
function_selector(std::string const& ap){
  if (ap=="ice_ring_test") return ice_ring_test;
  if (ap=="resolution_test") return resolution_test;
  if (ap=="resolution_test_nztt") return resolution_test_nztt;
  if (ap=="lo_pass_resolution_test") return lo_pass_resolution_test;
  if (ap=="lo_pass_resolution_test_nztt") return lo_pass_resolution_test_nztt;
  if (ap=="modal_test") return modal_test;
  if (ap=="low_skew") return low_skew_test;
  if (ap=="intensity_inlier") return intensity_test;
SCITBX_CHECK_POINT;
  return NULL; //non-void function must have return value
}

void
di::SpotFilterAgent::precompute_resolution(spot_list_t masterlist,
                             double const& twotheta,
                             af::shared<double> rotation_axis,
                             af::shared<double> camera_convention){
    SCITBX_ASSERT(rotation_axis[0]==0.0);
    SCITBX_ASSERT(rotation_axis[1]==1.0);
    SCITBX_ASSERT(rotation_axis[2]==0.0);
    double x = rotation_axis[0];
    double y = rotation_axis[1];
    double z = rotation_axis[2];
    point camera_beam = point(xbeam,ybeam,distance);//probably only valid for E=(0,1,0)
    matrix W = std::cos(twotheta) * matrix(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0) +
        (1.- std::cos(twotheta)) * matrix(x*x,x*y,x*z,y*x,y*y,y*z,z*x,z*y,z*z) +
        std::sin(twotheta) * matrix(0,-z,y,z,0,-x,-y,x,0);
    point camera_normal = point(0,0,1.0);
    matrix camera_to_film=matrix(camera_convention[0],camera_convention[1],0., // map 2-D matrix given by
                   camera_convention[2],camera_convention[3],0.,        // user into 3-D matrix
                   0.,0.,1.);
    matrix film_2_camera = camera_to_film.inverse();

    //Update the w_spot array with calculated resolutions.
    for (spot_list_t::iterator spotptr = masterlist.begin();
       spotptr != masterlist.end(); ++spotptr) {

       point observation_frame_X_f = point(spotptr->ctr_mass_x()*pixel_size,
                                      spotptr->ctr_mass_y()*pixel_size,
                                      0.0);
       //insert here, transformed spot =
       //   SXYC.select(rawspot,labelit_commands.spot_convention)
          //Commentary:  if this is too time consuming, formula can be
          // abbreviated by rewriting the matrix multiplications in 2D,
          // omitting the direction parallel to the spindle.  But this
          // would limit the generality of the formulae.
       observation_frame_X_f = camera_beam + film_2_camera * observation_frame_X_f;
       point laboratory_frame_X_f = W * observation_frame_X_f;
       //point X_c = ( distance / (camera_normal * laboratory_frame_X_f) ) *
       //       laboratory_frame_X_f;
       //double unrotated_radius = std::sqrt(X_c[0] * X_c[0] + X_c[1] * X_c[1]);
       //double theta = std::atan2(unrotated_radius,distance)/2.0;
       double spot_two_theta = std::acos((camera_normal * laboratory_frame_X_f)/ laboratory_frame_X_f.length() );
       spotptr->resolution = wavelength/(2.0 * std::sin(spot_two_theta/2.));
    }
}

af::shared<int>
di::SpotFilterAgent::filter(spot_list_t masterlist,
                            af::shared<int> parentselection,
                            std::string test){
  af::shared<int>z; z.reserve(parentselection.size());

  funcptr function = function_selector(test);

  for (int x=0; x<parentselection.size(); ++x){
    if (function(masterlist[parentselection[x]],this)) {
      z.push_back(parentselection[x]);
    }
  }
  return z;
}
struct cmp
{
    di::spot_list_t mstr;
    cmp(di::spot_list_t mstr):mstr(mstr){
    }
    bool operator()(const int& a1, const int& a2) const
    {
        return mstr[a1].resolution > mstr[a2].resolution;
    }
};

af::shared<int>
di::SpotFilterAgent::resolution_sort(spot_list_t masterlist,
                                     af::shared<int> parentselection){
  //First update the w_spot array with calculated resolutions.
  for (spot_list_t::iterator spotptr = masterlist.begin();
       spotptr != masterlist.end(); ++spotptr) {
       double spotradius = radius(*spotptr,this);
       double theta = std::atan2(spotradius,distance)/2.0;
       spotptr->resolution = wavelength/(2.0 * std::sin(theta));
  }
  cmp mycmp(masterlist);
  //sort according to resolution
  af::shared<int> to_be_sorted = parentselection.deep_copy();
  std::sort<PtrTyp,cmp>(to_be_sorted.begin(),to_be_sorted.end(),mycmp);
  return to_be_sorted;
}

af::shared<int>
di::SpotFilterAgent::resolution_sort_nztt(spot_list_t masterlist,
                                     af::shared<int> parentselection){
  cmp mycmp(masterlist);
  //sort according to resolution
  af::shared<int> to_be_sorted = parentselection.deep_copy();
  std::sort<PtrTyp,cmp>(to_be_sorted.begin(),to_be_sorted.end(),mycmp);
  return to_be_sorted;
}

af::shared<double>
di::SpotFilterAgent::get_resolution(spot_list_t masterlist,
                                     af::shared<int> parentselection){
  af::shared<double> z;
  for (af::shared<int>::iterator iptr = parentselection.begin();
       iptr != parentselection.end(); ++iptr) {
       z.push_back(masterlist[*iptr].resolution);
  }
  return z;
}

af::shared<double>
di::SpotFilterAgent::get_property(spot_list_t masterlist,
                              af::shared<int> parentselection,
                              std::string property){
  //This duplicates the functionality of the get_resolution method (i.e,
      // there should be a general function that does both resolution & xy
  //first implementation, property is always the x & y coordinates of spot
  af::shared<double> results;
  if (property=="maximum_pixel") {
    for (int i=0; i<parentselection.size(); ++i){
      results.push_back(masterlist[parentselection[i]].ctr_mass_x());
      results.push_back(masterlist[parentselection[i]].ctr_mass_y());
    }
  } else if (property=="center_of_mass") {
    for (int i=0; i<parentselection.size(); ++i){
      results.push_back(masterlist[parentselection[i]].max_pxl_x());
      results.push_back(masterlist[parentselection[i]].max_pxl_y());
    }
  } else if (property=="skewness") {
    for (int i=0; i<parentselection.size(); ++i){
      results.push_back(masterlist[parentselection[i]].skewness());
    }
  }
  return results;
}

void
di::SpotFilterAgent::set_arguments(af::ref<double>const& args){
  arguments = af::shared<double>();
  for (af::ref<double>::const_iterator i = args.begin();
       i != args.end(); ++i) {
     arguments.push_back(*i);
  }
}

struct cmp_to_be_generalized
{
    di::spot_list_t mstr;
    cmp_to_be_generalized(di::spot_list_t mstr):mstr(mstr){
    }
    bool operator()(const int& a1, const int& a2) const
    {
        return mstr[a1].peakintensity > mstr[a2].peakintensity;
    }
};

af::shared<int>
di::SpotFilterAgent::order_by(spot_list_t masterlist,
                              af::shared<int> parentselection,
                              std::string criterion){
  //first implementation, criterion is always intensity
  cmp_to_be_generalized mycmp(masterlist);
  //sort according to resolution
  af::shared<int> to_be_sorted = parentselection.deep_copy();
  std::sort<PtrTyp,cmp_to_be_generalized>(to_be_sorted.begin(),to_be_sorted.end(),mycmp);
  return to_be_sorted;
}
