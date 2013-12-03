#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>
#include <algorithm>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include <rstbx/indexing_api/indexing_api.h>
#include <rstbx/dps_core/direction.h>
#include <cctbx/math/mod.h>

namespace af = scitbx::af;
namespace constants = scitbx::constants;
using rstbx::Direction;
using rstbx::Directional_FFT;
using rstbx::kvalcmp;

af::shared<int>
rstbx::indexing_api::cpp_absence_test(af::shared<cctbx::miller::index<> > hkl,
                     const int& mod, cctbx::miller::index<> vecrep){
  af::shared<int> cum;
  for (std::size_t i=0; i<mod; ++i) {cum.push_back(0);}

  typedef af::shared<cctbx::miller::index<> >::const_iterator For;
  For b = hkl.begin();
  For e = hkl.end();
  for (;b!=e;++b){
    int pattern_sum = (*b)[0]*vecrep[0] + (*b)[1]*vecrep[1] + (*b)[2]*vecrep[2];
    cum[cctbx::math::mod_positive(pattern_sum, mod)]+=1;
  }

  return cum;
}

rstbx::indexing_api::dps_extended::dps_extended():rstbx::dps_core(){}

void
rstbx::indexing_api::dps_extended::setData(const pointlist& raw){
  rawdata=pointlistmm();
  rawdata.reserve(raw.size());
  for (std::size_t i = 0; i< raw.size(); ++i) {
      rawdata.push_back(raw[i]);
  }
}

rstbx::Direction
rstbx::indexing_api::dps_extended::refine_direction(const rstbx::Direction& candidate,
                               const double& current_grid,
                               const double& target_grid) const {
  double incr = current_grid/4.0; //recursively smaller grid with neighbor overlap
  if (incr < target_grid) { return candidate; }
  rstbx::Direction direction0(candidate.dvec);
  rstbx::Direction direction1(candidate.psi + incr, candidate.phi);
  scitbx::vec3<double> rotaxisi = direction0.dvec.cross(direction1.dvec).normalize();
  scitbx::vec3<double> rotaxisj = rotaxisi.cross(direction1.dvec).normalize();

  rstbx::Direction best_refined = candidate;
  for (int i = -2; i < 3; ++i) {
    scitbx::vec3<double> rotate1 = direction0.dvec.unit_rotate_around_origin(rotaxisi,-i*incr);
    for (int j = -2; j < 3; ++j) {
      scitbx::vec3<double> rotate2 = rotate1.unit_rotate_around_origin(rotaxisj,-j*incr);
      SCITBX_ASSERT(std::abs(1.0 - rotate2.length())<0.00001); //rotate2 is unit vector
      rstbx::Direction new_candidate(rotate2);
      rstbx::fftptr dfft( fft_factory(new_candidate) );

      if ( dfft->kval() > best_refined.kval ) {
        new_candidate.extract_directional_properties(dfft);
        best_refined=new_candidate;
      }
    }
  }
  return refine_direction(best_refined,incr,target_grid);
}

double
rstbx::indexing_api::dps_extended::high()const{
  //Go through all film coordinates and get high resolution limit
  // This is NOT a safe function if rawdata.size()==0 due to
  // dereferencing s.end().  This bug causes a painful infinite loop
  // in some runtime environments.  See the protecting test in classify_spots().
  // The real fix will be to split this
  // class into components, with the high() function being a member
  // of the class containing rawdata.

  SCITBX_ASSERT (xyzdata.size() > 0);
  af::shared<double> s(xyzdata.size());
  for (std::size_t i = 0; i<xyzdata.size(); ++i) {
    s[i] =         1./xyzdata[i].length();
  }
  double dmax = *(std::max_element<af::shared<double>::const_iterator>(s.begin(),s.end()));
  //dmax is maximum squared distance (mm^2) from the beam center
  return dmax;
}

af::shared< scitbx::vec3<double> >
rstbx::indexing_api::raw_spot_positions_mm_to_reciprocal_space_xyz(
  pointlist raw_spot_input,
  dxtbx::model::Detector const& detector, double const& inverse_wave,
  scitbx::vec3<double> const& S0_vector, scitbx::vec3<double> const& axis, af::shared<int> panelID
){

  // it is not clear how we will plug in to this function if and when
  // we need to use a derived detector class with a different implemented
  // behavior of get_lab_coord().

  af::shared< scitbx::vec3<double> > reciprocal_space_vectors;
  // with a single panel only, all the convenience functions work for us
  for (int n = 0; n != raw_spot_input.size(); ++n) {
    int pid = panelID[n];

    // tile surface to laboratory transformation
    scitbx::vec2<double> raw_spot(raw_spot_input[n][0],raw_spot_input[n][1]);
    scitbx::vec3<double> lab_direct = detector[pid].get_lab_coord(raw_spot);

    // laboratory direct to reciprocal space xyz transformation
    scitbx::vec3<double> lab_recip = (lab_direct.normalize() * inverse_wave) - S0_vector;

      // raw_spot_input[n][2] MUST be given in degrees.
    reciprocal_space_vectors.push_back( lab_recip.rotate_around_origin(
      axis, -raw_spot_input[n][2] * scitbx::constants::pi_180) );
  }
  return reciprocal_space_vectors;

}

af::shared< scitbx::vec3<double> >
rstbx::indexing_api::raw_spot_positions_mm_to_reciprocal_space_xyz(
  pointlist raw_spot_input,
  dxtbx::model::Detector const& detector, double const& inverse_wave,
  scitbx::vec3<double> const& S0_vector, af::shared<int> panelID
){

  // it is not clear how we will plug in to this function if and when
  // we need to use a derived detector class with a different implemented
  // behavior of get_lab_coord().

  af::shared< scitbx::vec3<double> > reciprocal_space_vectors;
  // with a single panel only, all the convenience functions work for us
  for (int n = 0; n != raw_spot_input.size(); ++n) {
    int pid = panelID[n];

    // tile surface to laboratory transformation
    scitbx::vec2<double> raw_spot(raw_spot_input[n][0],raw_spot_input[n][1]);
    scitbx::vec3<double> lab_direct = detector[pid].get_lab_coord(raw_spot);

    // laboratory direct to reciprocal space xyz transformation
    scitbx::vec3<double> lab_recip = (lab_direct.normalize() * inverse_wave) - S0_vector;

      // raw_spot_input[n][2] MUST be given in degrees.
    reciprocal_space_vectors.push_back(lab_recip);
  }
  return reciprocal_space_vectors;

}
