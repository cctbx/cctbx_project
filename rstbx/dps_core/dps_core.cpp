#include <algorithm>

#include <rstbx/dps_core/dps_core.h>

namespace af = scitbx::af;
namespace pd = rstbx;

double dps_ai_round(const double &x) {
  return std::floor(x+0.5);
}

scitbx::vec3<double>
dps_fractional_to_miller(const scitbx::vec3<double>& p) {
  return scitbx::vec3<double>(
                       dps_ai_round(p[0]),
                       dps_ai_round(p[1]),
                       dps_ai_round(p[2]));
}

pd::dps_core::dps_core():
   granularity(5),outliers_marked(false){}
   /*! granularity of the internal projection of observed reciprocal space
      vectors onto the directional axis.  granularity value is a fixed
      parameter of the Rossmann DPS algorithm */

void
pd::dps_core::setMaxcell(const double& a){
  amax=a;
}

pd::fftptr
pd::dps_core::fft_factory(const pd::Direction& angle) const {
  fftptr rv(new Directional_FFT(angle, xyzdata, granularity, amax));
  return rv;
}

void
pd::dps_core::setSolutions(af::shared<pd::Direction> nwsoln){
  //sort
  typedef af::shared<pd::Direction>::iterator Ran;
  Ran  first = nwsoln.begin();
  Ran last = nwsoln.end();
  std::sort<Ran,kvalcmp>(first,last,kvalcmp());
  hemisphere_solutions = nwsoln;

  //sztype is unsigned, so return now if array size is zero, otherwise
  // get a segmentation fault
  if (hemisphere_solutions.size()==0) {return;}

  //pop out solutions that are linearly dependent
  for (sztype i = hemisphere_solutions.size()-1; i>0; --i) {
    for (sztype j = 0; j<i; ++j) {
      if (hemisphere_solutions[i].is_nearly_collinear(hemisphere_solutions[j])){
        hemisphere_solutions[i].kval=0.0; //mark for later removal
        break;
      }
    }
  }
  first = hemisphere_solutions.begin(); last=hemisphere_solutions.end();
  std::sort<Ran,kvalcmp>(first,last,kvalcmp());
  while (hemisphere_solutions[hemisphere_solutions.size()-1].kval==0.0) {
    hemisphere_solutions.pop_back();}
}

void
pd::dps_core::set_presorted_solutions(af::shared<pd::Direction> nwsoln){
  typedef af::shared<pd::Direction>::iterator Ran;

  hemisphere_solutions = nwsoln;
  //sztype is unsigned, so return now if array size is zero, otherwise
  // get a segmentation fault
  if (hemisphere_solutions.size()==0) {return;}

  //pop out solutions that are linearly dependent
  for (sztype i = hemisphere_solutions.size()-1; i>0; --i) {
    for (sztype j = 0; j<i; ++j) {
      if (hemisphere_solutions[i].is_nearly_collinear(hemisphere_solutions[j])){
        hemisphere_solutions[i].kval=-10.0; //mark for later removal
        break;
      }
    }
  }

  Ran first = hemisphere_solutions.begin(); Ran last=hemisphere_solutions.end();
  std::sort<Ran,kvalcmp>(first,last,kvalcmp());
  while (hemisphere_solutions[hemisphere_solutions.size()-1].kval==-10.0) {
    hemisphere_solutions.pop_back();}
}

af::shared<pd::Direction>
pd::dps_core::getSolutions() const {return hemisphere_solutions;}

void pd::dps_core::setOrientation(const Orientation& input_o){
  orientation = input_o;
  this->classify_spots();
}

void pd::dps_core::set_orientation_direct_matrix(
  const scitbx::mat3<double>& dm){
  orientation = Orientation(dm,false);
  this->classify_spots();
}

void pd::dps_core::set_orientation_reciprocal_matrix(
  const scitbx::mat3<double>& rm){
  orientation = Orientation(rm, true);
  this->classify_spots();
}

void
pd::dps_core::classify_spots() {
  reset_spots();
}

void
pd::dps_core::reset_spots() {
  scitbx::mat3<double> Ainv = orientation.direct_matrix();
  // initial construction of status
  if (!outliers_marked) {
    status=dps_statuslist();
    status.reserve(xyzdata.size());
    Estatus=dps_statuslist();
    Estatus.reserve(xyzdata.size());
    obsdata = pointlistmm();
    obsdata.reserve(xyzdata.size());
    hkldata = pointlistmm();
    hkldata.reserve(xyzdata.size());
    for (std::size_t i = 0; i< xyzdata.size(); ++i) {
      status.push_back(GOOD);
      Estatus.push_back(NONE);
      obsdata.push_back(Ainv*xyzdata[i]);
      hkldata.push_back(dps_fractional_to_miller(obsdata[i]));
    }
  }
  // prevent 'OUTLIER' from being changed
  else {
    for (std::size_t i = 0; i< xyzdata.size(); ++i) {
      if (status[i] != OUTLIER) {
        status[i] = GOOD;
      }
      Estatus[i] = NONE;
      obsdata[i] = Ainv*xyzdata[i];
      hkldata[i] = dps_fractional_to_miller(obsdata[i]);
    }
  }
}
double
pd::dps_core::rmsdev() const {
  int N = 0;
  double sumsq = 0;
  pointlistmm obs = observed(); //builtin filter for GOOD spots only
  pointlistmm::const_iterator b = obs.begin();
  pointlistmm::const_iterator e = obs.end();
  for (; b!=e; ++b) {
      N+=1;
      point deviation((*b)[0] - dps_ai_round((*b)[0]),
                      (*b)[1] - dps_ai_round((*b)[1]),
                      (*b)[2] - dps_ai_round((*b)[2]));
      sumsq+=deviation*deviation;
  }
  return std::sqrt(sumsq/N);
}

pd::pointlistmm
pd::dps_core::observed() const {
  pointlistmm observed_f;
  observed_f.reserve(xyzdata.size());
  for (sztype i = 0; i < xyzdata.size(); ++i) {
    if (status[i]==GOOD) { // New design; never get all spots, just GOOD ones
      observed_f.push_back(obsdata[i]);}
  }
  return observed_f;
}

pd::hkllistmm
pd::dps_core::hklobserved() const {
  return hklobserved( observed() );
}

pd::hkllistmm
pd::dps_core::hklobserved(const pointlistmm & spots) const {
  hkllistmm observed_miller;

  observed_miller.reserve(spots.size());

  pointlistmm::const_iterator b = spots.begin();
  pointlistmm::const_iterator e = spots.end();

  for (; b!=e; ++b) {
    observed_miller.push_back(dps_fractional_to_miller(*b));
  }
  return observed_miller;
}
