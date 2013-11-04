#include <scitbx/constants.h>
#include <rstbx/dps_core/direction.h>
#include <rstbx/dps_core/directional_fft.h>

namespace pd = rstbx;
namespace constants = scitbx::constants;

pd::Direction::Direction(const double & psi, const double & phi):
  psi(psi),phi(phi) {
  initialize();
  dvec = point(std::sin(psi)*std::cos(phi),
               std::sin(psi)*std::sin(phi),
               std::cos(psi));
}

pd::Direction::Direction(const point & inputdvec):
  dvec(inputdvec) {
  initialize();
  uc_length = dvec.length();
  dvec = dvec.normalize();
  phi = std::atan2(dvec[1],dvec[0]); //formulae only valid for unit vectors
  psi = std::acos(dvec[2]);          //formulae only valid for unit vectors
}

pd::Direction::Direction():
  dvec(scitbx::vec3<double>(0.0,0.0,0.0)),
  psi(0.0),
  phi(0.0){initialize();}

void
pd::Direction::initialize(){
  kmax = 1000000;
  kval = 1.0e300;
  m = kmax;
  delta_p = kval;
  pmin = kval;
  uc_length = kval;
}

bool
pd::Direction::is_nearly_collinear(const Direction & other) const {
  //only valid if these are unit vectors
  double dot = this->dvec * other.dvec;
  double inner_angle = std::acos(std::min(1.0,(double)std::abs(dot)));
  return inner_angle < 10.0*constants::pi_180 ||
         inner_angle > 170.0*constants::pi_180;  //10 degrees
}

bool pd::kvalcmp::operator() (const pd::Direction& a, const pd::Direction& b){
  if (a.kval==b.kval) {return false;}
  return (a.kval > b.kval);
}

void
pd::Direction::extract_directional_properties(fftptr dfft,const bool PS){
  kmax = dfft->kmax();
  kval = dfft->kval();
  kval0 = dfft->kval0();
  kval2 = dfft->kval2();
  kval3 = dfft->kval3();
  m = dfft->m;
  delta_p = dfft->delta_p;
  pmin = dfft->pmin;
  uc_length = dfft->kmax()/ (m*delta_p);
  if (PS) {
    ff = dfft->power_spectrum();
  }
}
