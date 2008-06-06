#include <scitbx/constants.h>
#include <rstbx/dps_core/direction.h>

namespace pd = labelit::dptbx;
namespace constants = scitbx::constants;

pd::Direction::Direction(const double & psi, const double & phi):
  psi(psi),phi(phi) {
  dvec = point(std::sin(psi)*std::cos(phi),
               std::sin(psi)*std::sin(phi),
               std::cos(psi));
}

pd::Direction::Direction(const point & dvec):
  dvec(dvec) {
  phi = std::atan2(dvec[1],dvec[0]); //formulae only valid for unit vectors
  psi = std::acos(dvec[2]);          //formulae only valid for unit vectors
}

pd::Direction::Direction(){}

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
