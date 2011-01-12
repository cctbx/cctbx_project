#ifndef RSTBX_SPOT_POSITION_PARTIAL_H
#define RSTBX_SPOT_POSITION_PARTIAL_H
#include <rstbx/diffraction/ewald_sphere.h>

namespace af = scitbx::af;
namespace rstbx {

//! Class to assist calculation of fractional Miller index positions
/*! This class has all the methods of class rotation_angles, but
    also calculates the partial derivatives of the spot position (that is,
    the rotation angle of the spot) with respect to
    variation in the Miller index (i.e., fractional H,K,L values).
    This is convenient for calculating finite-width reciprocal space layers
    within labelit.precession photo.
*/
class partial_spot_position_partial_H: public rotation_angles {
 protected: //member data
  point diangle_dH[2];

 public: //member functions

  //! Constructor using state variables.
  /*! see documentation for the base class (rotation angles) for details.
   */
  partial_spot_position_partial_H(const double& R,
                   const matrix& m,
                   const double& w,
                   const point& axial_direction);
  bool
  operator()(scitbx::vec3<double>const& H);
  point dangle_(const int& idx)const{return diangle_dH[idx];}
};

} //namespace rstbx
#endif //RSTBX_SPOT_POSITION_PARTIAL_H
