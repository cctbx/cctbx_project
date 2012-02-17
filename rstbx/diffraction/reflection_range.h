#ifndef RSTBX_REFLECTION_RANGE_H
#define RSTBX_REFLECTION_RANGE_H
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/crystal_orientation.h>

namespace rstbx {

typedef double Angle;
typedef double AngleRange;

//! Class determines the Bragg spot reflecting range under rotation photography
/*! Expressions are given by:
    Greenhough TJ & Helliwell JR (1982), Oscillation
      Camera Data Processing: Reflecting Range and Prediction of Partiality.
      1. Conventional X-ray Sources. J. Appl. Cryst. 15, 338-351.
    Earlier references are given in that paper.
*/

class reflection_range {

 scitbx::vec3<double> unit_axis;
 Angle mosaicity_rad;
 cctbx::uctbx::unit_cell uc;
 double wavelength;
 scitbx::vec3<double> s0;

 //variables change with every diffracted ray (every call to operator())
 double d_star_sq_unitless, epsilon, Xi_sq;

 public:
  bool use_gh1982a;

  inline
  reflection_range( const scitbx::vec3<double>& rotation_axis,
                    const scitbx::vec3<double>& sample_to_source,
                    const scitbx::mat3<double>& UB
                  ): unit_axis(rotation_axis.normalize()),
                     uc(cctbx::crystal_orientation(UB,true).unit_cell()),
                     wavelength(1./std::sqrt(sample_to_source*sample_to_source)),
                     s0(sample_to_source),use_gh1982a(false)
                     {
  }

  inline void
  set_mosaicity( const double& mos ){
    mosaicity_rad = mos;
  }

  inline void
  set_rocking_curve( std::string const& rc ){
    if (rc=="gh1982a") {use_gh1982a=true;}
  }

  //! Determine if a Miller index is 100% stimulated by rotation
  /*!
      Evaluate GH eq.(I.6).  If 4 * Xi**2 > B**2 then the reciprocal-lattice volume of
      the Bragg spot fully passes through the Ewald sphere.
   */
  inline bool operator()(scitbx::vec3<double> const& q,// normalized diffracted ray vector
                         scitbx::vec3<double> const& s){ //essentially the d-star vector
    if (!use_gh1982a) {return true;}
    d_star_sq_unitless = wavelength * wavelength * s * s;
    double d_star = std::sqrt(d_star_sq_unitless);
    // zeta is the projection of d* onto the rotation axis, per GH Fig. 1
    double zeta = wavelength * s * unit_axis;
    // xi is the component of D* perpendicular to the rotation axis.
    Xi_sq = d_star_sq_unitless - zeta*zeta;

    //Use the approximation for reciprocal lattice volume in GH eqn.(V.6)
    //  the effective mosaic spread, full diameter, is Delta = mosaicity_rad
    double cos_two_theta = wavelength * (q * s0);
    double cos_theta = std::sqrt(0.5 * (1. + cos_two_theta));
    epsilon = 0.5 * mosaicity_rad * d_star * cos_theta;
    double B = d_star_sq_unitless - epsilon*epsilon + 2.*epsilon;
    return 4. * Xi_sq >= B*B;
  }

  inline double lorentz_factor(){
    return 2./std::sqrt( 4. * Xi_sq - d_star_sq_unitless * d_star_sq_unitless );
  }

  inline double get_full_width(){
    return 2. * epsilon * lorentz_factor();
  }
};
} //namespace rstbx
#endif //RSTBX_REFLECTION_RANGE_H
