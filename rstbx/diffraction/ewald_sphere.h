#ifndef RSTBX_EWALD_SPHERE_H
#define RSTBX_EWALD_SPHERE_H
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/crystal_orientation.h>

namespace af = scitbx::af;
namespace rstbx {

typedef double Angle;
typedef double AngleRange;

struct ewald_sphere_base_model{
  typedef scitbx::vec3<double>                 point;
  typedef scitbx::mat3<double>                 matrix;
  double R;           //R=the high resolution limit in Angstrom
  matrix orientation; //the A reciprocal matrix as defined in Rossmann & LABELIT papers
  double wavelength,srsq;//srsq turns out to be 1/wavelength squared
  point e_axial_direction;
  point spherecenter; //Center of Ewald sphere in Rossmann xyz frame
  point H;            //the miller indices
  double inv_maxdsq;
    //rejection criteria from Rossmann(1979)

//member functions:
  ewald_sphere_base_model(const double& R,
                       const matrix& m,
                       const double& w,
                       const point& axial_direction);
  void setH(const point& inH);
  void setH(const cctbx::miller::index<>& inH);
};

//! Class for determining if a Miller index diffracts
/*! Class rotation_angles determines if a given Miller index
    meets the geometric reflection conditions, given an
    orientation, rotation axis vector, limiting resolution,
    and wavelength.
*/
class rotation_angles: public ewald_sphere_base_model {

 public: //member functions

  //! Constructor using state variables.
  /*!
      R = limiting resolution.  Units can be any unit of length, but must be
          consistent throughout.  Typically given in Angstroms or nanometers.
      <P>
      m = orientation matrix.   This reciprocal space matrix is a generalization
          of the fractionalization matrix:
      <pre>
                   / A*x B*x C*x \
      Matrix A* =  | A*y B*y C*y |
                   \ A*z B*z C*z /
      </pre>
      A*, B*, and C* are the reciprocal unit vectors, and components along the
      x,y, and z laboratory axes are given in units of inverse length.
      <P>
      w = wavelength of incident X-rays, given in units of length
      <P>
      axial direction = vector defining the right-handed rotation axis of
      the crystal in laboratory space.
   */
  rotation_angles(const double& R,
                   const matrix& m,
                   const double& w,
                   const point& axial_direction);

  //! Constructor using a base class instance.
  /*!
      The base class, ewald_sphere_base_model, can be constructed using
      same four state variables listed in the above constructor.
   */
  rotation_angles(const ewald_sphere_base_model&); //copy from base class
  //! Determine if a Miller index reflects.
  /*!
      Function returns a bool indicating whether the Miller index can
      be brought into reflecting conditions, under any rotation about
      the axial direction.  If true, the angles of intersection are
      cached and can be recovered by the get_intersection_angles()
      function.  If false, the behavior of a subsequent call to
      the get_intersection_angles() function is undefined.
   */
  bool operator()(scitbx::vec3<double> const&); //evaluate a Miller index

 protected: //member data
  double a_dot_s; //cache value of unit rotation axis (dot) ewald sphere center
  point intersections[2]; //cache points where rotated H intersects Ewald sphere
  Angle iangle[2]; //phi angular rotation values where H intersects Ewald sphere
  Angle angle_(const int& idx)const{return iangle[idx];}

 public: //access member data
  inline point axis() const {return e_axial_direction;}
  inline double offsetdot() const {return a_dot_s;}

  //! Get the rotation angles where a Miller index reflects.
  /*!
      Rotation angles where a Miller index reflects are given, in
      radians.  The return type is a vector containing two values,
      corresponding to the entry and exit of the spot relative to the
      Ewald sphere.  If both values are equal, then the Miller index
      touches the surface of the Ewald sphere.

      <P> A call to this function must be preceeded by a call to
      operator(HKL), where HKL is the Miller index of interest,
      AND the bool returned must be true.

      <P> For example,

      scitbx::vec3<double> miller(3,4,5);
      if (rotation_angles_instance(miller)){
          scitbx::vec2<Angle> refl = rotation_angles_instance.get_intersection_angles();
          printf("Intersection angles (radians) are %10.7f and %10.7f\n",refl[0],refl[1]);
      }
  */
  inline scitbx::vec2<Angle>
  get_intersection_angles()const{//scitbx::vec3<double> const& miller){
    //this->operator()(miller);
    return scitbx::vec2<Angle>(angle_(0), angle_(1));
  }
};

struct scattering_list {
  scitbx::af::shared<scitbx::vec3<double> >mm_coord_result;
  scitbx::af::shared<cctbx::miller::index<> >reflections_result;

  scattering_list(scitbx::af::shared<cctbx::miller::index<> > refl,
                           const cctbx::crystal_orientation& Ori,
                           scitbx::vec3<double> beam_vector_B,
                           scitbx::vec2<double> full_pass,
                           const double& resolution,
                           const double& detector_distance);

  inline
  scitbx::af::shared<scitbx::vec3<double> >
  mm_coord()const {return mm_coord_result;}

  inline
  scitbx::af::shared<cctbx::miller::index<> >
  reflections()const {return reflections_result;}

};


} //namespace rstbx
#endif //RSTBX_EWALD_SPHERE_H
