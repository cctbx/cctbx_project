#ifndef RSTBX_EWALD_SPHERE_H
#define RSTBX_EWALD_SPHERE_H
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/crystal_orientation.h>
#include <rstbx/diffraction/reflection_range.h>
#include <rstbx/bpcx/detector_model/sensor.h>

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

  inline scattering_list(scitbx::af::shared<scitbx::vec3<double> >a,
                         scitbx::af::shared<cctbx::miller::index<> >hkl):
                         mm_coord_result(a),reflections_result(hkl)
  {/* for use as a general spot positions container, no unit of measure or ref frame implied*/}

  inline
  scitbx::af::shared<scitbx::vec3<double> >
  mm_coord()const {return mm_coord_result;}

  inline
  scitbx::af::shared<cctbx::miller::index<> >
  reflections()const {return reflections_result;}

};

/* additional code for reflection prediction - a class which takes the
   rotation axis, beam vector, UB matrix, detector origin, detector fast
   and slow axes, the extent of the detector in these directions w.r.t.
   the origin. assumes rectangular detector, edges colinear with fast and
   slow directions, returns positions in fast, slow direction if observed.

   here is the corresponding Python class:

class reflection_prediction:
    def __init__(self, axis, s0, ub, detector_origin,
                 detector_fast, detector_slow,
                 f_min, f_max, s_min, s_max):
        self._axis = axis
        self._s0 = s0
        self._ub = ub
        self._detector_origin = detector_origin
        self._detector_fast = detector_fast
        self._detector_slow = detector_slow
        self._limits = f_min, f_max, s_min, s_max

        return

    def predict(self, indices, angles):

        detector_normal = self._detector_fast.cross(self._detector_slow)
        distance = self._detector_origin.dot(detector_normal)

        observed_reflection_positions = []

        for hkl, angle in zip(indices, angles):
            s = (self._ub * hkl).rotate_around_origin(self._axis, angle)
            q = (s + self._s0).normalize()

            # check if diffracted ray parallel to detector face

            q_dot_n = q.dot(detector_normal)

            if q_dot_n == 0:
                continue

            r = (q * distance / q_dot_n) - self._detector_origin

            x = r.dot(self._detector_fast)
            y = r.dot(self._detector_slow)

            if x < self._limits[0] or y < self._limits[2]:
                continue
            if x > self._limits[1] or y > self._limits[3]:
                continue

            observed_reflection_positions.append((hkl, x, y, angle))

        return observed_reflection_positions

   initial API: constructor sets everything up, operator() tells you whether
   it was observed (and computes the intersection point if so) and
   get_prediction() returns this as a pair of indices, fast and slow.
   N.B. s0 has length 1.0 / wavelength, and is the vector from the source
   towards the sample.

*/

class reflection_prediction : public reflection_range {

 public:

  typedef rstbx::detector_model::sensor sensor_type;

  reflection_prediction(const scitbx::vec3<double> & axis,
                        const scitbx::vec3<double> & s0,
                        const scitbx::mat3<double> & ub,
                        const sensor_type & sensor);

  bool operator()(scitbx::vec3<double> const & hkl,
                  const double & angle);

  //aim to make intersect virtual and implemented by derived classes
  //see e.g. detector_model in bpcx_regression
  bool intersect(scitbx::vec3<double> const & ray);

  scitbx::vec2<double> get_prediction();
  scitbx::vec3<double> get_s();

 protected:

  scitbx::vec3<double> axis;
  scitbx::vec3<double> s0;
  scitbx::mat3<double> ub;
  //scitbx::vec3<double> origin;
  //scitbx::vec3<double> fast;
  //scitbx::vec3<double> slow;
  //scitbx::vec3<double> normal;
  sensor_type sensor;
  //double limits[4];
  //double distance;

  scitbx::vec2<double> prediction;
  scitbx::vec3<double> s;

};

} //namespace rstbx
#endif //RSTBX_EWALD_SPHERE_H
