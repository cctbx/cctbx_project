#ifndef RSTBX_EWALD_SPHERE_H
#define RSTBX_EWALD_SPHERE_H
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>

namespace af = scitbx::af;
namespace rstbx {

typedef double Angle;
typedef double AngleRange;

struct EwaldSphereBaseModel{
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
  EwaldSphereBaseModel(const double& R,
                       const matrix& m,
                       const double& w,
                       const point& axial_direction);
  void setH(const point& inH);
  void setH(const cctbx::miller::index<>& inH);
};

class RotationAngles: public EwaldSphereBaseModel {

 public: //member functions
  RotationAngles(const double& R,
                   const matrix& m,
                   const double& w,
                   const point& axial_direction);
  RotationAngles(const EwaldSphereBaseModel&); //copy from base class
  bool operator()(scitbx::vec3<double> const&); //evaluate a Miller index

 private: //member data
  double Sr; //radius of Ewald sphere
  double a_dot_s; //cache value of unit rotation axis (dot) ewald sphere center
  point intersections[2]; //cache points where rotated H intersects Ewald sphere
  Angle iangle[2]; //phi angular rotation values where H intersects Ewald sphere
  Angle angle_(const int& idx)const{return iangle[idx];}

 public: //access member data
  inline point axis() const {return e_axial_direction;}
  inline double offsetdot() const {return a_dot_s;}
  inline scitbx::vec2<Angle>
  get_intersection_angles(){//scitbx::vec3<double> const& miller){
    //this->operator()(miller);
    return scitbx::vec2<Angle>(angle_(0), angle_(1));
  }
};


} //namespace rstbx
#endif //RSTBX_EWALD_SPHERE_H
