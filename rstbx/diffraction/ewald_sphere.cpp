#include <rstbx/diffraction/ewald_sphere.h>

rstbx::ewald_sphere_base_model::ewald_sphere_base_model(const double& R,
                   const matrix& m,
                   const double& w,
                   const point& axial_direction):
    R(R),
    orientation(m),
    wavelength(w),
    e_axial_direction(axial_direction){
    spherecenter = point(0.0,0.0,-1./wavelength);
    srsq = spherecenter.length_sq();
    inv_maxdsq = 1./(R*R);
  }

void
rstbx::ewald_sphere_base_model::setH(const point& inH){
    H = inH;
  }

void
rstbx::ewald_sphere_base_model::setH(const cctbx::miller::index<>& inH){
    H = point(inH[0],inH[1],inH[2]);
  }

/*
Use ideas suggested by Mike Hohn 12/8/2005:
given the sphere S and circle Co in plane P, how about this sequence:
    1. find the sphere / plane intersection to get a circle Cn, point, or
       nothing.
    2. for a circle Cn, find intersection of Co and Cn (which are
       coplanar)

  For 1. you can use the distance D from the center of S to the plane P, e.g.,
    http://mathworld.wolfram.com/Point-PlaneDistance.html
  where the plane is easily expressed using
    http://mathworld.wolfram.com/HessianNormalForm.html
  and use D to compute Cn.

  For 2., using a shifted version of this would work:
    http://mathworld.wolfram.com/Circle-CircleIntersection.html

  Not closed form, but easily checked at every step.
*/

rstbx::rotation_angles::rotation_angles(const ewald_sphere_base_model& ew):
    ewald_sphere_base_model(ew){

    //make sure that the axial direction is of unit length
    e_axial_direction = e_axial_direction.normalize();

    //cache the dot product of axial direction and Ewald sphere center
    a_dot_s = e_axial_direction * spherecenter;
}
    // point spherecenter: center of Ewald sphere
    // double Sr: radius of Ewald sphere
    // point e_axial_direction: normal to the plane of rotation

rstbx::rotation_angles::rotation_angles(const double& R,
                   const matrix& m,
                   const double& w,
                   const point& axial_direction):
    ewald_sphere_base_model(R,m,w,axial_direction){
    e_axial_direction = e_axial_direction.normalize();
    a_dot_s = e_axial_direction * spherecenter;
}

bool
rstbx::rotation_angles::operator()(scitbx::vec3<double>const& H){
  setH(H);
  //these are mandatory criteria~should be factored into the base class later
  point Hxyz = orientation * H; //H coordinates in lab orthonormal frame
  double len_sq = (Hxyz).length_sq();

  //condition 1.  Miller index is past the user-defined limiting resolution.
  if (len_sq > inv_maxdsq) {return false;}

  //condition 2.  Miller index is 000.
  if (H[0]==0 && H[1]==0 && H[2]==0) {return false;}

  //get_projection_of_H_onto the rotation_axis (could be positive or negative)
  double DplaneH = Hxyz * e_axial_direction;

  //Express the plane P in which H rotates in Hessian Normal Form
  // P is the locus of points x such that e_axial_direction . x = -DplaneH

  //The distance from a point X0 to the plane P is given by
  // D =  e_axial_direction . X0 + DplaneH
  // In this case we are interested in the distance of the center of the
  // Ewald sphere (this distance can be positive or negative)
  //  Correct sign for a_dot_s provided by Graeme Winter 8/6/2009:
  double D = DplaneH - a_dot_s;

  //Remove this assertion for the Pringle_Shen gonio; spindle not perpendicular to beam
  //SCITBX_ASSERT(DplaneH * D >= 0.); //sphere centers Co and Cn must be on same side of origin

  //condition 3.  Miller index is outside the Ewald Sphere diameter.
  if  ( std::abs(D)>=spherecenter.length() ) {return false;}

  //Consider the intersection of two circles:
  //  Co, the circle of rotation of H.
      point  Cocenter = e_axial_direction * DplaneH;
      double Coradius = (Cocenter - Hxyz).length();

  //  Cn, the circular section through the Ewald sphere.
      point  Cncenter = spherecenter + D * e_axial_direction;
      double Cnradius_squared = (srsq-D*D);

  //  Taking a page from mathworld.wolfram.com, calculate the distance d
  //  between the centers of Co and Cn,
      double d = (Cocenter - Cncenter).length();
      //Remove this assertion for the Pringle_Shen gonio; spindle not perpendicular to beam
      //SCITBX_ASSERT (d-1./wavelength < 1.e-14);

  //  The chord of intersection between Co and Cn lies a
  //  distance x along the (Cocenter - Cncenter) vector
      double fac = d*d - (Coradius*Coradius) + (Cnradius_squared);
      double x = 0.5*(fac)/d;
      double discriminant = 4.*d*d*Cnradius_squared - fac*fac;

      //condition 4.  Miller index is near spindle; never intersects Ewald sphere
      if (discriminant <= 0.0) {return false;}

  //  Calculate half-length of the chord of intersection
      double y = 0.5*std::sqrt(discriminant)/d;
      point chord_direction =
        (e_axial_direction.cross(Cocenter - Cncenter)).normalize();

  //  Two intersection points

      intersections[0]=(Cncenter + x*(Cocenter - Cncenter).normalize() + y*chord_direction)
                       -Cocenter;
      //The intersection point has to lie on the circle Co
      SCITBX_ASSERT (std::abs(intersections[0].length()-Coradius)<1.e-10);
      intersections[1]=(Cncenter + x*(Cocenter - Cncenter).normalize() - y*chord_direction)
                       -Cocenter;

      point Hxyzunit = (Hxyz-Cocenter).normalize();
      //point third_axis = Hxyzunit.cross(e_axial_direction);
      point third_axis = e_axial_direction.cross(Hxyzunit);

      iangle[0]= std::atan2(intersections[0]*third_axis,intersections[0]*Hxyzunit);
      iangle[1]= std::atan2(intersections[1]*third_axis,intersections[1]*Hxyzunit);

  return true;
}
