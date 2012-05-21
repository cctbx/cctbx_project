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

rstbx::scattering_list::scattering_list(scitbx::af::shared<cctbx::miller::index<> > reflections,
                           const cctbx::crystal_orientation& Ori,
                           scitbx::vec3<double> beam_vector_B,
                           scitbx::vec2<double> full_pass,
                           const double& resolution,
                           const double& detector_distance){

    scitbx::vec3<double> XTD(0.,0.,detector_distance);//crystal to detector vector
    scitbx::vec3<double> DX = XTD + scitbx::vec3<double>(1.,0.,0.); //detector x
    scitbx::vec3<double> DY = XTD + scitbx::vec3<double>(0.,1.,0.); //detector y
    double numerator = scitbx::mat3<double>(DX[0],DY[0],0.,
                                            DX[1],DY[1],0.,
                                            DX[2],DY[2],0.).determinant() -
                       scitbx::mat3<double>(XTD[0],DY[0],0.,
                                            XTD[1],DY[1],0.,
                                            XTD[2],DY[2],0.).determinant() +
                       scitbx::mat3<double>(XTD[0],DX[0],0.,
                                            XTD[1],DX[1],0.,
                                            XTD[2],DX[2],0.).determinant() -
                       scitbx::mat3<double>(XTD[0],DX[0],DY[0],
                                            XTD[1],DX[1],DY[1],
                                            XTD[2],DX[2],DY[2]).determinant();

    for (int x = 0; x < reflections.size(); ++x){
      cctbx::miller::index<> hkl = reflections[x];
      scitbx::vec3<double> hkld (hkl[0],hkl[1],hkl[2]);
      scitbx::vec3<double> H = Ori.reciprocal_matrix()*hkld;
      if (H.length()==0.0) { continue; }
      if (1./H.length() < resolution) { continue; }//resolution cutoff

      double t1 = 0.5 * (H*H) / (-beam_vector_B*H);
      if (t1 <= 0) { continue; }

      scitbx::vec3<double>C = t1 * -beam_vector_B;//actual vector to center Ewald Sphere

      double Clen = C.length();
      if (1./Clen < full_pass[1] || 1./Clen > full_pass[0]) { continue; }

      scitbx::vec3<double> H1 =  H - C;

      double denominator = scitbx::mat3<double>(DX[0],DY[0],H1[0],
                                                DX[1],DY[1],H1[1],
                                                DX[2],DY[2],H1[2]).determinant() -
                           scitbx::mat3<double>(XTD[0],DY[0],H1[0],
                                                XTD[1],DY[1],H1[1],
                                                XTD[2],DY[2],H1[2]).determinant() +
                           scitbx::mat3<double>(XTD[0],DX[0],H1[0],
                                                XTD[1],DX[1],H1[1],
                                                XTD[2],DX[2],H1[2]).determinant();

      double t = numerator/denominator;
      mm_coord_result.push_back( -t*H1 );
      reflections_result.push_back( hkl );
    }
}

/* reflection prediction implementation - detailed in header file */

rstbx::reflection_prediction::reflection_prediction(
  const scitbx::vec3<double> & _axis,
  const scitbx::vec3<double> & _s0,
  const scitbx::mat3<double> & _ub,
  //const scitbx::vec3<double> & _origin,
  //const scitbx::vec3<double> & _fast,
  //const scitbx::vec3<double> & _slow,
  //const double & f_min,
  //const double & f_max,
  //const double & s_min,
  //const double & s_max
  const sensor_type & _sensor):
    reflection_range( _axis, _s0,  _ub )
{
  axis = _axis;
  s0 = _s0;
  ub = _ub;
  sensor = _sensor;
}

bool rstbx::reflection_prediction::operator()(scitbx::vec3<double> const & hkl,
                                              const double & angle)
{
  scitbx::vec3<double> s, q;

  s = (ub * hkl).rotate_around_origin(axis, angle);
  q = (s + s0).normalize();

  if (!this->reflection_range::operator()(q,s)) {
    //printf("Too close to the spindle for full measurement. %4d %4d %4d angle %7.4f\n",
    //  int(hkl[0]),int(hkl[1]),int(hkl[2]),angle*180./scitbx::constants::pi);
    /* come back to this later.  These spots are not fully recorded since they are too close
       to the spindle, but we may want to remember their location for future creation
       of overlap masks
    */
    return false;
  }

  return intersect(q); // sets prediction if true
}

bool rstbx::reflection_prediction::intersect(scitbx::vec3<double> const & ray)
{
  scitbx::vec3<double> r;
  double ray_dot_n, x, y;

  ray_dot_n = ray * sensor.get_normal();
  if (ray_dot_n == 0) return false;

  r = (ray * sensor.get_distance() / ray_dot_n) - sensor.get_origin();

  x = r * sensor.get_dir1();
  y = r * sensor.get_dir2();

  scitbx::vec2<double> lim1 = sensor.get_lim1();
  scitbx::vec2<double> lim2 = sensor.get_lim2();

  if (x < lim1[0]) return false;
  if (y < lim2[0]) return false;
  if (x > lim1[1]) return false;
  if (y > lim2[1]) return false;

  prediction[0] = x;
  prediction[1] = y;
  return true;
}

scitbx::vec2<double> rstbx::reflection_prediction::get_prediction()
{
  return scitbx::vec2<double>(prediction[0], prediction[1]);
}
