#include <rstbx/diffraction/partial_spot_position_partial_H.h>
#include <scitbx/array_family/shared.h>

rstbx::partial_spot_position_partial_H::partial_spot_position_partial_H(const double& R,
                   const matrix& m,
                   const double& w,
                   const point& axial_direction):
    rotation_angles(R,m,w,axial_direction){
}

rstbx::ewald_sphere_base_model::point
derivative_of_normalized_vector(rstbx::ewald_sphere_base_model::point const& u,
                                rstbx::ewald_sphere_base_model::point const& du){
  double v_sq = u*u;
  rstbx::ewald_sphere_base_model::point v = std::sqrt(v_sq);
  rstbx::ewald_sphere_base_model::point dv = (1./v) * (u * du);
  return (du * v - u * dv)/v_sq;
}

bool
rstbx::partial_spot_position_partial_H::operator()
  (scitbx::vec3<double>const& H){

  bool result = rstbx::rotation_angles::operator()(H);

  if (result) {
    // Instead of tedious coding to calculate analytical derivatives,
    // just take finite differences, which may be more time
    // consuming but are much easier to implement.  Later if this
    // feature becomes rate limiting, come back and do the calculus.

    //increments of 0.01 Miller index serve as good approximants for
    // finite difference derviatives
    double finite = 0.01;
    af::shared<point> dH;
    dH.push_back(point(finite,0.0,0.0));
    dH.push_back(point(0.0,finite,0.0));
    dH.push_back(point(0.0,0.0,finite));

    rstbx::rotation_angles diff = rstbx::rotation_angles(*this);
    for (int i = 0; i < 3; ++i){
      diff(H+dH[i]);
      scitbx::vec2<Angle> refl = diff.get_intersection_angles();
      for (int idx = 0; idx < 2; ++idx){
        diangle_dH[idx][i] = (refl[idx] - this->angle_(idx))/finite;
      }
    }
  }

  return result;

}
