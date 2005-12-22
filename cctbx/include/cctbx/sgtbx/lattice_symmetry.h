#ifndef CCTBX_SGTBX_LATTICE_SYMMETRY_H
#define CCTBX_SGTBX_LATTICE_SYMMETRY_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/math/utils.h>
#include <scitbx/constants.h>

namespace cctbx { namespace sgtbx { namespace lattice_symmetry {

  inline
  uc_mat3
  two_fold_matrix_from_axis_direction(uc_vec3 const& ev_cart)
  {
    double f = 2. / ev_cart.length_sq();
    double x = ev_cart[0];
    double y = ev_cart[1];
    double z = ev_cart[2];
    return uc_mat3(f*x*x-1., f*x*y,    f*x*z,
                   f*y*x,    f*y*y-1., f*y*z,
                   f*z*x,    f*z*y,    f*z*z-1.);
  }

  inline
  uc_mat3
  n_fold_operator_from_axis_direction(
    uc_vec3 const& cart,
    int n,
    int sense=1)
  {
    if (n == 1) return uc_mat3(1,0,0,0,1,0,0,0,1);
    if (n == 2) return two_fold_matrix_from_axis_direction(cart);
    CCTBX_ASSERT(sense == 1 || sense == -1);
    CCTBX_ASSERT(n == 1 || n == 2 || n == 3 || n == 4 || n == 6);
    uc_vec3 v = cart.normalize();
    double angle = scitbx::constants::two_pi / (sense * n);
    double c = std::cos(angle);
    double s = std::sin(angle);
    double d = 1 - c;
    return uc_mat3(
      c+d*v[0]*v[0],      d*v[0]*v[1]-s*v[2], d*v[0]*v[2]+s*v[1],
      d*v[0]*v[1]+s*v[2], c+d*v[1]*v[1],      d*v[1]*v[2]-s*v[0],
      d*v[0]*v[2]-s*v[1], d*v[1]*v[2]+s*v[0], c+d*v[2]*v[2]);
  }

  inline
  sg_mat3
  as_integer_matrix(uc_mat3 const& m)
  {
    sg_mat3 result;
    for(std::size_t i=0;i<m.size();i++) {
      result[i] = scitbx::math::iround(m[i]);
    }
    return result;
  }

  struct potential_axis_t
  {
    potential_axis_t() {}

    potential_axis_t(
      sg_vec3 const& u_,
      sg_vec3 const& h_,
      int abs_uh_)
    :
      u(u_),
      h(h_),
      abs_uh(abs_uh_)
    {}

    sg_vec3 u;
    sg_vec3 h;
    int abs_uh;
  };

  /*! Reference:
        Y. Le Page
        The derivation of the axes of the conventional unit cell from the
        dimensions of the Buerger-reduced cell
        J. Appl. Cryst. (1982). 15, 255-259
  */
  class group_search
  {
    public:
      group_search(int modulus=2)
      :
        modulus_(modulus)
      {}

      std::size_t
      n_potential_axes();

      space_group
      operator()(
        uctbx::unit_cell const& niggli_cell,
        double max_delta=3.,
        bool enforce_max_delta_for_generated_two_folds=true);

    private:
      int modulus_;
      std::vector<potential_axis_t> potential_axes_;

      void
      compute_potential_axes();
  };

  double
  find_max_delta(
    uctbx::unit_cell const& niggli_cell,
    space_group const& group);

}}} // namespace cctbx::sgtbx::lattice_symmetry

#endif // CCTBX_SGTBX_LATTICE_SYMMETRY_H
