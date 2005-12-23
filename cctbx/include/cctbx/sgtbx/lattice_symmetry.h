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

  double
  find_max_delta(
    uctbx::unit_cell const& reduced_cell,
    space_group const& group);

  struct reduced_cell_two_fold_info
  {
    sg_mat3 r;
    sg_vec3 u;
    sg_vec3 h;
  };

  // Reformatted output of cctbx/cctbx/examples/reduced_cell_two_folds.py
  static const reduced_cell_two_fold_info reduced_cell_two_folds[] = {
    sg_mat3(-1,-1,-1,0,0,1,0,1,0),sg_vec3(-1,1,1),sg_vec3(0,1,1),
    sg_mat3(-1,-1,0,0,1,0,0,-1,-1),sg_vec3(1,-2,1),sg_vec3(0,1,0),
    sg_mat3(-1,-1,0,0,1,0,0,0,-1),sg_vec3(-1,2,0),sg_vec3(0,1,0),
    sg_mat3(-1,-1,0,0,1,0,0,1,-1),sg_vec3(-1,2,1),sg_vec3(0,1,0),
    sg_mat3(-1,-1,1,0,0,-1,0,-1,0),sg_vec3(1,-1,1),sg_vec3(0,-1,1),
    sg_mat3(-1,0,-1,0,-1,-1,0,0,1),sg_vec3(-1,-1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,-1,0,-1,0,0,0,1),sg_vec3(-1,0,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,-1,0,-1,1,0,0,1),sg_vec3(-1,1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,-1,0,-1,1,-1,0),sg_vec3(0,-1,1),sg_vec3(1,-1,1),
    sg_mat3(-1,0,0,-1,0,1,-1,1,0),sg_vec3(0,1,1),sg_vec3(-1,1,1),
    sg_mat3(-1,0,0,-1,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,-2,1),
    sg_mat3(-1,0,0,-1,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,2,0),
    sg_mat3(-1,0,0,-1,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,2,1),
    sg_mat3(-1,0,0,0,-1,-1,0,0,1),sg_vec3(0,-1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,0,-1,0,-1,-1,1),sg_vec3(0,0,1),sg_vec3(-1,-1,2),
    sg_mat3(-1,0,0,0,-1,0,-1,0,1),sg_vec3(0,0,1),sg_vec3(-1,0,2),
    sg_mat3(-1,0,0,0,-1,0,-1,1,1),sg_vec3(0,0,1),sg_vec3(-1,1,2),
    sg_mat3(-1,0,0,0,-1,0,0,-1,1),sg_vec3(0,0,1),sg_vec3(0,-1,2),
    sg_mat3(-1,0,0,0,-1,0,0,0,1),sg_vec3(0,0,1),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,0,-1,0,0,1,1),sg_vec3(0,0,1),sg_vec3(0,1,2),
    sg_mat3(-1,0,0,0,-1,0,1,-1,1),sg_vec3(0,0,1),sg_vec3(1,-1,2),
    sg_mat3(-1,0,0,0,-1,0,1,0,1),sg_vec3(0,0,1),sg_vec3(1,0,2),
    sg_mat3(-1,0,0,0,-1,0,1,1,1),sg_vec3(0,0,1),sg_vec3(1,1,2),
    sg_mat3(-1,0,0,0,-1,1,0,0,1),sg_vec3(0,1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,0,0,-1,0,-1,0),sg_vec3(0,-1,1),sg_vec3(0,-1,1),
    sg_mat3(-1,0,0,0,0,1,0,1,0),sg_vec3(0,1,1),sg_vec3(0,1,1),
    sg_mat3(-1,0,0,0,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,-2,1),
    sg_mat3(-1,0,0,0,1,0,0,-1,-1),sg_vec3(0,-2,1),sg_vec3(0,1,0),
    sg_mat3(-1,0,0,0,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,1,0),
    sg_mat3(-1,0,0,0,1,0,0,1,-1),sg_vec3(0,2,1),sg_vec3(0,1,0),
    sg_mat3(-1,0,0,0,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,2,1),
    sg_mat3(-1,0,0,1,0,-1,-1,-1,0),sg_vec3(0,-1,1),sg_vec3(-1,-1,1),
    sg_mat3(-1,0,0,1,0,1,1,1,0),sg_vec3(0,1,1),sg_vec3(1,1,1),
    sg_mat3(-1,0,0,1,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,-2,1),
    sg_mat3(-1,0,0,1,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,2,0),
    sg_mat3(-1,0,0,1,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,2,1),
    sg_mat3(-1,0,1,0,-1,-1,0,0,1),sg_vec3(1,-1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,1,0,-1,0,0,0,1),sg_vec3(1,0,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,1,0,-1,1,0,0,1),sg_vec3(1,1,2),sg_vec3(0,0,1),
    sg_mat3(-1,1,-1,0,0,-1,0,-1,0),sg_vec3(-1,-1,1),sg_vec3(0,-1,1),
    sg_mat3(-1,1,0,0,1,0,0,-1,-1),sg_vec3(-1,-2,1),sg_vec3(0,1,0),
    sg_mat3(-1,1,0,0,1,0,0,0,-1),sg_vec3(1,2,0),sg_vec3(0,1,0),
    sg_mat3(-1,1,0,0,1,0,0,1,-1),sg_vec3(1,2,1),sg_vec3(0,1,0),
    sg_mat3(-1,1,1,0,0,1,0,1,0),sg_vec3(1,1,1),sg_vec3(0,1,1),
    sg_mat3(0,-1,-1,-1,0,1,0,0,-1),sg_vec3(-1,1,0),sg_vec3(-1,1,1),
    sg_mat3(0,-1,-1,0,-1,0,-1,1,0),sg_vec3(-1,0,1),sg_vec3(-1,1,1),
    sg_mat3(0,-1,0,-1,0,0,-1,1,-1),sg_vec3(-1,1,1),sg_vec3(-1,1,0),
    sg_mat3(0,-1,0,-1,0,0,0,0,-1),sg_vec3(-1,1,0),sg_vec3(-1,1,0),
    sg_mat3(0,-1,0,-1,0,0,1,-1,-1),sg_vec3(1,-1,1),sg_vec3(-1,1,0),
    sg_mat3(0,-1,1,-1,0,-1,0,0,-1),sg_vec3(-1,1,0),sg_vec3(1,-1,1),
    sg_mat3(0,-1,1,0,-1,0,1,-1,0),sg_vec3(1,0,1),sg_vec3(1,-1,1),
    sg_mat3(0,0,-1,-1,-1,1,-1,0,0),sg_vec3(-1,1,1),sg_vec3(-1,0,1),
    sg_mat3(0,0,-1,0,-1,0,-1,0,0),sg_vec3(-1,0,1),sg_vec3(-1,0,1),
    sg_mat3(0,0,-1,1,-1,-1,-1,0,0),sg_vec3(-1,-1,1),sg_vec3(-1,0,1),
    sg_mat3(0,0,1,-1,-1,-1,1,0,0),sg_vec3(1,-1,1),sg_vec3(1,0,1),
    sg_mat3(0,0,1,0,-1,0,1,0,0),sg_vec3(1,0,1),sg_vec3(1,0,1),
    sg_mat3(0,0,1,1,-1,1,1,0,0),sg_vec3(1,1,1),sg_vec3(1,0,1),
    sg_mat3(0,1,-1,0,-1,0,-1,-1,0),sg_vec3(-1,0,1),sg_vec3(-1,-1,1),
    sg_mat3(0,1,-1,1,0,-1,0,0,-1),sg_vec3(1,1,0),sg_vec3(-1,-1,1),
    sg_mat3(0,1,0,1,0,0,-1,-1,-1),sg_vec3(-1,-1,1),sg_vec3(1,1,0),
    sg_mat3(0,1,0,1,0,0,0,0,-1),sg_vec3(1,1,0),sg_vec3(1,1,0),
    sg_mat3(0,1,0,1,0,0,1,1,-1),sg_vec3(1,1,1),sg_vec3(1,1,0),
    sg_mat3(0,1,1,0,-1,0,1,1,0),sg_vec3(1,0,1),sg_vec3(1,1,1),
    sg_mat3(0,1,1,1,0,1,0,0,-1),sg_vec3(1,1,0),sg_vec3(1,1,1),
    sg_mat3(1,-1,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,1,1),
    sg_mat3(1,-1,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,1,0),
    sg_mat3(1,-1,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,-1,1),
    sg_mat3(1,0,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,0,1),
    sg_mat3(1,0,0,-1,-1,0,-1,0,-1),sg_vec3(-2,1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,-1,-1,0,0,0,-1),sg_vec3(-2,1,0),sg_vec3(1,0,0),
    sg_mat3(1,0,0,-1,-1,0,1,0,-1),sg_vec3(2,-1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,0,-1,0,-1,0,-1),sg_vec3(-2,0,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(1,0,0),
    sg_mat3(1,0,0,0,-1,0,1,0,-1),sg_vec3(2,0,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,1,-1,0,-1,0,-1),sg_vec3(-2,-1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,1,-1,0,0,0,-1),sg_vec3(2,1,0),sg_vec3(1,0,0),
    sg_mat3(1,0,0,1,-1,0,1,0,-1),sg_vec3(2,1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,0,1),
    sg_mat3(1,1,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,-1,1),
    sg_mat3(1,1,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,1,0),
    sg_mat3(1,1,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,1,1)
  };

  /*! References:
        Y. Le Page
        The derivation of the axes of the conventional unit cell from the
        dimensions of the Buerger-reduced cell
        J. Appl. Cryst. (1982). 15, 255-259

        Andrey A. Lebedev, Alexei A. Vagin & Garib N. Murshudov
        Acta Cryst. (2006). D62, 83-95.
        http://journals.iucr.org/d/issues/2006/01/00/ba5089/index.html
        Appendix A1. Algorithms used in the determination of twinning
          operators and their type of merohedry
  */
  space_group
  group(
    uctbx::unit_cell const& reduced_cell,
    double max_delta=3.,
    bool enforce_max_delta_for_generated_two_folds=true);

}}} // namespace cctbx::sgtbx::lattice_symmetry

#endif // CCTBX_SGTBX_LATTICE_SYMMETRY_H
