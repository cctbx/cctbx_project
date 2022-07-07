#include <smtbx/refinement/constraints/shared.h>
#include <scitbx/math/euler_angles.h>
#include <boost/lambda/lambda.hpp>
#include <scitbx/array_family/tiny_algebra.h>

namespace smtbx { namespace refinement { namespace constraints {
  // shared u_star
  void
  shared_u_star
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    u_star_parameter *u = reference();
    value = u->value;

    if (jacobian_transpose == 0) {
      return;
    }
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int i=0; i<6; i++) {
      jt.col(index() + i) = jt.col(u->index() + i);
    }
  }

  // shared rotated u_star
  void
  shared_rotated_u_star
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    // rotation vector
    cart_t rv = direction()->direction(unit_cell);
    const u_star_parameter *ref = reference();
    independent_scalar_parameter *ang = angle();
    const double
      ca = cos(ang->value),
      sa = sin(ang->value),
      t = 1.0-ca;
    const tensor_rank_2_t u_c =
      cctbx::adptbx::u_star_as_u_cart(unit_cell, ref->value);
    // build the rotation matrix and rotate the u_cart
    const scitbx::mat3<double> rm(
      t*rv[0]*rv[0] + ca,       t*rv[0]*rv[1] + sa*rv[2], t*rv[0]*rv[2] - sa*rv[1],
      t*rv[0]*rv[1] - sa*rv[2], t*rv[1]*rv[1] + ca,       t*rv[1]*rv[2] + sa*rv[0],
      t*rv[0]*rv[2] + sa*rv[1], t*rv[2]*rv[1] - sa*rv[0], t*rv[2]*rv[2] + ca
    );
    scitbx::mat3<double> u_t = rm*u_c*rm.transpose();
    // update the value
    value = cctbx::adptbx::u_cart_as_u_star(unit_cell,
      tensor_rank_2_t(u_t[0], u_t[4], u_t[8], u_t[1], u_t[2], u_t[5]));

    if (jacobian_transpose == 0) {
      return;
    }
    // convenience accessor array for the symmetric matrix...
    static const int sym_acs[] = {0,4,8,1,2,5};
    sparse_matrix_type &jt = *jacobian_transpose;
    // transforms for the jacobian values
    scitbx::mat3<double> jtm = unit_cell.fractionalization_matrix()*
      rm*unit_cell.orthogonalization_matrix(), jtm_t = jtm.transpose();
    for (int i=0; i<jt.n_rows(); i++)  {
      tensor_rank_2_t t;
      bool zero = true;
      for (int j=0; j<6; j++) {
        if ((t[j] = jt(i, ref->index() + j)) != 0) {
          zero = false;
        }
      }
      if (zero) {
        for (int j = 0; j < 6; j++) {
          jt(i, index() + j) = 0;
        }
      }
      else {
        scitbx::mat3<double> x = jtm*t*jtm_t;
        for (int j = 0; j < 6; j++) {
          jt(i, index() + j) = x[sym_acs[j]];
        }
      }
    }
    if (ang->is_variable()) {
      // derivative of the rotation matrix by angle
      const scitbx::mat3<double> rmd(
        sa*rv[0]*rv[0] - sa,       sa*rv[0]*rv[1] + ca*rv[2], sa*rv[0]*rv[2] - ca*rv[1],
        sa*rv[0]*rv[1] - ca*rv[2], sa*rv[1]*rv[1] - sa,       sa*rv[1]*rv[2] + ca*rv[0],
        sa*rv[0]*rv[2] + ca*rv[1], sa*rv[1]*rv[2] - ca*rv[0], sa*rv[2]*rv[2] - sa
      );
      // d(m*u*mt)/da = dm/da*u*mt + m*u*(dm/da)t
      scitbx::mat3<double> dm = rm*u_c*rmd.transpose() + rmd*u_c*rm.transpose();
      dm = unit_cell.fractionalization_matrix()*dm*
        unit_cell.fractionalization_matrix().transpose();
      for (int k = 0; k < 6; k++) {
        jt(ang->index(), index() + k) = dm[sym_acs[k]];
      }
    }
  }

  // shared u_iso
  void
  shared_u_iso
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    scalar_parameter *u_iso = reference();
    value = u_iso->value;

    if (jacobian_transpose == 0) {
      return;
    }
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.col(index()) = jt.col(u_iso->index());
  }

  // site
  void
  shared_site
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    site_parameter *site = reference();
    value = site->value;

    if (jacobian_transpose == 0) {
      return;
    }
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int i = 0; i < 3; i++) {
      jt.col(index() + i) = jt.col(site->index() + i);
    }
  }

  // shared fp
  void
    shared_fp
    ::linearise(uctbx::unit_cell const &unit_cell,
      sparse_matrix_type *jacobian_transpose)
  {
    scalar_parameter *fp = reference();
    value = fp->value;

    if (jacobian_transpose == 0) {
      return;
    }
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.col(index()) = jt.col(fp->index());
  }

  // shared fdp
  void
    shared_fdp
    ::linearise(uctbx::unit_cell const &unit_cell,
      sparse_matrix_type *jacobian_transpose)
  {
    scalar_parameter *fdp = reference();
    value = fdp->value;

    if (jacobian_transpose == 0) {
      return;
    }
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.col(index()) = jt.col(fdp->index());
  }

  // shared rotated u_star
  void
    shared_rotating_u_star
    ::linearise(uctbx::unit_cell const& unit_cell,
      sparse_matrix_type* jacobian_transpose)
  {
    const scalar_parameter* angles[] = {
      alpha(), beta(), gamma()
    };
    const scalar_parameter* sc = this->scale();
    scitbx::mat3<double> rmd[3];
    scitbx::mat3<double> rm = scitbx::math::euler_angles::xyz_matrix_rad(
      angles[0]->value, angles[1]->value, angles[2]->value, &rmd[0]),
      rmt = rm.transpose();
    const u_star_parameter* ref = reference();
    tensor_rank_2_t u_c =
      cctbx::adptbx::u_star_as_u_cart(unit_cell, ref->value);
    tensor_rank_2_t u_c_s = u_c * sc->value;
    scitbx::mat3<double> u_t = rm * u_c_s * rmt;
    // update the value
    value = cctbx::adptbx::u_cart_as_u_star(unit_cell,
      tensor_rank_2_t(
        u_t[0], u_t[4], u_t[8],
        u_t[1], u_t[2], u_t[5]));

    if (jacobian_transpose == 0) {
      return;
    };
    // convenience accessor array for the symmetric matrix...
    static const int sym_acs[] = { 0,4,8,1,2,5 };
    sparse_matrix_type& jt = *jacobian_transpose;
    // transforms for the jacobian values
    scitbx::mat3<double> jtm = unit_cell.fractionalization_matrix() * rm,
      jtm_t = jtm.transpose();
    for (int i = 0; i < jt.n_rows(); i++) {
      tensor_rank_2_t t;
      bool zero = true;
      for (int j = 0; j < 6; j++) {
        if ((t[j] = jt(i, ref->index() + j)) != 0) {
          zero = false;
        }
      }
      if (zero) {
        for (int j = 0; j < 6; j++) {
          jt(i, index() + j) = 0;
        }
      }
      else {
        tensor_rank_2_t t_c =
          cctbx::adptbx::u_star_as_u_cart(unit_cell, t);
        scitbx::mat3<double> x = jtm * (t_c * sc->value) * jtm_t;
        for (int j = 0; j < 6; j++) {
          jt(i, index() + j) = j < 3 ? x[sym_acs[j]] : 2*x[sym_acs[j]];
        }
      }
    }

    for (int j = 0; j < 3; j++) {
      if (!angles[j]->is_variable()) {
        continue;
      }
      scitbx::mat3<double> dm = rm * u_c_s * rmd[j].transpose() +
        rmd[j] * u_c_s * rmt;
      dm = unit_cell.fractionalization_matrix() * dm *
        unit_cell.fractionalization_matrix().transpose();
      for (int k = 0; k < 6; k++) {
        jt(angles[j]->index(), index() + k) = k <  3 ? dm[sym_acs[k]] : 2*dm[sym_acs[k]];
      }
    }

    // expansion/contraction
    if (sc->is_variable()) {
      for (int j = 0; j < 6; j++) {
        double v = value[j] / sc->value;
        jt(sc->index(), index()+j) = j < 3 ? v : 2*v;
      }
    }

  }

  void shared_rotating_u_star::validate() {
    if (scale()->value < 0) {
      scale()->value = 1e-3;
    }
    double tp = 2 * scitbx::constants::pi;
    //double tp = 360;
    while (alpha()->value < 0) {
      alpha()->value += tp;
    }
    while (alpha()->value > tp) {
      alpha()->value -= tp;
    }
    while (beta()->value < 0) {
      beta()->value += tp;
    }
    while (beta()->value > tp) {
      beta()->value -= tp;
    }
    while (gamma()->value < 0) {
      gamma()->value += tp;
    }
    while (gamma()->value > tp) {
      gamma()->value -= tp;
    }
  }


}}}
