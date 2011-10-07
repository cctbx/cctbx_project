#include <smtbx/refinement/constraints/geometrical_hydrogens.h>
#include <scitbx/sparse/io.h>

#include <boost/lambda/lambda.hpp>
#include <scitbx/array_family/tiny_reductions.h>

namespace smtbx { namespace refinement { namespace constraints {

  //*** Base class ***

  template <int n_hydrogens>
  index_range
  geometrical_hydrogen_sites<n_hydrogens>
  ::component_indices_for(scatterer_type const *scatterer) const
  {
    using boost::lambda::_1;
    boost::optional<std::size_t>
    i_sc = af::first_index(hydrogen, _1 == scatterer);
    return i_sc ? index_range(index() + 3*(*i_sc), 3)
                : index_range();
  }

  template <int n_hydrogens>
  void geometrical_hydrogen_sites<n_hydrogens>
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    using boost::lambda::_1;
    boost::optional<std::size_t>
    i_sc = af::first_index(hydrogen, _1 == scatterer);
    if (i_sc) {
      scatterer_type const *h = hydrogen[*i_sc];
      output << h->label << ".x,";
      output << h->label << ".y,";
      output << h->label << ".z,";
    }
  }

  //*** CH3, NH2, OH ***

  template <int n_hydrogens, bool staggered>
  void
  terminal_tetrahedral_xhn_sites<n_hydrogens, staggered>
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter
    *pivot           = dynamic_cast<site_parameter *>(this->argument(0)),
    *pivot_neighbour = dynamic_cast<site_parameter *>(this->argument(1));
    scalar_parameter *azimuth, *length;
    site_parameter *stagger;
    if (staggered) stagger = dynamic_cast<site_parameter *>(this->argument(2));
    else           azimuth = dynamic_cast<scalar_parameter *>(this->argument(2));
    length  = dynamic_cast<scalar_parameter *>(this->argument(3));

    // Local frame
    cart_t x_p = unit_cell.orthogonalize(pivot->value),
           x_pn = unit_cell.orthogonalize(pivot_neighbour->value);
    if (staggered) {
      cart_t x_s = unit_cell.orthogonalize(stagger->value);
      e_zero_azimuth = x_s - x_pn;
    }
    af::tiny<cart_t, 3>
    e = scitbx::math::orthonormal_basis(x_p - x_pn, e_zero_azimuth);

    double l = length->value;
    double cos_phi, sin_phi;
    if (!staggered) {
      double phi = azimuth->value;
      cos_phi = std::cos(phi);
      sin_phi = std::sin(phi);
    }
    else {
      cos_phi = -1.;
      sin_phi = 0.;
    }

    // Loop over the Hydrogen atoms
    for (int k=0; k < n_hydrogens; ++k) {

      // Cosine and Sine of the azimutal angle of the k-th Hydrogen
      /* Mathematica:
       Table[TrigExpand[Cos[\[Phi] + n Pi/3]], {n, {2, 4}}]
       Table[TrigExpand[Sin[\[Phi] + n Pi/3]], {n, {2, 4}}]
       */
      double c, s;
      switch (k) {
        case 0:
          // 1st Hydrogen: azimuthal angle = phi
          c = cos_phi;
          s = sin_phi;
          break;
        case 1:
          // 2nd Hydrogen: azimuthal angle = phi + 2pi/3
          c = -0.5        *cos_phi - half_sqrt_3*sin_phi;
          s =  half_sqrt_3*cos_phi -         0.5*sin_phi;
          break;
        case 2:
          // 3rd Hydrogen: azimuthal angle = phi + 4pi/3
          c = -0.5        *cos_phi + half_sqrt_3*sin_phi;
          s = -half_sqrt_3*cos_phi -         0.5*sin_phi;
        default:
          break;
      }

      // Site of k-th Hydrogen
      cart_t u = sin_tetrahedral_angle*(c*e[1] + s*e[2]) + e[0]/3.;
      this->x_h[k] = x_p + l*u;

      // Derivatives
      if (!jacobian_transpose) continue;
      sparse_matrix_type &jt = *jacobian_transpose;
      std::size_t const j_h = this->index() + 3*k;

      // Riding
      for (int i=0; i<3; ++i) {
        jt.col(j_h + i) = jt.col(pivot->index() + i);
      }

      /** We take advantage of the fact that azimuth and length are
          independent variables. So jt.col(azimuth->index()) is either
          zero or is a column of the identity matrix.
       */

      // Rotation
      if (!staggered && azimuth->is_variable()) {
        cart_t grad_c = l*sin_tetrahedral_angle*(-s*e[1] + c*e[2]);
        frac_t grad_f = unit_cell.fractionalize(grad_c);
        for (int i=0; i<3; ++i) jt(azimuth->index(), j_h + i) = grad_f[i];
      }

      // Bond stretching
      if (length->is_variable()) {
        frac_t grad_f = unit_cell.fractionalize(u);
        for (int i=0; i<3; ++i) jt(length->index(), j_h + i) = grad_f[i];
      }
    }
  }

  template class terminal_tetrahedral_xhn_sites<1, /*staggered=*/false>;
  template class terminal_tetrahedral_xhn_sites<2, /*staggered=*/false>;
  template class terminal_tetrahedral_xhn_sites<3, /*staggered=*/false>;

  template class terminal_tetrahedral_xhn_sites<1, /*staggered=*/true>;
  template class terminal_tetrahedral_xhn_sites<2, /*staggered=*/true>;
  template class terminal_tetrahedral_xhn_sites<3, /*staggered=*/true>;

  // angle parameter

  void angle_parameter::linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose)
  {
    site_parameter *sites[] = {
      dynamic_cast<site_parameter *>(argument(0)),
      dynamic_cast<site_parameter *>(argument(1)),
      dynamic_cast<site_parameter *>(argument(2))
    };
    cart_t crds[] = {
      unit_cell.orthogonalize(sites[0]->value),
      unit_cell.orthogonalize(sites[1]->value),
      unit_cell.orthogonalize(sites[2]->value)
    };
    value = (crds[0]-crds[1]).angle(crds[2]-crds[1]);
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    // http://salilab.org/modeller/8v0/manual/node248.html
    cart_t ij = (crds[0] - crds[1]);
    cart_t kj = (crds[2] - crds[1]);
    const double ij_l = ij.length(), kj_l = kj.length();
    ij = ij.normalize();
    kj = kj.normalize();
    const double ca = ij*kj;
    //if( std::abs(ca) >= 1.0-1e-16 )
    const double oos = 1./std::sqrt(1-ca*ca);
    cart_t grad[3];
    grad[0] = (ij*ca - kj)*oos/ij_l;  // d(angle)/d(left)
    grad[2] = (kj*ca - ij)*oos/kj_l;  // d(angle)/d(right)
    grad[1] = -(grad[0] + grad[2]);   // d(angle)/d(center)
    for (int i=0; i < 3; i++)  {
      frac_t grad_f = unit_cell.fractionalize(grad[i]);
      for (int j=0; j < 3; j++) {
        jt(sites[i]->index()+j, index()) = grad_f[j];
      }
    }
  }

  // X-CH2-Y

  void secondary_xh2_sites::linearise(uctbx::unit_cell const &unit_cell,
                                      sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter
    *pivot             = dynamic_cast<site_parameter *>(argument(0)),
    *pivot_neighbour_0 = dynamic_cast<site_parameter *>(argument(1)),
    *pivot_neighbour_1 = dynamic_cast<site_parameter *>(argument(2));
    scalar_parameter *length = dynamic_cast<scalar_parameter *>(argument(3));
    scalar_parameter *h_c_h = dynamic_cast<scalar_parameter *>(argument(4));

    bool is_angle_param = dynamic_cast<angle_parameter *>(argument(4)) != 0;
    // Local frame
    /* (C, e0, e1) is the bisecting plane of the angle X-C-Y
        with e0 bisecting X-C-Y
     */
    cart_t x_p  = unit_cell.orthogonalize(pivot->value);
    cart_t x_pn_0 = unit_cell.orthogonalize(pivot_neighbour_0->value),
           x_pn_1 = unit_cell.orthogonalize(pivot_neighbour_1->value);
    cart_t u_pn_0 = (x_p - x_pn_0).normalize(),
           u_pn_1 = (x_p - x_pn_1).normalize();
    cart_t e0 = (u_pn_1 + u_pn_0).normalize();
    cart_t e2 = (u_pn_1 - u_pn_0).normalize();
    cart_t e1 = e2.cross(e0);
    double l = length->value,
      theta = h_c_h->value,
      cos_phi = 0;
    const double theta_to_phi_const = 0.0698;
    // calculate new theta/2 = 0.9678 + 0.0698 cos phi as in shelxl
    if (is_angle_param) {
      cos_phi = u_pn_0*u_pn_1;
      theta = 2*(0.9678 + theta_to_phi_const*cos_phi);
    }
    // Hydrogen sites
    double c = std::cos(theta/2), s = std::sin(theta/2);
    af::tiny<cart_t, 2> u_h(c*e0 + s*e1, c*e0 - s*e1);
    for (int k=0; k<2; ++k) x_h[k] = x_p + l*u_h[k];

    // Derivatives
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    af::tiny<std::size_t, 2> j_h(index(), index() + 3);

    // Riding
    for (int k=0; k<2; ++k) {
      for (int i=0; i<3; ++i)
        jt.col(j_h[k] + i) = jt.col(pivot->index() + i);
    }

    // Bond stretching
    if (length->is_variable()) {
      for (int k=0; k<2; ++k) {
        frac_t grad_f = unit_cell.fractionalize(u_h[k]);
        for (int i=0; i<3; ++i) jt(length->index(), j_h[k] + i) = grad_f[i];
      }
    }

    // H-C-H flapping
    if (h_c_h->is_variable() && !is_angle_param) {
      af::tiny<cart_t, 2> grad_c(l/2*(-s*e0 + c*e1), l/2*(-s*e0 - c*e1));
      for (int k=0; k<2; ++k) {
        frac_t grad_f = unit_cell.fractionalize(grad_c[k]);
        for (int i=0; i<3; ++i) jt(h_c_h->index(), j_h[k] + i) = grad_f[i];
      }
    }
    else if (is_angle_param) {
      /*
      d(cos(theta/2))/d(phi) = d(cos(theta/2)/d(theta)*d(theta)/d(phi) =
        sin(theta/2)*theta_to_phi_const*sin(phi)
      d(sin(theta/2))/d(phi) = d(sin(theta/2)/d(theta)*d(theta)/d(phi) =
        cos(theta/2)*theta_to_phi_const*cos(phi)
      */
      double sin_phi = std::sqrt(1-cos_phi*cos_phi),
        k = l*theta_to_phi_const/2;
      af::tiny<cart_t, 2> grad_c(
        k*(s*sin_phi*e0 + c*cos_phi*e1),
        k*(s*sin_phi*e0 - c*cos_phi*e1));
      for (int i=0; i < 2; i++) {
        frac_t grad_f = unit_cell.fractionalize(grad_c[i]);
        for (int j=0; j < 3; j++)  {
          jt(j_h[i]+j, h_c_h->index()) = grad_f[j];
        }
      }
    }
  }

  /***    H
          |
       X0-C-X1
          |
          X2
   */
  void tertiary_xh_site::linearise(uctbx::unit_cell const &unit_cell,
                                   sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter *pivot = dynamic_cast<site_parameter *>(argument(0));
    af::tiny<site_parameter *, 3> pivot_neighbour;
    for (int k=0; k<3; ++k) {
      pivot_neighbour[k] = dynamic_cast<site_parameter *>(argument(k+1));
    }
    scalar_parameter *length = dynamic_cast<scalar_parameter *>(argument(4));

    // Local frame
    cart_t x_p = unit_cell.orthogonalize(pivot->value);
    af::tiny<cart_t, 3> u_cn; // Directions C->Xi
    for (int k=0; k<3; ++k) {
      cart_t x = unit_cell.orthogonalize(pivot_neighbour[k]->value);
      u_cn[k] = (x_p - x).normalize();
    }
    cart_t u = u_cn[0] - u_cn[1], v = u_cn[1] - u_cn[2];
    cart_t e0 = u.cross(v).normalize();
    if (e0*(u_cn[0] + u_cn[1] + u_cn[2]) < 0) e0 = -e0;
    double l = length->value;

    // Hydrogen site
    x_h[0] = x_p + l*e0;

    // Derivatives
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    std::size_t j_h = index();

    // Riding
    for (int i=0; i<3; ++i) {
      jt.col(j_h + i) = jt.col(pivot->index() + i);
    }

    // Bond stretching
    if (length->is_variable()) {
      frac_t grad_f = unit_cell.fractionalize(e0);
      for (int i=0; i<3; ++i) jt(length->index(), j_h + i) = grad_f[i];
    }
  }

  /// aromatic or amide Y-XH-Z
  void secondary_planar_xh_site::linearise(uctbx::unit_cell const &unit_cell,
                                           sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter *pivot = dynamic_cast<site_parameter *>(argument(0));
    site_parameter
    *pivot_neighbour_0 = dynamic_cast<site_parameter *>(argument(1)),
    *pivot_neighbour_1 = dynamic_cast<site_parameter *>(argument(2));
    scalar_parameter *length = dynamic_cast<scalar_parameter *>(argument(3));

    // Local frame
    cart_t x_p = unit_cell.orthogonalize(pivot->value);
    cart_t
    u_yx = (x_p - unit_cell.orthogonalize(pivot_neighbour_0->value)).normalize(),
    u_zx = (x_p - unit_cell.orthogonalize(pivot_neighbour_1->value)).normalize();
    cart_t e0 = (u_yx + u_zx).normalize();
    double l = length->value;

    // Hydrogen site
    x_h[0] = x_p + l*e0;

    // Jacobian
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    std::size_t j_h = index();

    // Riding
    for (int i=0; i<3; ++i) jt.col(j_h + i) = jt.col(pivot->index() + i);

    // Bond stretching
    if (length->is_variable()) {
      frac_t grad_f = unit_cell.fractionalize(e0);
      for (int i=0; i<3; ++i) jt(length->index(), j_h + i) = grad_f[i];
    }
  }

  // Terminal Z-Y=XH2
  void terminal_planar_xh2_sites::
  linearise(uctbx::unit_cell const &unit_cell,
            sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter *pivot = dynamic_cast<site_parameter *>(argument(0));
    site_parameter
    *pivot_neighbour = dynamic_cast<site_parameter *>(argument(1));
    site_parameter
    *pivot_neighbour_substituent = dynamic_cast<site_parameter *>(argument(2));
    scalar_parameter *length = dynamic_cast<scalar_parameter *>(argument(3));

    // Local frame
    cart_t p = unit_cell.orthogonalize(pivot->value);
    cart_t y = unit_cell.orthogonalize(pivot_neighbour->value);
    cart_t z = unit_cell.orthogonalize(pivot_neighbour_substituent->value);
    cart_t e0 = (p - y).normalize();
    cart_t u_yz = z - y;
    cart_t e1 = (e0 - 1/(e0*u_yz) * u_yz).normalize();
    double l = length->value;

    // Hydrogen sites
    af::tiny<cart_t, 2> u_h(0.5*e0 + half_sqrt_3*e1,
                            0.5*e0 - half_sqrt_3*e1);
    for (int k=0; k<2; ++k) x_h[k] = p + l*u_h[k];

    // Jacobian
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    std::size_t j_h = index();

    // Riding
    for (int k=0; k<2; ++k) for (int i=0; i<3; ++i) {
      jt.col(j_h + 3*k + i) = jt.col(pivot->index() + i);
    }

    // Bond stretching
    if (length->is_variable()) {
      for (int k=0; k<2; ++k) {
        frac_t grad_f = unit_cell.fractionalize(u_h[k]);
        for (int i=0; i<3; ++i) jt(length->index(), j_h + 3*k + i) = grad_f[i];
      }
    }
  }


  // Acetylenic X-CH
  void terminal_linear_ch_site
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter *pivot = dynamic_cast<site_parameter *>(argument(0));
    site_parameter *pivot_neighbour = dynamic_cast<site_parameter *>(argument(1));
    scalar_parameter *length = dynamic_cast<scalar_parameter *>(argument(2));

    // Local frame
    cart_t p = unit_cell.orthogonalize(pivot->value);
    cart_t x = unit_cell.orthogonalize(pivot_neighbour->value);
    cart_t e0 = (p - x).normalize();
    double l = length->value;

    // Hydrogen site
    x_h[0] = p + l*e0;

    // Jacobian
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    std::size_t j_h = index();

    // Riding
    for (int i=0; i<3; ++i) jt.col(j_h + i) = jt.col(pivot->index() + i);

    // Bond stretching
    if (length->is_variable()) {
      frac_t grad_f = unit_cell.fractionalize(e0);
      for (int i=0; i<3; ++i) jt(length->index(), j_h + i) = grad_f[i];
    }
  }

  // boron cage B(n)-H
  void polyhedral_bh_site
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter *pivot = dynamic_cast<site_parameter *>(argument(0));
    scalar_parameter *length = dynamic_cast<scalar_parameter *>(argument(1));
    const cart_t c = unit_cell.orthogonalize(pivot->value);
    cart_t p(0,0,0);
    for (int i=2; i<n_arguments(); i++)  {
      p = p + (unit_cell.orthogonalize(
        dynamic_cast<site_parameter *>(argument(i))->value) - c).normalize();
    }
    p = -p.normalize();
    double l = length->value;
    // Hydrogen site
    x_h[0] = c + p*l;

    // Jacobian
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    std::size_t j_h = index();

    // Riding
    for (int i=0; i<3; i++)
      jt.col(j_h + i) = jt.col(pivot->index() + i);

    // Bond stretching
    if (length->is_variable()) {
      frac_t grad_f = unit_cell.fractionalize(p);
      for (int i=0; i<3; i++)
        jt(length->index(), j_h + i) = grad_f[i];
    }
  }

}}}
