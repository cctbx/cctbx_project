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
    site_parameter *pivot = (site_parameter *)this->argument(0),
                   *pivot_neighbour = (site_parameter *)this->argument(1);
    independent_scalar_parameter *azimuth, *length;
    site_parameter *stagger;
    if (staggered) stagger = (site_parameter *)this->argument(2);
    else           azimuth = (independent_scalar_parameter *)this->argument(2);
    length  = (independent_scalar_parameter *)this->argument(3);

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
      cos_phi = 1.;
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

  // X-CH2-Y

  void secondary_ch2_sites::linearise(uctbx::unit_cell const &unit_cell,
                                      sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter *pivot             = (site_parameter *)argument(0),
                   *pivot_neighbour_0 = (site_parameter *)argument(1),
                   *pivot_neighbour_1 = (site_parameter *)argument(2);
    independent_scalar_parameter
    *length = (independent_scalar_parameter *)argument(3);
    angle_starting_tetrahedral
    *h_c_h = (angle_starting_tetrahedral *)argument(4);

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
    double l = length->value, theta = h_c_h->value;

    // Hydrogen sites
    double c = std::cos(theta/2), s = std::sin(theta/2);
    af::tiny<cart_t, 2> u_h(c*e0 + s*e1, c*e0 - s*e1);
    for (int k=0; k<2; ++k) x_h[k] = x_p + l*u_h[k];

    // Derivatives
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    af::tiny<std::size_t, 2> j_h(index(), index() + 3);

    // Riding
    for (int k=0; k<2; ++k) for (int i=0; i<3; ++i) {
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
    if (h_c_h->is_variable()) {
      af::tiny<cart_t, 2> grad_c(l/2*(-s*e0 + c*e1), l/2*(-s*e0 - c*e1));
      for (int k=0; k<2; ++k) {
        frac_t grad_f = unit_cell.fractionalize(grad_c[k]);
        for (int i=0; i<3; ++i) jt(h_c_h->index(), j_h[k] + i) = grad_f[i];
      }
    }
  }

  /***    H
          |
       X0-C-X1
          |
          X2
   */
  void tertiary_ch_site::linearise(uctbx::unit_cell const &unit_cell,
                                   sparse_matrix_type *jacobian_transpose)
  {
    using namespace constants;
    site_parameter *pivot = (site_parameter *)argument(0);
    af::tiny<site_parameter *, 3> pivot_neighbour;
    for (int k=0; k<3; ++k) {
      pivot_neighbour[k] = (site_parameter *)argument(k+1);
    }
    independent_scalar_parameter
    *length = (independent_scalar_parameter *)argument(4);

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
    site_parameter *pivot = (site_parameter *)argument(0);
    site_parameter *pivot_neighbour_0 = (site_parameter *)argument(1),
                   *pivot_neighbour_1 = (site_parameter *)argument(2);
    independent_scalar_parameter
    *length = (independent_scalar_parameter *)argument(3);

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
    site_parameter *pivot = (site_parameter *)argument(0);
    site_parameter *pivot_neighbour = (site_parameter *)argument(1);
    site_parameter *pivot_neighbour_substituent = (site_parameter *)argument(2);
    independent_scalar_parameter
    *length = (independent_scalar_parameter *)argument(3);

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
    for (int k=0; k<2; ++k) x_h[k] = p[k] + l*u_h[k];

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
    site_parameter *pivot = (site_parameter *)argument(0);
    site_parameter *pivot_neighbour = (site_parameter *)argument(1);
    independent_scalar_parameter
    *length = (independent_scalar_parameter *)argument(2);

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

}}}
