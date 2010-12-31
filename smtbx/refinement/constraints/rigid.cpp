#include <smtbx/refinement/constraints/rigid.h>

namespace smtbx { namespace refinement { namespace constraints {

  index_range
  rigid_group_base
  ::component_indices_for(scatterer_type const *scatterer) const {
    for (int i=0; i < scatterers_.size(); i++)  {
      if( scatterers_[i] == scatterer )
        return index_range(index() + 3*i, 3);
    }
    return index_range();
  }

  void
  rigid_group_base
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    for (int i=0; i < scatterers_.size(); i++)  {
      if( scatterers_[i] == scatterer )  {
        output << scatterers_[i]->label << ".x,";
        output << scatterers_[i]->label << ".y,";
        output << scatterers_[i]->label << ".z,";
        return;
      }
    }
  }

  void
  rigid_group_base
  ::InitCoordinates(uctbx::unit_cell const &unit_cell,
                    fractional<double> const& pivot)
  {
    if( crd_initialised )  return;
    original_pivot_crd = unit_cell.orthogonalize(pivot);
    for (int i=0; i < scatterers_.size(); i++)  {
      co_s[i] = unit_cell.orthogonalize(scatterers_[i]->site);
      rotation_center += co_s[i];
    }
    rotation_center = rotation_center/scatterers_.size();
    for (int i=0; i < scatterers_.size(); i++)
      co_s[i] = co_s[i] - rotation_center;
    crd_initialised = true;
  }


  // pivoted rotable...
  void
  pivoted_rotable_group
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    site_parameter
      *pivot = dynamic_cast<site_parameter *>(this->argument(0)),
      *pivot_neighbour = dynamic_cast<site_parameter *>(this->argument(1));
    scalar_parameter
      *azimuth = dynamic_cast<scalar_parameter *>(this->argument(2));

    const cart_t
      x_p = unit_cell.orthogonalize(pivot->value),
      x_pn = unit_cell.orthogonalize(pivot_neighbour->value),
      rv = (x_p - x_pn).normalize();
    const double
      angle = azimuth->value,
      ca = cos(angle),
      sa = sin(angle),
      t = 1.0-ca;
    // rotation matrix
    const scitbx::mat3<double> rm(
      t*rv[0]*rv[0] + ca,       t*rv[0]*rv[1] + sa*rv[2], t*rv[0]*rv[2] - sa*rv[1],
      t*rv[0]*rv[1] - sa*rv[2], t*rv[1]*rv[1] + ca,       t*rv[1]*rv[2] + sa*rv[0],
      t*rv[0]*rv[2] + sa*rv[1], t*rv[2]*rv[1] - sa*rv[0], t*rv[2]*rv[2] + ca
    );
    // derivative of the rotation matrix by angle
    const scitbx::mat3<double> rmd(
      sa*rv[0]*rv[0] - sa,       sa*rv[0]*rv[1] + ca*rv[2], sa*rv[0]*rv[2] - ca*rv[1],
      sa*rv[0]*rv[1] - ca*rv[2], sa*rv[1]*rv[1] - sa,       sa*rv[1]*rv[2] + ca*rv[0],
      sa*rv[0]*rv[2] + ca*rv[1], sa*rv[1]*rv[2] - ca*rv[0], sa*rv[2]*rv[2] - sa
    );
    // rotation happens around the geometrical center of the group...
    InitCoordinates(unit_cell, pivot->value);
    const cart_t center = rotation_center +
       unit_cell.orthogonalize(pivot->value) - original_pivot_crd;
    // Loop over the scatterers
    for (int i=0; i < scatterers_.size(); i++) {
      // update site of i-th scatterers
     fx_s[i] = unit_cell.fractionalize(co_s[i]*rm + center);

      // Derivatives
      if (!jacobian_transpose) continue;
      sparse_matrix_type &jt = *jacobian_transpose;
      std::size_t const j_s = this->index() + 3*i;

      // Riding
      for (int j=0; j<3; j++)
        jt.col(j_s + j) = jt.col(pivot->index() + j);

      // Rotating
      if (azimuth->is_variable()) {
        cart_t grad_c = co_s[i]*rmd;
        frac_t grad_f = unit_cell.fractionalize(grad_c);
        for (int j=0; j<3; j++)
          jt(azimuth->index(), j_s + j) = grad_f[j];
      }
    }
  }

  // spherical rotable expandable...
  void
  rotable_expandable_group
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    site_parameter
      *pivot = dynamic_cast<site_parameter *>(this->argument(0));
    scalar_parameter
      *size = dynamic_cast<scalar_parameter *>(argument(1));
    const scalar_parameter * angles[] = {
      dynamic_cast<scalar_parameter *>(argument(2)),
      dynamic_cast<scalar_parameter *>(argument(3)),
      dynamic_cast<scalar_parameter *>(argument(4))
    };
    const double
      a_v = angles[0]->value, c_a = cos(a_v), s_a = sin(a_v),
      b_v = angles[1]->value, c_b = cos(b_v), s_b = sin(b_v),
      g_v = angles[2]->value, c_g = cos(g_v), s_g = sin(g_v),
      size_value = size->value;
    // rotation matrix, for ref: http://en.wikipedia.org/wiki/Rotation_matrix
    const scitbx::mat3<double> rm(
      c_b*c_g, -c_b*s_g,  s_b,
      s_a*s_b*c_g + c_a*s_g, -s_a*s_b*s_g+c_a*c_g, -s_a*c_b,
      -c_a*s_b*c_g + s_a*s_g, c_a*s_b*s_g + s_a*c_g, c_a*c_b
    );
    // derivative of the rotation matrix by angle
    const scitbx::mat3<double> rmd[3] = {
      scitbx::mat3<double>(
        0,                      0,                      0,
        c_a*s_b*c_g - s_a*s_g, -c_a*s_b*s_g - s_a*c_g, -c_a*c_b,
        s_a*s_b*c_g + c_a*s_g, -s_a*s_b*s_g + c_a*c_g, -s_a*c_b
      ),
      scitbx::mat3<double>(
        -c_g*s_b,     s_g*s_b,      c_b,
        s_a*c_g*c_b, -s_a*s_g*c_b,  s_a*s_b,
        -c_a*c_g*c_b, c_a*s_g*c_b, -c_a*s_b
      ),
      scitbx::mat3<double>(
        -c_b*s_g,               -c_b*c_g,               0,
        -s_a*s_b*s_g + c_a*c_g, -s_a*s_b*c_g - c_a*s_g, 0,
         c_a*s_b*s_g + s_a*c_g,  c_a*s_b*c_g - s_a*s_g, 0
      )
    };
    // calculate the geometrical center
    InitCoordinates(unit_cell, pivot->value);
    const cart_t center = rotation_center +
       unit_cell.orthogonalize(pivot->value) - original_pivot_crd;
    // expansion/contraction happens from/to the center
    for (int i=0; i < scatterers_.size(); i++) {
      // update site of i-th atoms
      const cart_t len_vec = size_value*co_s[i];
      fx_s[i] = unit_cell.fractionalize(len_vec*rm + center);

      // Derivatives
      if (!jacobian_transpose) continue;
      sparse_matrix_type &jt = *jacobian_transpose;
      std::size_t const j_s = this->index() + 3*i;

      // Riding
      for (int j=0; j<3; j++)
        jt.col(j_s + j) = jt.col(pivot->index() + j);

      // Rotating
      for (int j=0; j<3; j++ )  {
        if (!angles[j]->is_variable())  continue;
        cart_t grad_c = len_vec*rmd[j];
        frac_t grad_f = unit_cell.fractionalize(grad_c);
        for (int k=0; k<3; k++)
          jt(angles[j]->index(), j_s + k) = grad_f[k];
      }

      // bond stretching
      if (size->is_variable())  {
        frac_t grad_f = unit_cell.fractionalize(co_s[i]*rm);
        for (int j=0; j<3; j++)
          jt(size->index(), j_s + j) = grad_f[j];
      }
    }
  }

}}}
