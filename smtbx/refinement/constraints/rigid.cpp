#include <smtbx/refinement/constraints/rigid.h>

namespace smtbx { namespace refinement { namespace constraints {

  index_range
  pivoted_rotable_group
  ::component_indices_for(scatterer_type const *scatterer) const {
    for (int i=0; i < scatterers_.size(); i++)  {
      if( scatterers_[i] == scatterer )
        return index_range(index() + 3*i, 3);
    }
    return index_range();
  }

  void
  pivoted_rotable_group
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    for (int i=0; i < scatterers_.size(); i++)  {
      if( scatterers_[i] == scatterer )  {
        output << scatterers_[i]->label << ".x,";
        output << scatterers_[i]->label << ".y,";
        output << scatterers_[i]->label << ".z,";
      }
    }
  }

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

    cart_t
      x_p = unit_cell.orthogonalize(pivot->value),
      x_pn = unit_cell.orthogonalize(pivot_neighbour->value),
      rv = (x_p - x_pn).normalize();
    const double
      angle = azimuth->value,
      ca = cos(angle),
      sa = sin(angle),
      t = 1.0-ca;
    scitbx::mat3<double> rm(
      t*rv[0]*rv[0] + ca,       t*rv[0]*rv[1] + sa*rv[2], t*rv[0]*rv[2] - sa*rv[1],
      t*rv[0]*rv[1] - sa*rv[2], t*rv[1]*rv[1] + ca,       t*rv[1]*rv[2] + sa*rv[0],
      t*rv[0]*rv[2] + sa*rv[1], t*rv[2]*rv[1] - sa*rv[0], t*rv[2]*rv[2] + ca
    );
    scitbx::mat3<double> rmd(
      sa*rv[0]*rv[0] - sa,       sa*rv[0]*rv[1] + ca*rv[2], sa*rv[0]*rv[2] - ca*rv[1],
      sa*rv[0]*rv[1] - ca*rv[2], sa*rv[1]*rv[1] - sa,       sa*rv[1]*rv[2] + ca*rv[0],
      sa*rv[0]*rv[2] + ca*rv[1], sa*rv[1]*rv[2] - ca*rv[0], sa*rv[2]*rv[2] - sa
    );
    cart_t center;
    if( original_crd_initialised )  {
      for (int i=0; i < scatterers_.size(); i++)
        center += unit_cell.orthogonalize(scatterers_[i]->site);
    }
    else  {
      for (int i=0; i < scatterers_.size(); i++)  {
        co_s[i] = unit_cell.orthogonalize(scatterers_[i]->site);
        center += co_s[i];
      }
      original_crd_initialised = true;
    }
    center = center/scatterers_.size();
    // Loop over the scatterers
    for (int i=0; i < scatterers_.size(); i++) {
      // update site of i-th scatterers
     cx_s[i] = (co_s[i]-center)*rm + center;
     fx_s[i] = unit_cell.fractionalize(cx_s[i]);

      // Derivatives
      if (!jacobian_transpose) continue;
      sparse_matrix_type &jt = *jacobian_transpose;
      std::size_t const j_s = this->index() + 3*i;

      // Riding
      for (int j=0; j<3; j++)
        jt.col(j_s + j) = jt.col(pivot->index() + j);

      // Rotating
      if (azimuth->is_variable()) {
        cart_t grad_c = (co_s[i]-center)*rmd;
        frac_t grad_f = unit_cell.fractionalize(grad_c);
        for (int j=0; j<3; j++)
          jt(azimuth->index(), j_s + j) = grad_f[j];
      }
    }
  }

}}}
