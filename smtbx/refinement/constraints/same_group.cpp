#include <smtbx/refinement/constraints/same_group.h>

namespace smtbx { namespace refinement { namespace constraints {

  index_range
  same_group_xyz
  ::component_indices_for(scatterer_type const *scatterer) const {
    for (int i=0; i < scatterers_.size(); i++) {
      if (scatterers_[i] == scatterer)
        return index_range(index() + i*3, 3);
    }
    return index_range();
  }

  index_range
  same_group_u_iso
  ::component_indices_for(scatterer_type const *scatterer) const {
    int offset = 0;
    for (int i=0; i < scatterers_.size(); i++) {
      if (scatterers_[i] == scatterer)
        return index_range(index() + i, 1);
    }
    return index_range();
  }

  index_range
  same_group_u_star
  ::component_indices_for(scatterer_type const *scatterer) const {
    for (int i=0; i < scatterers_.size(); i++) {
      if (scatterers_[i] == scatterer)
        return index_range(index() + i*6, 6);
    }
    return index_range();
  }

  void
  same_group_xyz
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    for (int i=0; i < scatterers_.size(); i++)  {
      if( scatterers_[i] == scatterer )  {
        output << scatterers_[i]->label << ".x,"
          << scatterers_[i]->label << ".y,"
          << scatterers_[i]->label << ".z,";
        return;
      }
    }
  }

  void
  same_group_u_iso
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    for (int i=0; i < scatterers_.size(); i++)  {
      if( scatterers_[i] == scatterer )  {
        output << scatterers_[i]->label << ".uiso,";
        return;
      }
    }
  }

  void
  same_group_u_star
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    for (int i=0; i < scatterers_.size(); i++)  {
      if( scatterers_[i] == scatterer )  {
        output << scatterers_[i]->label << ".u11,"
          << scatterers_[i]->label << ".u22,"
          << scatterers_[i]->label << ".u33,"
          << scatterers_[i]->label << ".u12,"
          << scatterers_[i]->label << ".u13,"
          << scatterers_[i]->label << ".u23,";
        return;
      }
    }
  }



  void
  same_group_xyz
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<6> *values =
      dynamic_cast<independent_small_vector_parameter<6> *>(argument(0));
    const double
      a_v = values->value[3], c_a = cos(a_v), s_a = sin(a_v),
      b_v = values->value[4], c_b = cos(b_v), s_b = sin(b_v),
      g_v = values->value[5], c_g = cos(g_v), s_g = sin(g_v);

    cart_t shift = cart_t(
      values->value[0], values->value[1], values->value[2]);

    // rotation matrix, for ref: http://en.wikipedia.org/wiki/Rotation_matrix
    scitbx::mat3<double> rm = alignment_matrix * scitbx::mat3<double>(
        c_b*c_g, -c_b*s_g,  s_b,
        s_a*s_b*c_g + c_a*s_g, -s_a*s_b*s_g+c_a*c_g, -s_a*c_b,
        -c_a*s_b*c_g + s_a*s_g, c_a*s_b*s_g + s_a*c_g, c_a*c_b
      );
    const scitbx::mat3<double> rm_t = rm.transpose();
    // derivative of the rotation matrix by angle
    const scitbx::mat3<double> rmd[3] = {
      alignment_matrix*scitbx::mat3<double>(
        0,                      0,                      0,
        c_a*s_b*c_g - s_a*s_g, -c_a*s_b*s_g - s_a*c_g, -c_a*c_b,
        s_a*s_b*c_g + c_a*s_g, -s_a*s_b*s_g + c_a*c_g, -s_a*c_b
      ),
      alignment_matrix*scitbx::mat3<double>(
        -c_g*s_b,     s_g*s_b,      c_b,
        s_a*c_g*c_b, -s_a*s_g*c_b,  s_a*s_b,
        -c_a*c_g*c_b, c_a*s_g*c_b, -c_a*s_b
      ),
      alignment_matrix*scitbx::mat3<double>(
        -c_b*s_g,               -c_b*c_g,               0,
        -s_a*s_b*s_g + c_a*c_g, -s_a*s_b*c_g - c_a*s_g, 0,
         c_a*s_b*s_g + s_a*c_g,  c_a*s_b*c_g - s_a*s_g, 0
      )
    };
    const frac_t grad_t[3] = {
      unit_cell.fractionalize(cart_t(1,0,0)),
      unit_cell.fractionalize(cart_t(0,1,0)),
      unit_cell.fractionalize(cart_t(0,0,1))
    };
    // calculate the rotation center
    af::shared<cart_t> co_s(scatterers_.size());
    cart_t rot_cnt(0,0,0);
    for (int i=0; i < scatterers_.size(); i++)  {
      co_s[i] = unit_cell.orthogonalize(
        dynamic_cast<site_parameter*>(argument(i+1))->value);
      rot_cnt += co_s[i];
    }
    rot_cnt = rot_cnt/scatterers_.size();
    scitbx::mat3<double> jtm = unit_cell.orthogonalization_matrix()*
      rm*unit_cell.fractionalization_matrix();
    for (int i=0; i < scatterers_.size(); i++) {
      // update site of i-th atoms
      co_s[i] -= rot_cnt;
      fx_s[i] = unit_cell.fractionalize(co_s[i]*rm + shift + rot_cnt);
      // derivatives
      if (jacobian_transpose != NULL) {
        sparse_matrix_type &jt = *jacobian_transpose;
        std::size_t const j_s = this->index() + i*3;
        // translation
        for (int j=0; j<3; j++) {
          for (int k=0; k<3; k++)
            jt(values->index()+j, j_s+k) = grad_t[j][k];
        }
        // transform the jacobian for shifts
        for (int j=0; j<jt.n_rows(); j++)  {
          cart_t t;
          for (int k=0; k<3; k++)
            t[k] = jt(j, argument(i+1)->index()+k);
          cart_t x = t*jtm;
          for (int k=0; k<3; k++)
            jt(j, j_s+k) = x[k];
        }
        // rotation
        for (int j=3; j<6; j++ )  {
          frac_t grad_f = unit_cell.fractionalize(co_s[i]*rmd[j-3]);
          for (int k=0; k<3; k++)
            jt(values->index()+j, j_s+k) = grad_f[k];
        }
      }
    }
  }

  void
  same_group_u_iso
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    for (int i=0; i < scatterers_.size(); i++) {
      u_isos[i] = dynamic_cast<scalar_parameter*>(argument(i+1))->value;
      if (jacobian_transpose != NULL) {
        sparse_matrix_type &jt = *jacobian_transpose;
        jt.col(index()) = jt.col(argument(i+1)->index());
      }
    }
  }

  void
  same_group_u_star
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<6> *values =
      dynamic_cast<independent_small_vector_parameter<6> *>(argument(0));
    const double
      a_v = values->value[3], c_a = cos(a_v), s_a = sin(a_v),
      b_v = values->value[4], c_b = cos(b_v), s_b = sin(b_v),
      g_v = values->value[5], c_g = cos(g_v), s_g = sin(g_v);

    scitbx::mat3<double> rm = alignment_matrix * scitbx::mat3<double>(
        c_b*c_g, -c_b*s_g,  s_b,
        s_a*s_b*c_g + c_a*s_g, -s_a*s_b*s_g+c_a*c_g, -s_a*c_b,
        -c_a*s_b*c_g + s_a*s_g, c_a*s_b*s_g + s_a*c_g, c_a*c_b
        );
    const scitbx::mat3<double> rm_t = rm.transpose();
    static const int sym_acs[] = {0,4,8,1,2,5};
    // transforms for the jacobian values
    scitbx::mat3<double> jtm = unit_cell.fractionalization_matrix()*
      rm*unit_cell.orthogonalization_matrix(), jtm_t = jtm.transpose();
    for (int i=0; i < scatterers_.size(); i++) {
      const tensor_rank_2_t u_c =
        cctbx::adptbx::u_star_as_u_cart(unit_cell,
        dynamic_cast<u_star_parameter*>(argument(i+1))->value);
      scitbx::mat3<double> u_t = rm*u_c*rm_t;
      u_stars[i] = cctbx::adptbx::u_cart_as_u_star(unit_cell,
        tensor_rank_2_t(
          u_t[sym_acs[0]], u_t[sym_acs[1]], u_t[sym_acs[2]],
          u_t[sym_acs[3]], u_t[sym_acs[4]], u_t[sym_acs[5]])
        );
      if (!jacobian_transpose) continue;
      sparse_matrix_type &jt = *jacobian_transpose;
      std::size_t const j_u = this->index() + i*6;
      // update jacobian for the rotation of ADP
      for (int j=0; j<jt.n_rows(); j++)  {
        tensor_rank_2_t t;
        bool zero = true;
        for (int k=0; k < 6; k++) {
          if ((t[k] = jt(j, argument(i+1)->index()+k)) != 0 )
            zero = false;
        }
        if (zero) {
          for (int k=0; k < 6; k++)
            jt(j, j_u+k) = 0;
        }
        else {
          scitbx::mat3<double> x = jtm*t*jtm_t;
          for (int k=0; k < 6; k++)
            jt(j, j_u+k) = x[sym_acs[k]];
        }
      }
    }
  }
}}}
