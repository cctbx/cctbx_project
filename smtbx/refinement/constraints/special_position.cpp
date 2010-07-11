#include <smtbx/refinement/constraints/special_position.h>

namespace smtbx { namespace refinement { namespace constraints {

  /**** Sites ****/

  std::size_t special_position_site::size() const { return 3; }

  void special_position_site::linearise(uctbx::unit_cell const &unit_cell,
                                        sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<3> const &p = independent_params();

    site = site_constraints.all_params(p.value);

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    af::const_ref<double, af::mat_grid>
    compact_jt = site_constraints.gradient_sum_matrix();
    for (int j=0; j<3; ++j) {
      for (int i=0; i<site_constraints.n_independent_params(); ++i) {
        jt(p.index() + i, index() + j) = compact_jt(i, j);
      }
    }
  }

  void special_position_site::store(uctbx::unit_cell const &unit_cell) const {
    scatterer->site = site;
  }

  /**** ADP ****/

  std::size_t special_position_cartesian_adp::size() const { return 6; }

  void special_position_cartesian_adp
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<6> const &p = independent_params();

    u_cart = adp_constraints.all_params(p.value);

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    sparse_matrix_type const &jt_u_cart = adp_constraints.jacobian();
    for (int j=0; j<6; ++j) {
      for (sparse_matrix_type::const_row_iterator i=jt_u_cart.col(j).begin();
           i != jt_u_cart.col(j).end();
           ++i)
      {
        jt(j, i.index()) = *i;
      }
    }
  }

  void special_position_cartesian_adp
  ::store(uctbx::unit_cell const &unit_cell) const
  {
    scatterer->u_star = adptbx::u_cart_as_u_star(unit_cell, u_cart);
  }

}}}
