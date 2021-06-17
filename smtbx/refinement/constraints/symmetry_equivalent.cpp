#include <smtbx/refinement/constraints/symmetry_equivalent.h>
#include <smtbx/refinement/constraints/special_position.h>
#include <scitbx/array_family/versa_matrix.h>

namespace smtbx { namespace refinement { namespace constraints {

  // Symmetry equivalent site

  symmetry_equivalent_site_parameter
    ::symmetry_equivalent_site_parameter(site_parameter *site,
    sgtbx::rt_mx const &op)
    : parameter(1),
    special_position(dynamic_cast<special_position_site_parameter *>(site)),
    op(op)
  {
    set_arguments(site);
    if (special_position == 0) {
      scitbx::mat3<double> r_t = op.r().as_double().transpose();
      local_jt = af::ref<double, af::mat_grid>(r_t.begin(), af::mat_grid(3, 3));
    }
    else {
      af::const_ref<double, af::mat_grid> gm =
        special_position->get_site_constraints().gradient_sum_matrix();
      af::const_ref<double, af::mat_grid> rm =
        af::const_ref<double, af::mat_grid>(op.r().as_double().begin(), 3, 3);
      af::shared<double> res(3 * gm.n_rows());
      local_jt = af::ref<double, af::mat_grid>(res.begin(), gm.n_rows(), 3);
      af::multiply(gm, rm, local_jt);
    }
  }

  void symmetry_equivalent_site_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    site_parameter *site = original();

    value = op*site->value;

    if (!jacobian_transpose) {
      return;
    }
    sparse_matrix_type &jt = *jacobian_transpose;
    if (special_position == 0) {
      jt.assign_block(local_jt, site->index(), index());
    }
    else {
      jt.assign_block(local_jt,
        special_position->independent_params().index(), index());
    }
  }

}}}
