#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SYMMETRY_EQUIVALENT_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SYMMETRY_EQUIVALENT_H

#include <cctbx/sgtbx/rt_mx.h>
#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

  /// A symmetry equivalent site of a site in the asu.
  class symmetry_equivalent_site_parameter : public site_parameter
  {
  public:
    symmetry_equivalent_site_parameter(site_parameter *site,
                                       sgtbx::rt_mx const &op)
    : parameter(1),
      op(op), local_jt(3, 3)
    {
      set_arguments(site);
      scitbx::mat3<double> r_t = op.r().as_double().transpose();
      af::const_ref<double, af::mat_grid> rr_t(r_t.begin(), af::mat_grid(3,3));
      local_jt.assign_block(rr_t, 0, 0);
    }

    /// The original site in the asu this is a symmetry equivalent of
    asu_site_parameter *original() const {
      return dynamic_cast<asu_site_parameter *>(argument(0));
    }

    /// The symmetry moving original to this
    sgtbx::rt_mx const &motion() { return op; }

    virtual void linearise(uctbx::unit_cell const &unit_cell,
                           sparse_matrix_type *jacobian_transpose);
  private:
    sgtbx::rt_mx op;

    // Transpose of the Jacobian of the transform x -> op*x
    scitbx::sparse::matrix<double> local_jt;
  };


}}}

#endif // GUARD
