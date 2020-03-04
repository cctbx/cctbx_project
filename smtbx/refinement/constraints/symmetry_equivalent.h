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
      sgtbx::rt_mx const &op);

    /// The original site in the asu this is a symmetry equivalent of
    asu_site_parameter *original() const {
      return dynamic_cast<asu_site_parameter *>(argument(0));
    }

    /// The symmetry moving original to this
    sgtbx::rt_mx const &motion() { return op; }

    virtual void linearise(uctbx::unit_cell const &unit_cell,
                           sparse_matrix_type *jacobian_transpose);
  private:
    class special_position_site_parameter *special_position;
    sgtbx::rt_mx op;
    // Transpose of the Jacobian of the transform x -> op*x
    af::ref<double, af::mat_grid> local_jt;
  };


}}}

#endif // GUARD
