#ifndef SMTBX_REFINEMENT_CONSTRAINTS_U_ISO_DEPENDENT_U_ISO_H
#define SMTBX_REFINEMENT_CONSTRAINTS_U_ISO_DEPENDENT_U_ISO_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

/** u_iso parameter that is proportional to the u_iso of another scatterer
 */
class u_iso_proportional_to_pivot_u_iso : public asu_u_iso_parameter
{
public:
  u_iso_proportional_to_pivot_u_iso(scalar_parameter *pivot_u_iso,
                                   double multiplier,
                                   scatterer_type *scatterer)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer),
    multiplier(multiplier)
  {
    this->set_arguments(pivot_u_iso);
  }

  scalar_parameter *pivot_u_iso() const {
    return dynamic_cast<scalar_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  double multiplier;
};


}}}

#endif // GUARD
