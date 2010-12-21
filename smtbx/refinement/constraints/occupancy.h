#ifndef SMTBX_REFINEMENT_CONSTRAINTS_OCCUPANCY_H
#define SMTBX_REFINEMENT_CONSTRAINTS_OCCUPANCY_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

/** occupancy of one site depends on the occupancy of the other site
 */
class dependent_occupancy : public asu_occupancy_parameter
{
public:
/** if _as_one is true - this_occu = _multiplier*original_occu, else
    this_occu = _multiplier*(1 - original_occu)
 */
  dependent_occupancy(scalar_parameter *original_occu,
                double _multiplier,
                bool _as_one,
                scatterer_type *scatterer)
  : parameter(1),
    multiplier(_multiplier),
    as_one(_as_one),
    single_asu_scatterer_parameter(scatterer)
  {
    this->set_arguments(original_occu);
  }

  scalar_parameter *reference() const {
    return dynamic_cast<scalar_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
  double multiplier;
  bool as_one;
};


}}}

#endif // GUARD
