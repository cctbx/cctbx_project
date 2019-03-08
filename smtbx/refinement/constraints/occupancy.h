#ifndef SMTBX_REFINEMENT_CONSTRAINTS_OCCUPANCY_H
#define SMTBX_REFINEMENT_CONSTRAINTS_OCCUPANCY_H

#include <smtbx/refinement/constraints/reparametrisation.h>
#include <smtbx/refinement/constraints/affine.h>

namespace smtbx { namespace refinement { namespace constraints {

/** A specialisation of `affine_scalar_parameter` for the occupancy of
 *  a scatterer in the a.s.u.
 */
class affine_asu_occupancy_parameter : public affine_scalar_parameter,
                                       public virtual asu_occupancy_parameter

{
public:
  affine_asu_occupancy_parameter(scalar_parameter *u_0, double a_0,
                                 double b,
                                 scatterer_type *scatterer)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer),
    affine_scalar_parameter(u_0, a_0, b)
  {}

  affine_asu_occupancy_parameter(scalar_parameter *u_0, double a_0,
                                 scalar_parameter *u_1, double a_1,
                                 double b,
                                 scatterer_type *scatterer)
  : parameter(2),
    single_asu_scatterer_parameter(scatterer),
    affine_scalar_parameter(u_0, a_0, u_1, a_1, b)
  {}

  affine_asu_occupancy_parameter(af::shared<scalar_parameter *> const &u,
                                 af::shared<double> const &a,
                                 double b,
                                 scatterer_type *scatterer)
  : parameter(u.size()),
    single_asu_scatterer_parameter(scatterer),
    affine_scalar_parameter(u, a, b)
  {}
};

/** occupancy of one site depends on the occupancy of the other site
 */
class dependent_occupancy : public asu_occupancy_parameter
{
public:
/** if _as_one is true - this_occu = _multiplier*original_occu, else
    this_occu = _multiplier*(1 - original_occu)
 */
  dependent_occupancy(scalar_parameter *original_occu,
                double original_multiplier,
                double _multiplier,
                bool _as_one,
                scatterer_type *scatterer)
  : parameter(1),
    multiplier(_multiplier),
    original_multiplier(original_multiplier),
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
  double multiplier, original_multiplier;
  bool as_one;
};


}}}

#endif // GUARD
