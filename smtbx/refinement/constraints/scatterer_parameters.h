#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SCATTERER_PARAMETERS_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SCATTERER_PARAMETERS_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

/// The constraints::parameters corresponding to a single scatterer
struct scatterer_parameters
{
  parameter *site, *occupancy, *u;

  scatterer_parameters()
    : site(0), occupancy(0), u(0)
  {}

  scatterer_parameters(parameter *site, parameter *occupancy, parameter *u)
    : site(site), occupancy(occupancy), u(u)
  {}
};


}}}
#endif // GUARD
