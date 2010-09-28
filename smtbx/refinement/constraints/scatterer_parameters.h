#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SCATTERER_PARAMETERS_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SCATTERER_PARAMETERS_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

/// The constraints::parameters corresponding to a single scatterer
struct scatterer_parameters
{
  typedef crystallographic_parameter::scatterer_type scatterer_type;

  scatterer_type const *scatterer;
  crystallographic_parameter *site, *occupancy, *u;

  scatterer_parameters() {}

  scatterer_parameters(scatterer_type const *scatterer)
    : scatterer(scatterer),
      site(0), occupancy(0), u(0)
  {}

  scatterer_parameters(scatterer_type const *scatterer,
                       crystallographic_parameter *site,
                       crystallographic_parameter *occupancy,
                       crystallographic_parameter *u)
    : scatterer(scatterer),
      site(site), occupancy(occupancy), u(u)
  {}

};


}}}
#endif // GUARD
