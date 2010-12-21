#ifndef SMTBX_REFINEMENT_CONSTRAINTS_U_EQUALS_U_H
#define SMTBX_REFINEMENT_CONSTRAINTS_U_EQUALS_U_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

/** u_star parameter equals to the u_star of some other scatterer.
 */
class shared_u_star : public asu_u_star_parameter
{
public:
  shared_u_star(u_star_parameter *original_u,
             scatterer_type *scatterer)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer)
  {
    this->set_arguments(original_u);
  }

  u_star_parameter *original() const {
    return dynamic_cast<u_star_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

};

/** u_iso parameter equals to the u_iso of some other scatterer.
 */
class shared_u_iso : public asu_u_iso_parameter
{
public:
  shared_u_iso(scalar_parameter *original_u_iso,
             scatterer_type *scatterer)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer)
  {
    this->set_arguments(original_u_iso);
  }

  scalar_parameter *original() const {
    return dynamic_cast<scalar_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

};

/** site is shared by two scatterers
 */
class shared_site : public asu_site_parameter
{
public:
  shared_site(site_parameter *original_site,
             scatterer_type *scatterer)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer)
  {
    this->set_arguments(original_site);
  }

  site_parameter *original() const {
    return dynamic_cast<site_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

};

}}}

#endif // GUARD
