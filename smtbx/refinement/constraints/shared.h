#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SHARED_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SHARED_H

#include <smtbx/refinement/constraints/reparametrisation.h>
#include <smtbx/refinement/constraints/direction.h>

namespace smtbx { namespace refinement { namespace constraints {

/** u_star parameter equals to the reference u_star
 */
class shared_u_star : public asu_u_star_parameter {
public:
  shared_u_star(scatterer_type *scatterer,
    u_star_parameter *reference)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer)
  {
    this->set_arguments(reference);
  }

  u_star_parameter *reference() const {
    return dynamic_cast<u_star_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

/** u_star equals to reference u_star but the two are related by a
rotation around a vector defined by direction_from and direction_to
*/
class shared_rotated_u_star : public asu_u_star_parameter {
  direction_base* direction_;
public:
  shared_rotated_u_star(scatterer_type *scatterer,
    u_star_parameter *reference,
    direction_base *direction,
    independent_scalar_parameter *angle)
  : parameter(2),
    single_asu_scatterer_parameter(scatterer),
    direction_(direction)
  {
    set_arguments(reference, angle);
  }

  u_star_parameter *reference() const {
    return dynamic_cast<u_star_parameter *>(this->argument(0));
  }

  direction_base *direction() const {  return direction_;  }

  independent_scalar_parameter *angle() const {
    return dynamic_cast<independent_scalar_parameter *>(this->argument(1));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

/** u_iso parameter equals to the reference u_iso/scalar value.
 */
class shared_u_iso : public asu_u_iso_parameter {
public:
  shared_u_iso(scatterer_type *scatterer,
    scalar_parameter *reference)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer)
  {
    this->set_arguments(reference);
  }

  scalar_parameter *reference() const {
    return dynamic_cast<scalar_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

/** site is shared by two scatterers
 */
class shared_site : public asu_site_parameter {
public:
  shared_site(scatterer_type *scatterer,
    site_parameter *reference)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer)
  {
    this->set_arguments(reference);
  }

  site_parameter *reference() const {
    return dynamic_cast<site_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

}}}

#endif // GUARD
