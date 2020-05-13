#include <smtbx/refinement/constraints/reparametrisation.h>
namespace smtbx { namespace refinement { namespace constraints {

class scalar_scaled_u_star_parameter : public asu_u_star_parameter
{
  const scitbx::sym_mat3<double> u_star_ref_;
public:
  scalar_scaled_u_star_parameter(
    independent_scalar_parameter *scalar,
    scatterer_type *scatterer)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer),
    u_star_ref_(scatterer->u_star)
  {
    set_arguments(scalar);
  }

  scitbx::sym_mat3<double> const *reference() {
    return &u_star_ref_;
  }

  independent_scalar_parameter *scalar() const {
    return dynamic_cast<independent_scalar_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

};

class scalar_scaled_u_iso_parameter : public asu_u_iso_parameter
{
  const double u_iso_ref_;
public:
  scalar_scaled_u_iso_parameter(
    independent_scalar_parameter *scalar,
    scatterer_type *scatterer)
  : parameter(1),
    single_asu_scatterer_parameter(scatterer),
    u_iso_ref_(scatterer->u_iso)
  {
    set_arguments(scalar);
  }

  double const reference() {
    //this returns by value (unlike other reference() member functions) so we
    //can expose it in python, which wouldn't work if it returned double*
    return u_iso_ref_;
  }

  independent_scalar_parameter *scalar() const {
    return dynamic_cast<independent_scalar_parameter *>(this->argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

};


}}}
