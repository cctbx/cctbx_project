#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

  // parameter

  parameter::~parameter() { delete[] arg; }

  bool parameter::is_variable() const { return variable; }

  void parameter::set_variable(bool f) { variable = f; }

  double *parameter::components() { return 0; }

 // independent_scalar_parameter

  std::size_t independent_scalar_parameter::size() const { return 1; }

  void independent_scalar_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {}

  double *independent_scalar_parameter::components() { return &value; }

  // independent_small_vector_parameter

  template class independent_small_vector_parameter<3>;
  template class independent_small_vector_parameter<6>;

  // site_parameter

  std::size_t site_parameter::size() const { return 3; }

  void site_parameter::store(uctbx::unit_cell const &unit_cell) const {
    scatterer->site = value;
  }

  // independent_site_parameter

  void independent_site_parameter::set_variable(bool f) {
    scatterer->flags.set_grad_site(f);
  }

  bool independent_site_parameter::is_variable() const {
    return scatterer->flags.grad_site();
  }

  void independent_site_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    value = scatterer->site;
  }

  double *independent_site_parameter::components() {
    return value.begin();
  }

  // ADP

  std::size_t cartesian_adp::size() const { return 6; }

  void cartesian_adp::store(uctbx::unit_cell const &unit_cell) const {
    scatterer->u_star = adptbx::u_cart_as_u_star(unit_cell, value);
  }

  // independent ADP

  void independent_cartesian_adp::set_variable(bool f) {
    if (f) scatterer->flags.set_use_u_aniso(true);
    scatterer->flags.set_grad_u_aniso(f);
  }

  bool independent_cartesian_adp::is_variable() const {
    return scatterer->flags.use_u_aniso() && scatterer->flags.grad_u_aniso();
  }

  void independent_cartesian_adp
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    value = adptbx::u_star_as_u_cart(unit_cell, scatterer->u_star);
  }

  double *independent_cartesian_adp::components() {
    return value.begin();
  }

}}}
