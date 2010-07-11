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

}}}
