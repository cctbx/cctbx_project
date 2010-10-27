#include <smtbx/refinement/constraints/reparametrisation.h>
#include <iostream>

namespace smtbx { namespace refinement { namespace constraints {

  // parameter

  parameter::~parameter() { delete[] arg; }

  bool parameter::is_variable() const { return variable; }

  void parameter::set_variable(bool f) { variable = f; }

  // scalar parameter

  af::ref<double> scalar_parameter::components() {
    return af::ref<double>(&value, 1);
  }

  // independent_scalar_parameter

  void independent_scalar_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {}

  // independent_small_vector_parameter

  template class independent_small_vector_parameter<3>;
  template class independent_small_vector_parameter<6>;

  // single_scatterer_parameter

  asu_parameter::scatterer_sequence_type
  single_asu_scatterer_parameter::scatterers() const {
    return scatterer_sequence_type(&scatterer, 1);
  }

  index_range
  single_asu_scatterer_parameter
  ::component_indices_for(scatterer_type const *scatterer) const
  {
    return scatterer == this->scatterer ? index_range(index(), size())
                                        : index_range();
  }

  // site_parameter

  af::ref<double> site_parameter::components() { return value.ref(); }

  // asu_site_parameter

  void asu_site_parameter::set_variable(bool f) {
    scatterer->flags.set_grad_site(f);
  }

  bool asu_site_parameter::is_variable() const {
    return scatterer->flags.grad_site();
  }

  void asu_site_parameter
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    if (scatterer == this->scatterer) {
      output << scatterer->label << ".x,"
             << scatterer->label << ".y,"
             << scatterer->label << ".z,";
    }
  }

  void asu_site_parameter::store(uctbx::unit_cell const &unit_cell) const {
    scatterer->site = value;
  }

  // independent_site_parameter

  void independent_site_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {}

  // ADP

  af::ref<double> u_star_parameter::components() { return value.ref(); }

  // asu ADP

  void asu_u_star_parameter::set_variable(bool f) {
    if (f) scatterer->flags.set_use_u_aniso(true);
    scatterer->flags.set_grad_u_aniso(f);
  }

  bool asu_u_star_parameter::is_variable() const {
    return scatterer->flags.use_u_aniso() && scatterer->flags.grad_u_aniso();
  }

  void
  asu_u_star_parameter
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    if (scatterer == this->scatterer) {
      output << scatterer->label << ".u11,"
             << scatterer->label << ".u22,"
             << scatterer->label << ".u33,"
             << scatterer->label << ".u12,"
             << scatterer->label << ".u13,"
             << scatterer->label << ".u23,";
    }
  }

  void asu_u_star_parameter::store(uctbx::unit_cell const &unit_cell) const {
    scatterer->u_star = value;
  }

  // independent ADP

  void independent_u_star_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {}

  // Occupancy

  void asu_occupancy_parameter::set_variable(bool f) {
    scatterer->flags.set_grad_occupancy(f);
  }

  bool asu_occupancy_parameter::is_variable() const {
    return scatterer->flags.grad_occupancy();
  }

  void asu_occupancy_parameter
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    if (scatterer == this->scatterer) output << scatterer->label << ".occ,";
  }
  void asu_occupancy_parameter::store(uctbx::unit_cell const &unit_cell) const {
    scatterer->occupancy = value;
  }

  // independent Occupancy

  void independent_occupancy_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {}

  // u_iso

  void asu_u_iso_parameter::set_variable(bool f) {
    if (f) scatterer->flags.set_use_u_iso(true);
    scatterer->flags.set_grad_u_iso(f);
  }

  bool asu_u_iso_parameter::is_variable() const {
    return scatterer->flags.grad_u_iso();
  }

  void asu_u_iso_parameter
  ::write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const
  {
    if (scatterer == this->scatterer) output << scatterer->label << ".uiso,";
  }
  void asu_u_iso_parameter::store(uctbx::unit_cell const &unit_cell) const {
    scatterer->u_iso = value;
  }

  // independent u_iso

  void independent_u_iso_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {}

}}}
