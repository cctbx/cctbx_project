#include <smtbx/refinement/constraints/reparametrisation.h>
#include <iostream>
#include <scitbx/math/accumulators.h>
#include <scitbx/array_family/ref_reductions.h>

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

  // reparametrisation

  void reparametrisation
  ::analyse_variability() {
    /* Assign variability to each parameter.
     It also evaluates constant parameters once and for all.
     */
    variability_visitor var(unit_cell);
    accept(var);

    // Assign indices to parameters
    n_independents_ = n_intermediates_ = n_non_trivial_roots_ = 0;
    BOOST_FOREACH(parameter *p, all) {
      std::size_t s = p->size();
      if      (!p->is_variable())   n_intermediates_ += s;
      else if (p->is_independent()) n_independents_ += s;
      else if (p->is_root())        n_non_trivial_roots_ += s;
      else                          n_intermediates_ += s;
    }
    std::size_t i_independent = 0,
    i_intermediate = n_independents(),
    i_non_trivial_root = n_independents() + n_intermediates();
    BOOST_FOREACH(parameter *p, all) {
      std::size_t s = p->size();
      if      (!p->is_variable()) {
        p->set_index(i_intermediate);
        i_intermediate += s;
      }
      else if (!p->n_arguments()) {
        p->set_index(i_independent);
        i_independent += s;
      }
      else if (p->is_root()) {
        p->set_index(i_non_trivial_root);
        i_non_trivial_root += s;
      }
      else {
        p->set_index(i_intermediate);
        i_intermediate += s;
      }
    }

    // Initialise Jacobian transpose: [ dx_j/dx_i ]_ij
    /* The block of independent parameters is initialised to the identity matrix.
     Logically, it should be done in independent_xxxx_parameter::linearise,
     but it is more efficient to do it once and for all here.
     */
    sparse_matrix_type jt(n_independents(), n_components());
    for (std::size_t j=0; j<n_independents(); ++j) jt(j, j) = 1.;
    jacobian_transpose = jt;
  }

  void reparametrisation
  ::add(parameter *p) {
    typedef std::back_insert_iterator<std::vector<parameter *> >
    all_param_inserter_t;
    topologist<all_param_inserter_t> t(std::back_inserter(all));
    t.visit(p);
  }

  void reparametrisation
  ::finalise() {
    whiten(); // only time we need to call that explicitly
    analyse_variability();
  }

  reparametrisation::~reparametrisation() {
    BOOST_FOREACH(parameter *p, all) delete p;
  }

  void reparametrisation
  ::linearise() {
    // Initialise to zero Jacobian columns of intermediate and non trivial roots
    for (std::size_t j=n_independents(); j<n_components(); ++j) {
      jacobian_transpose.col(j).zero();
    }
    evaluator eval(unit_cell, &jacobian_transpose);
    accept(eval);
  }

  void reparametrisation
  ::apply_shifts(af::const_ref<double> const &shifts) {
    SMTBX_ASSERT(shifts.size() == n_independents());
    BOOST_FOREACH(parameter *p, all) {
      if (p->is_independent() && p->is_variable()) {
        double const *s = &shifts[p->index()];
        af::ref<double> x = p->components();
        for (std::size_t i=0; i<x.size(); ++i) x[i] += s[i];
      }
    }
  }

  double reparametrisation
  ::norm_of_independent_parameter_vector() {
    scitbx::math::accumulator::norm_accumulator<double> acc;
    BOOST_FOREACH(parameter *p, all) {
      if (p->is_independent() && p->is_variable()) {
        acc(af::sum_sq(p->components()));
      }
    }
    return acc.norm();
  }

  void reparametrisation
  ::store() {
    BOOST_FOREACH(parameter *p, all) {
      asu_parameter *cp = dynamic_cast<asu_parameter *> (p);
      if (cp) cp->store(unit_cell);
    }
  }

  void reparametrisation
  ::whiten() {
    BOOST_FOREACH(parameter *p, all) p->set_colour(white);
  }


}}}
