#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SYMMETRY_EQUIVALENT_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SYMMETRY_EQUIVALENT_H

#include <cctbx/sgtbx/rt_mx.h>
#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

  /** Multiple inheritance would be the clean way to do it but that would
      necessarily involve virtual inheritance and with the provision of
      the Boost.Python hierarchy extension to come, we shied away from
      that solution and settled for the pragmatic macros that follow
   */

  #define SMTBX_CONSTRAINTS_SYMMETRY_EQUIVALENT_DECLARE_BOILER_PLATE           \
    virtual index_range                                                        \
    component_indices_for(                                                     \
      crystallographic_parameter::scatterer_type const *scatterer) const;      \
                                                                               \
  virtual void                                                                 \
  write_component_annotations_for(                                             \
    crystallographic_parameter::scatterer_type const *scatterer,               \
    std::ostream &output) const;                                               \
                                                                               \
  virtual void store(uctbx::unit_cell const &unit_cell) const;


  #define SMTBX_CONSTRAINTS_SYMMETRY_EQUIVALENT_DEFINE_BOILER_PLATE(klass)     \
    index_range klass                                                          \
    ::component_indices_for(                                                   \
      crystallographic_parameter::scatterer_type const *scatterer) const       \
    {                                                                          \
      throw SMTBX_ERROR("scatterer outside the asu");                          \
    }                                                                          \
                                                                               \
    void klass                                                                 \
    ::write_component_annotations_for(                                         \
      crystallographic_parameter::scatterer_type const *scatterer,             \
      std::ostream &output) const                                              \
    {                                                                          \
      throw SMTBX_ERROR("scatterer outside the asu");                          \
    }                                                                          \
                                                                               \
    void klass                                                                 \
    ::store(uctbx::unit_cell const &unit_cell) const                           \
    {                                                                          \
      throw SMTBX_ERROR("scatterer outside the asu");                          \
    }


  /// A symmetry equivalent site of a site in the asu.
  class symmetry_equivalent_site_parameter : public site_parameter
  {
  public:
    symmetry_equivalent_site_parameter(site_parameter *site,
                                       sgtbx::rt_mx const &op)
    : site_parameter(static_cast<scatterer_type *>(0), 1),
      op(op), local_jt(3, 3)
    {
      set_arguments(site);
      scitbx::mat3<double> r_t = op.r().as_double().transpose();
      af::const_ref<double, af::mat_grid> rr_t(r_t.begin(), af::mat_grid(3,3));
      local_jt.assign_block(rr_t, 0, 0);
    }

    SMTBX_CONSTRAINTS_SYMMETRY_EQUIVALENT_DECLARE_BOILER_PLATE

    /// The original site in the asu this is a symmetry equivalent of
    site_parameter *original() const {
      return (site_parameter *)this->argument(0);
    }

    /// The symmetry moving original to this
    sgtbx::rt_mx const &motion() { return op; }

    virtual std::size_t size() const;

    virtual void linearise(uctbx::unit_cell const &unit_cell,
                           sparse_matrix_type *jacobian_transpose);

  private:
    sgtbx::rt_mx op;

    // Transpose of the Jacobian of the transform x -> op*x
    scitbx::sparse::matrix<double> local_jt;

    fractional<double> value;
  };


}}}

#endif // GUARD
