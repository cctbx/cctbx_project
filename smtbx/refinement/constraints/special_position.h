#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SPECIAL_POSITION_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SPECIAL_POSITION_H

#include <scitbx/array_family/small.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/site_symmetry.h>
#include <smtbx/refinement/constraints/reparametrisation.h>
#include <smtbx/import_cctbx.h>

namespace smtbx { namespace refinement { namespace constraints {

/// Site constrained to be on a special position
/** Parameter components are the fractional coordinates */
class special_position_site_parameter : public site_parameter
{
public:
  special_position_site_parameter(sgtbx::site_symmetry_ops const &site_symmetry,
                        scatterer_type *scatterer)
    : site_parameter(scatterer, 1),
      site_constraints(site_symmetry.site_constraints())
  {
    value = site_symmetry.special_op()*scatterer->site;
    set_arguments(new independent_small_vector_parameter<3>(
      site_constraints.independent_params(value),
      scatterer->flags.grad_site()));
  }

  independent_small_vector_parameter<3> const &independent_params() {
    return *(independent_small_vector_parameter<3> *)argument(0);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

private:
  sgtbx::site_constraints<double> site_constraints;
};


/// Anisotropic displacement constrained by the symmetry of a special position
/** Parameter components are those of the tensor in fractional coordinates
 */
class special_position_u_star_parameter : public u_star_parameter
{
public:
  typedef sgtbx::tensor_rank_2::constraints<double>
          adp_constraints_t;

  special_position_u_star_parameter(sgtbx::site_symmetry_ops const &site_symmetry,
                                    scatterer_type *scatterer)
    : u_star_parameter(scatterer, 1),
      adp_constraints(site_symmetry.adp_constraints())
  {
    value = site_symmetry.average_u_star(scatterer->u_star);
    set_arguments(new independent_small_vector_parameter<6>(
      adp_constraints.independent_params(value),
      scatterer->flags.use_u_aniso() && scatterer->flags.grad_u_aniso()));
  }

  independent_small_vector_parameter<6> const &independent_params() {
    return *(independent_small_vector_parameter<6> *)argument(0);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

private:
  adp_constraints_t adp_constraints;
};



}}}

#endif // GUARD
