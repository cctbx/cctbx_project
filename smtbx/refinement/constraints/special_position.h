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
class special_position_site : public site_parameter
{
public:
  special_position_site(sgtbx::site_symmetry const &site_symmetry,
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
/** Parameter components are those of the tensor in Cartesian coordinates
 */
class special_position_cartesian_adp : public cartesian_adp
{
public:
  typedef sgtbx::tensor_rank_2::cartesian_constraints<double>
          adp_constraints_t;

  special_position_cartesian_adp(sgtbx::site_symmetry const &site_symmetry,
                                 uctbx::unit_cell const &unit_cell,
                                 scatterer_type *scatterer)
    : cartesian_adp(scatterer, 1),
      adp_constraints(site_symmetry.cartesian_adp_constraints(unit_cell))
  {
    tensor_rank_2_t u_star = site_symmetry.average_u_star(scatterer->u_star);
    value = adptbx::u_star_as_u_cart(unit_cell, u_star);
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
