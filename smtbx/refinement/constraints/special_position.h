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
class special_position_site : public crystallographic_parameter
{
public:
  special_position_site(sgtbx::site_symmetry const &site_symmetry,
                        scatterer_type *scatterer)
    : crystallographic_parameter(1),
      site_constraints(site_symmetry.site_constraints()),
      scatterer(scatterer)
  {
    site = site_symmetry.special_op()*scatterer->site;
    set_arguments(new independent_small_vector_parameter<3>(
      site_constraints.independent_params(site),
      scatterer->flags.grad_site()));
  }

  independent_small_vector_parameter<3> const &independent_params() {
    return *(independent_small_vector_parameter<3> *)argument(0);
  }

  virtual std::size_t size() const;

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual void store(uctbx::unit_cell const &unit_cell) const;

private:
  frac_t site;
  sgtbx::site_constraints<double> site_constraints;
  scatterer_type * scatterer;
};


/// Anisotropic displacement constrained by the symmetry of a special position
/** Parameter components are those of the tensor in Cartesian coordinates
 */
class special_position_cartesian_adp : public crystallographic_parameter
{
public:
  typedef sgtbx::tensor_rank_2::cartesian_constraints<double>
          adp_constraints_t;

  special_position_cartesian_adp(sgtbx::site_symmetry const &site_symmetry,
                                 uctbx::unit_cell const &unit_cell,
                                 scatterer_type *scatterer)
    : crystallographic_parameter(1),
      adp_constraints(site_symmetry.cartesian_adp_constraints(unit_cell)),
      scatterer(scatterer)
  {
    tensor_rank_2_t u_star = site_symmetry.average_u_star(scatterer->u_star);
    u_cart = adptbx::u_star_as_u_cart(unit_cell, u_star);
    set_arguments(new independent_small_vector_parameter<6>(
      adp_constraints.independent_params(u_cart),
      scatterer->flags.use_u_aniso() && scatterer->flags.grad_u_aniso()));
  }

    independent_small_vector_parameter<6> const &independent_params() {
    return *(independent_small_vector_parameter<6> *)argument(0);
  }

  virtual std::size_t size() const;

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual void store(uctbx::unit_cell const &unit_cell) const;

private:
  tensor_rank_2_t u_cart;
  adp_constraints_t adp_constraints;
  scatterer_type * scatterer;
};



}}}

#endif // GUARD
