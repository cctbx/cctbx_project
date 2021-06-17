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
class special_position_site_parameter : public asu_site_parameter
{
public:
  special_position_site_parameter(sgtbx::site_symmetry_ops const &site_symmetry,
                                  scatterer_type *scatterer)
    : parameter(1),
      single_asu_scatterer_parameter(scatterer),
      site_constraints(site_symmetry.site_constraints())
  {
    value = site_symmetry.special_op()*scatterer->site;
    set_arguments(new independent_small_vector_parameter<3>(
      site_constraints.independent_params(value),
      scatterer->flags.grad_site()));
  }

  independent_small_vector_parameter<3> const &independent_params() {
    return *dynamic_cast<independent_small_vector_parameter<3> *>(argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
  const sgtbx::site_constraints<double> &get_site_constraints() const {
    return site_constraints;
  }
private:
  sgtbx::site_constraints<double> site_constraints;
};


/// Anisotropic displacement constrained by the symmetry of a special position
/** Parameter components are those of the tensor in fractional coordinates
 */
class special_position_u_star_parameter : public asu_u_star_parameter
{
public:
  typedef sgtbx::tensor_rank_2::constraints<double>
          adp_constraints_t;

  special_position_u_star_parameter(sgtbx::site_symmetry_ops const &site_symmetry,
                                    scatterer_type *scatterer)
    : parameter(1),
      single_asu_scatterer_parameter(scatterer),
      adp_constraints(site_symmetry.adp_constraints())
  {
    value = site_symmetry.average_u_star(scatterer->u_star);
    set_arguments(new independent_small_vector_parameter<6>(
      adp_constraints.independent_params(value),
      scatterer->flags.use_u_aniso() && scatterer->flags.grad_u_aniso()));
  }

  independent_small_vector_parameter<6> const &independent_params() {
    return *dynamic_cast<independent_small_vector_parameter<6> *>(argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

private:
  adp_constraints_t adp_constraints;
};

/// Anisotropic displacement constrained by the symmetry of a special position
/** Parameter components are those of the tensor in fractional coordinates
 */
class special_position_anharmonic_adp_parameter : public asu_anharmonic_adp_parameter
{
public:
  typedef sgtbx::site_symmetry_ops::tensor_rank_3_constraints_t
    tensor_r3_constraints_t;
  typedef sgtbx::site_symmetry_ops::tensor_rank_4_constraints_t
    tensor_r4_constraints_t;

  special_position_anharmonic_adp_parameter(
    sgtbx::site_symmetry_ops const &site_symmetry,
    scatterer_type *scatterer)
    : parameter(1),
    single_asu_scatterer_parameter(scatterer),
    tensor_r3_constraints(site_symmetry.tensor_rank_3_constraints()),
    tensor_r4_constraints(site_symmetry.tensor_rank_4_constraints())
  {
    for (size_t i = 0; i < 10; i++) {
      value[i] = scatterer->anharmonic_adp->C[i];
    }
    for (size_t i = 0; i < 15; i++) {
      value[i + 10] = scatterer->anharmonic_adp->D[i];
    }
    af::shared<double> params = tensor_r3_constraints.independent_params(
      scatterer->anharmonic_adp->C.data());
    c_count = params.size();
    af::shared<double> ip = tensor_r4_constraints.independent_params(
      scatterer->anharmonic_adp->D.data());
    for (size_t i = 0; i < ip.size(); i++) {
      params.push_back(ip[i]);
    }
    set_arguments(new independent_vector_parameter(
      params,
      scatterer->flags.use_u_aniso() && scatterer->flags.grad_u_aniso()));
  }

  independent_vector_parameter const &independent_params() {
    return *dynamic_cast<independent_vector_parameter *>(argument(0));
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
    sparse_matrix_type *jacobian_transpose);

private:
  size_t c_count;
  tensor_r3_constraints_t tensor_r3_constraints;
  tensor_r4_constraints_t tensor_r4_constraints;
};


}}}

#endif // GUARD
